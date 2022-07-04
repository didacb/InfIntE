
loadInfIntE <- function() {
  # package.location<-  system.file(package = "InfIntE-nets")
  package.location <- "/home/didac/Desktop/InfIntE/str"
  current_wd <- getwd()
  setwd(package.location)

  # reticulate::use_python(system2("which", "python3", stdout = TRUE))

  reticulate::use_python("/usr/bin/python3")

  reticulate::source_python(file.path(package.location, "pygolm_abduce.py"), envir = globalenv())
  setwd(current_wd)
  on.exit(setwd( current_wd))
}

abduce <- function(bottom_clauses, hypothesis, head_clauses = NULL) {
  if (is.null(head_clauses)) {
    head_clauses <- bottom_clauses$head
  }

  # Define abducibles
  abducible <- c("effect_up", "effect_down")

  # Source python function
  loadInfIntE()

  # Execute abduction using InfIntE
  coverage <- abduction(bottom_clauses$clauses, abducible, positive_example_list = head_clauses,
                        constant_set = bottom_clauses$const, meta_rule = hypothesis, metric = "predictive_power")

  # Extract compression values from python object
  compressions <- vapply(names(coverage), function(x) {
    coverage[[x]]
  }, FUN.VALUE = numeric(1))

  # Join compression with interaction names
  coverage <- cbind(names(coverage), unname(compressions))

  # Format the output into a table
  coverage <- formatOutput(coverage)

  return(coverage)
}

get_bottom_clause <- function(otu_data, head_clauses, body_clauses, search.depth = 2, cores = 1) {

  # Get constants
  const <- vapply(c(head_clauses, body_clauses), function(clau){
              last_particle<- unlist(strsplit(clau, ","))
              last_particle<- last_particle[length(last_particle)]
              last_particle<- gsub(").", "", last_particle)
  }, FUN.VALUE = character(1))
  const<- unique(unname(const))

  const <- c(const, rownames(otu_data$otu_tb))
  body1_clauses <- gsub("presence", "presence1", body_clauses)

  # Source InfIntE
  loadInfIntE()

  # Generate bottom clauses
  P <- generate_bottom_clause(c(body_clauses, body1_clauses), const, head_clauses,
                              NULL, container = "memory", depth = search.depth)

  # Format bottom clause
  P <- reticulate::dict(P, convert = TRUE)

  # Create and object with all the elements necessaty for abduction
  bottom <- list(P, head_clauses, body_clauses,  const, otu_data$otu_tb)
  names(bottom) <- c("clauses", "head", "body", "const", "otu_tb")
  return(bottom)
}

pulsar_infinte <- function(sub_otu_tb, lambda, bottom_clauses, hypothesis, exclusion, otu_data) {

  sub_otu_tb<- t(sub_otu_tb)
  # Get head clauses from sub sampled table

  comparisons<- get_comparsions(ncol(sub_otu_tb))
  sub_depth<- otu_data$depth[colnames(sub_otu_tb)]

  #Head
  head_clauses<- lapply(rownames(sub_otu_tb), function(otu){
    pos<- which(rownames(sub_otu_tb) == otu)
    abundances<- do.call(what = otu_data$abundance_function[pos],
                         args = list("otu_abundance"=sub_otu_tb[pos,,drop=FALSE],
                                      "comparisons"=comparisons, "depth"=sub_depth, "exclusion"=exclusion))
    return(abundances)
  })
  head_clauses<- unlist(head_clauses)

  #Abduce
  abduced <- abduce(bottom = bottom_clauses, hypothesis = hypothesis)

  #Get I values
  abduced.table<- get_I_values(abduced)

  snames<- rownames(otu_data$otu_tb)
  # Add non interacting asvs
  noi.asvs <- snames[!snames %in% unique(c(abduced.table$sp1, abduced.table$sp2))]
  noi.tab <- data.frame(noi.asvs, noi.asvs, rep("no_effect", length(noi.asvs)), rep(0, length(noi.asvs)))
  colnames(noi.tab) <- colnames(abduced.table)
  abduced.table <- rbind(abduced.table, noi.tab)

  # Construct adjacency matrix
  g <- igraph::graph_from_data_frame(abduced.table[, c(1, 2, 4)])
  ad <- igraph::get.adjacency(g, attr = "comp", sparse = TRUE)

  # For each lambda
  pt <- lapply(lambda, function(lam) {
    # Compression bigger than lambda
    tmp <- ad > lam
    ga <- igraph::graph_from_adjacency_matrix(tmp)
    tmp <- igraph::get.adjacency(ga, sparse = TRUE)

    diag(tmp) <- FALSE
    return(tmp)
  })
  # Return the path over lambda
  return(list(path = pt))
}



infinte <- function(otu_tb, hypothesis, thresh = 0.01, exclusion = FALSE, nperms = 50, search.depth = 2, qpcr = NULL, depth = NULL, absolute_abundance=NULL) {

  # Join absolute and compositional data in a table
  otu_data<- join_abundances(otu_tb, absolute_abundance, depth)

  # All possible pairs of samples
  comparisons<- get_comparsions(length(otu_data$samp_names))


  #Get head logic clauses
  head_clauses<- lapply(rownames(otu_data$otu_tb), function(otu){
    pos<- which(rownames(otu_data$otu_tb) == otu)
    abundances<- do.call(what = otu_data$abundance_function[pos],
                         args = list("otu_abundance"=otu_data$otu_tb[pos,,drop=FALSE],
                         "comparisons"=comparisons, "depth"=otu_data$depth, "exclusion"=exclusion))
    return(abundances)
  })
  head_clauses<- unlist(head_clauses)

  #Get Body logic clauses
  body_clauses<- get_presence(otu_data$otu_tb)

  # Produce bottom clause
  bottom_clauses <- get_bottom_clause(otu_data = otu_data, head_clauses = head_clauses, body_clauses = body_clauses)

  # Abduce effects
  abduced_effects <- abduce(bottom = bottom_clauses, hypothesis = hypothesis)

  #Get I values
  abduced_effects<- get_I_values(abduced_effects)

  # Length observations
  mx <- length(bottom_clauses$head)

  # Lambda distribution
  lambda <- pulsar::getLamPath(max = mx, min = 0, 50, FALSE)

  # Pulsar execution
  pulsar_output <- pulsar::pulsar(t(otu_data$otu_tb), fun = pulsar_infinte,
                       fargs = list(lambda = lambda, bottom_clauses = bottom_clauses, hypothesis = hypothesis, exclusion = exclusion, otu_data=otu_data),
                       rep.num = nperms, lb.stars = TRUE, ub.stars = TRUE, thresh = thresh)

  # Format output to dataframe
  fitted_model <- pulsar::refit(pulsar_output, criterion = "stars")
  interactions <- data.frame(igraph::get.edgelist(igraph::graph_from_adjacency_matrix(fitted_model$refit$stars)))

  # Take values from abduced effects dataframe
  interactions <- abduced_effects[paste0(abduced_effects[, 1], abduced_effects[, 2]) %in% paste0(interactions[, 1], interactions[, 2]), ]

  # Classify ad give bac original names
  interactions <- classifyInteraction(interactions)
  interactions<- return_names(interactions, otu_data$otu_names)

  # Prepare output object
  infinte_output <- list(selected_interactions = interactions, pulsar_result = pulsar_output, abduced_table = return_names(abduced_effects, otu_data$otu_names))

  return(infinte_output)
}

##############################################################################
############################# Classify interactions #########################
#############################################################################
classifyInteraction <- function(effect_table) {
  right <- paste0(effect_table[, 1], effect_table[, 2])
  oposite <- paste0(effect_table[, 2], effect_table[, 1])

  interactions <- vapply(right, function(x) {
    if (any(oposite == x)) {
      if (which(oposite == x) < which(right == x)) {
        int <- paste0(effect_table[which(right == x), 3], "/", effect_table[which(oposite == x), 3])
      } else {
        int <- "take_out"
      }
    } else {
      int <- effect_table[which(right == x), 3]
    }
  }, FUN.VALUE = character(1))

  interactions <- gsub("effect_down/effect_down", "competition", interactions)
  interactions <- gsub("effect_up/effect_up", "mutualism", interactions)
  interactions <- gsub("effect_up/effect_down", "predation", interactions)
  interactions <- gsub("effect_down/effect_up", "take_out", interactions)
  interactions <- gsub("^effect_down$", "amensalism", interactions)
  interactions <- gsub("^effect_up$", "commensalism", interactions)

  effect_table[, 3] <- interactions
  effect_table <- effect_table[effect_table[, 3] != "take_out", ]

  return(effect_table)
}

