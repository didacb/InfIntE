

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
  body_clauses<- get_presence(otu_data)

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
