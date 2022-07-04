
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