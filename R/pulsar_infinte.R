
#' InfInte function for pulsar
#' Function to be used by pulsar to perform the abduction from a subset of the OTU table and obtain a list of adjacency matrix along the lambda path
#' @param sub_otu_tb subset of the OTU table provided by pulsar
#' @param lambda numeric vector containing the lambda path obtained with the \code{\link{getLamPath}} function.
#' @param bottom_clauses list including a python dictionary with the bottom clauses together with all the information required for the abduction.
#' @param hypothesis character vector where each element relates the different observations produced by the get_observation functions.
#' @param exclusion logical indicating if exclusion cases have to be considered for generating the logic clauses. Default is FALSE.
#' @param otu_data list joining all the OTU information. Obtained using the join_abundances function.
#'
#' @return a list of p∗p adjacency matrices, one for each value of lambda.
#' @export
#'
#' @examples
#' pulsar_infinte(sub_otu_tb, lambda, bottom_clauses, hypothesis, exclusion, otu_data)
pulsar_infinte <- function(sub_otu_tb, lambda, bottom_clauses, hypothesis, exclusion, otu_data) {
  sub_otu_tb <- t(sub_otu_tb)
  # Get head clauses from sub sampled table

  comparisons <- get_comparsions(ncol(sub_otu_tb))
  sub_depth <- otu_data$depth[colnames(sub_otu_tb)]

  # Head
  head_clauses <- lapply(rownames(sub_otu_tb), function(otu) {
    pos <- which(rownames(sub_otu_tb) == otu)
    abundances <- do.call(
      what = otu_data$abundance_function[pos],
      args = list(
        "otu_abundance" = sub_otu_tb[pos, , drop = FALSE],
        "comparisons" = comparisons, "depth" = sub_depth, "exclusion" = exclusion
      )
    )
    return(abundances)
  })
  head_clauses <- unlist(head_clauses)

  # Abduce
  abduced.table <- abduce(bottom_clauses =  bottom_clauses, hypothesis = hypothesis, head_clauses = head_clauses)

  # Get I values
  abduced.table <- get_I_values(abduced.table)

  snames <- rownames(otu_data$otu_tb)
  
  # Add non interacting asvs
  noi.asvs <- snames[!snames %in% unique(c(abduced.table$sp1, abduced.table$sp2))]
  noi.tab <- data.frame(noi.asvs, noi.asvs, rep("no_effect",
                                        length(noi.asvs)), rep(0, length(noi.asvs)))
  colnames(noi.tab) <- colnames(abduced.table)
  abduced.table <- rbind(abduced.table, noi.tab)
  
  # Construct adjacency matrix
  g <- igraph::graph_from_data_frame(abduced.table[, c(1, 2, 4)])
  ad <- igraph::as_adj(g, attr = "comp")
  #ad[upper.tri(ad)]<- ad[lower.tri(ad)]
  
  # For each lambda
  pt <- lapply(lambda, function(lam) {
    # Compression bigger than lambda
    tmp <- ad > lam
    diag(tmp) <- FALSE
    
    tmp<- as(tmp, "lMatrix")
    return(tmp)
  })
  # Return the path over lambda
  return(list(path = pt))
}
