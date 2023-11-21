#' Get presence clauses
#' Generate logical clauses depicting the presence of each OTU in each sample
#' @param otu_data list joining all the OTU information. Obtained using the join_abundances function.
#'
#' @return character vector with the presence clauses describing the presence of the OTUs in the samples.
#' @export
#'
#' @examples
#' get_presence(otu_data)
get_presence <- function(otu_tb) {

  # Iterate by ASV table rows
  presence <- lapply(seq_len(nrow(otu_tb)), function(x) {

    # Check presence in each sample
    yes.no <- ifelse(as.numeric(otu_tb[x, ]) > 0, "yes", "no")

    # Build character vector
    presence <- paste0("presence(", colnames(otu_tb), ",", rownames(otu_tb)[x], ",", unname(yes.no), ").")
    return(presence)
  })
  return(unlist(presence))
}
