#' Get absolute observations
#' Generate logical clauses depicting the absolute abundance change between two samples for each OTU.
#' @param otu_abundance list joining all the OTU information. Obtained using the join_abundances function.
#' @param comparisons list with the position of the pairs of samples to compare.
#' @param depth numeric vector with the sequencing depth
#' @param exclusion logical indicating if exclusion cases have to be considered for generating the logic clauses. Default is FALSE.
#'
#' @return character vector with the abundance clauses describing the absolute abundance change between samples.
#' @export
#'
#' @examples
#' get_absolute_observations(otu_abundance, comparisons)
get_absolute_observations <- function(otu_abundance, comparisons, depth, exclusion = FALSE) {
  samps <- colnames(otu_abundance)
  spec <- rownames(otu_abundance)

  log_abundance <- log(otu_abundance)
  log_abundance[log_abundance == -Inf] <- 0

  # For each comparison
  abundance_clauses <- vapply(comparisons, function(x) {
    # Subset pair of samples
    pair_samps <- log_abundance[, x]
    abundance_unit <- as.character(NA)
    if (all(pair_samps != 0)) {
      # If log samp1 is smaller than samp2
      if (abs(pair_samps[1] - pair_samps[2]) > 0.5 & pair_samps[1] < pair_samps[2]) {
        abundance_unit <- paste0(
          "abundance(", samps[x[1]], ",",
          samps[x[2]], ",", spec, ",up)."
        )
      }
      # If log samp1 is 0.5 bigger than samp2
      if (abs(pair_samps[1] - pair_samps[2]) > 0.5 & pair_samps[1] > pair_samps[2]) {
        abundance_unit <- paste0(
          "abundance(", samps[x[1]], ",",
          samps[x[2]], ",", spec, ",down)."
        )
      }
    } else {
      if (exclusion) {
        if (pair_samps[1] == 0 & pair_samps[2] > 0) {
          abundance_unit <- paste0(
            "abundance(", samps[x[1]], ",",
            samps[x[2]], ",", spec, ",app)."
          )
        }
        if (pair_samps[1] > 0 & pair_samps[2] == 0) {
          abundance_unit <- paste0(
            "abundance(", samps[x[1]], ",",
            samps[x[2]], ",", spec, ",dis)."
          )
        }
      }
    }
    return(abundance_unit)
  }, FUN.VALUE = character(1))
  # Delete NAs
  abundance_clauses <- abundance_clauses[!is.na(abundance_clauses)]
  return(abundance_clauses)
}
