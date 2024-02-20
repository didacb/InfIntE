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
get_temporal_observations <- function(otu, t0, t1, exclusion = TRUE, zero = FALSE) {
  samps <- colnames(t0)
  spec <- otu
  
  # For each comparison
  abundance_clauses <- vapply(samps, function(x) {
    # Subset pair of samples
    pair_samps <- c(t0[otu, x],t1[otu, x])
    abundance_unit <- as.character(NA)
  
    if (all(pair_samps != 0)) {
      if(zero){
        abundance_unit <-  paste0("abundance(", x, ",", spec, ",zero).") 
      }
      # If t0 is smaller than t1
      if ((log(pair_samps[2]) - log(pair_samps[1])) / log(pair_samps[2]) > 0.1 & pair_samps[1] < pair_samps[2]) {
        abundance_unit <- paste0("abundance(", x, ",", spec, ",up).")
      }
      # If log samp1 is bigger than samp2
      if  ((log(pair_samps[1]) - log(pair_samps[2])) / log(pair_samps[2]) > 0.1 & pair_samps[1] > pair_samps[2]) {
        abundance_unit <- paste0("abundance(", x, ",", spec, ",down).")
      }
    } else {
      if (exclusion) {
        if (pair_samps[1] == 0 & pair_samps[2] > 0) {
          abundance_unit <- paste0(
            "abundance(", x, ",", spec, ",up)."
          )
        }
        if (pair_samps[1] > 0 & pair_samps[2] == 0) {
          abundance_unit <- paste0(
            "abundance(", x, ",", spec, ",down)."
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
