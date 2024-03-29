#' Perform Chi Square test
#' Generate a logical clause depicting the compositional abundance change between two samples.
#' @param comparison numeric vector with the position, in the OTU table , of the samples to compare
#' @param otu_abundance list joining all the OTU information. Obtained using the join_abundances function.
#' @param depth numeric vector with the sequencing depth
#' @param exclusion logical indicating if exclusion cases have to be considered for generating the logic clauses. Default is FALSE.
#'
#' @return character with the abundance clause describing the compositional abundance change between samples.
#' @export
#'
#' @examples
#' do_chi_abundance(comparison, otu_abundance, depth)
do_chi_abundance <- function(comparison, otu_abundance, depth, exclusion = FALSE) {

  # Detect the pair of values and depth to compare
  pair <- otu_abundance[comparison]
  samps <- colnames(otu_abundance)
  spec <- rownames(otu_abundance)
  pair_depth <- depth[comparison]

  # If both pair values no zero and the samples are not equal to depth
  if (any(pair != pair_depth)) {
    if (all(pair != 0)) {

      # Make table and add margins
      contingency_tab <- as.matrix(rbind(pair, pair_depth - pair))
      contingency_tab <- addmargins(
        A = contingency_tab,
        margin = seq_along(dim(contingency_tab))
      )

      # Do the cishq test
      chi_test <- suppressWarnings(chisq.test(contingency_tab))
      if (chi_test[3] < 0.05) {
        # If it is significant escontingency_tablish the up and down
        if (pair[1] / pair_depth[1] < pair[2] / pair_depth[2]) {
          up.down <- "up"
        } else {
          up.down <- "down"
        }
        # Build logical abundance input
        abundance_unit <- paste0(
          "abundance(", samps[comparison[1]], ",",
          samps[comparison[2]], ",", spec, ",", up.down, ")."
        )
      } else {
        abundance_unit <- as.character(NA)
      }
    } else {
      if (any(pair != 0) & exclusion) {
        # Make contingency_table and add margins
        contingency_tab <- as.matrix(rbind(pair, pair_depth - pair))
        contingency_tab <- addmargins(
          A = contingency_tab,
          margin = seq_along(dim(contingency_tab))
        )

        # Do the chisq test
        chi_test <- suppressWarnings(chisq.test(contingency_tab))
        if (chi_test[3] < 0.05) {
          # If it is significant escontingency_tablish the up and down
          if (pair[1] / pair_depth[1] < pair[2] / pair_depth[2]) {
            up.down <- "app"
          } else {
            up.down <- "dis"
          }
          # Build logic abundance input
          abundance_unit <- paste0(
            "abundance(", samps[comparison[1]], ",",
            samps[comparison[2]], ",", spec, ",", up.down, ")."
          )
        } else {
          abundance_unit <- as.character(NA)
        }
      } else {
        abundance_unit <- as.character(NA)
      }
    }
  } else {
    abundance_unit <- as.character(NA)
  }
  return(abundance_unit)
}


#' Get compositional observation
#' Generate logical clauses depicting the compositional abundance change between two samples.
#' @param otu_abundance list joining all the OTU information. Obtained using the join_abundances function.
#' @param comparisons list with the position of the pairs of samples to compare.
#' @param depth numeric vector with the sequencing depth
#' @param exclusion logical indicating if exclusion cases have to be considered for generating the logic clauses. Default is FALSE.
#'
#' @return character vector with the abundance clauses describing the compositional abundance change between samples.
#' @export
#'
#' @examples
#' get_compositional_observations(otu_abundance, comparisons, depth)
get_compositional_observations <- function(otu_abundance, comparisons, depth, exclusion = FALSE) {

  # Run pairwise
  abundance_clauses <- vapply(comparisons, do_chi_abundance,
    otu_abundance, depth, exclusion,
    FUN.VALUE = character(1)
  )
  # Delete NAs
  abundance_clauses <- abundance_clauses[!is.na(abundance_clauses)]

  return(abundance_clauses)
}
