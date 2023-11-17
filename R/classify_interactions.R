
#' Classify interactions
#' Automatic classification of significant effects to obtain the types of interacions
#' @param effect_table a data frame, containing the significant effects where each row contains a pair of OTU, the effect on the abundance and the I value supporting the effect.
#'
#' @return a data frame, similar to effect_table, where effect_up and effect_down has been changed by the types of interactions proposed by Derocles et al. (2018)
#' @export
#'
#' @examples
#' classify_interactions(effect_table)
classify_interactions <- function(effect_table) {
  right <- paste0(effect_table[, 1], effect_table[, 2])
  oposite <- paste0(effect_table[, 2], effect_table[, 1])

  interactions <- vapply(right, function(x) {
    if (any(oposite == x)) {
      if (which(oposite == x) < which(right == x)) {
        int <- paste0(
          effect_table[which(right == x), 3], "/",
          effect_table[which(oposite == x), 3]
        )
      } else {
        int <- "take_out"
      }
    } else {
      int <- effect_table[which(right == x), 3]
    }
  }, FUN.VALUE = character(1))

  interactions <- gsub("predation/predation", "competition", interactions)

  effect_table[, 3] <- interactions
  effect_table <- effect_table[effect_table[, 3] != "take_out", ]

  return(effect_table)
}
