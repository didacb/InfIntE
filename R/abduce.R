#' Abduce
#'
#' Uses PyGol abduction procedure to get the compression values of the different effects
#'
#' @param bottom_clauses list produced by the \code{\link{get_bottom_clause}} function.
#' @param hypothesis character vector where each element relates the different observations produced by the get_observation functions.
#' @param head_clauses character vector with new head clauses that override the ones present in the bottom clause list.
#'
#' @return a data frame where each row contains a pair of OTU, the effect on the abundance and the compression supporting the effect.
#' @export
#'
#' @examples
#' abduce(bottom_clause, hypothesis)
#'
abduce <- function(bottom_clauses, hypothesis, head_clauses = NULL, abducible = NULL) {
  if (is.null(head_clauses)) {
    head_clauses <- bottom_clauses$head
  }

  # Define abducibles
  if (is.null(abducible)){
    abducible <- c("predation", "no_predation")
  }
  
  # Source python function
  load_PyGol()

  # Execute abduction using InfIntE
  coverage <- abduction(bottom_clauses$clauses, abducible,
    positive_example_list = head_clauses,
    constant_set = bottom_clauses$const,
    meta_rule = hypothesis, metric = "predictive_power"
  )

  # Extract compression values from python object
  compressions <- vapply(names(coverage), function(x) {
    coverage[[x]]
  }, FUN.VALUE = numeric(1))

  # Join compression with interaction names
  coverage <- cbind(names(coverage), unname(compressions))

  # Format the output into a table
  abduced <- formatOutput(coverage)

  return(abduced)
}
