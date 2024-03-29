
#' Generate bottom clause
#' Uses PyGol to generate the bottom clauses from the head and body observations of the hypothesis of interaction.
#' @param otu_data list joining all the OTU information. Obtained using the join_abundances function.
#' @param head_clauses character vector with the clauses in the head position of the hypothesis of interaction.
#' @param body_clauses character vector  with the clauses in the body position of the hypothesis of interaction.
#' @param search.depth numeric selecting how deep has to search PyGol to generate the abduction clauses.
#'
#' @return list including a python dictionary with the bottom clauses together with all the information required for the abduction.
#' @export
#'
#' @examples
#' get_bottom_clause(otu_data, head_clauses, body_clauses)
get_bottom_clause <- function(otu_data, head_clauses, body_clauses, search.depth = 2) {

  # Get constants
  const <- vapply(c(head_clauses, body_clauses), function(clau) {
    last_particle <- unlist(strsplit(clau, ","))
    last_particle <- last_particle[length(last_particle)]
    last_particle <- gsub(").", "", last_particle)
  }, FUN.VALUE = character(1))
  const <- unique(unname(const))

  const <- c(const, rownames(otu_data$otu_tb))
  body1_clauses <- gsub("presence", "presence1", body_clauses)

  # Source InfIntE
  load_PyGol()

  # Generate bottom clauses
  P <- generate_bottom_clause(c(body_clauses, body1_clauses), const, head_clauses,
    NULL,
    container = "memory", depth = search.depth
  )

  # Format bottom clause
  P <- reticulate::dict(P, convert = TRUE)

  # Create and object with all the elements necessaty for abduction
  bottom <- list(P, head_clauses, body_clauses, const, otu_data$otu_tb)
  names(bottom) <- c("clauses", "head", "body", "const", "otu_tb")
  return(bottom)
}
