

#' InfIntE
#' InfIntE (INFerence of INTeractions using Explainable machine learning). This function contains the full InfIntE pipeline. It uses compositional and/or absolute OTU data to infer ecological interactions using a logical hypothesis of interaction.
#' The InfIne inference process parses the OTU data into logic clauses and performs the abduction of interactions using PyGol. Interactions are selected using pulsar.
#' @param otu_tb data frame with the OTUs counts with samples on columns and OTUs in rows.
#' @param thresh numeric threshold for StARS model selection in pulsar, default is 0.01
#' @param exclusion logical indicating if exclusion cases have to be considered for generating the logic clauses. Default is TRUE.
#' @param nperms  numeric number of subsamples in pulsar, default is 50
#' @param search.depth numeric search depth when producing constructing the bottom clauses, default is 2
#' @param depth numeric vector containing an alternative sequence depth in case the compositional OTU table is a subset of one with more OTUs
#' @param hypothesis character vector where each element relates the different observations produced by the get_observation functions.
#' @param absolute_abundance data frame with the absolute counts with samples on columns and OTUs in rows.
#' @param ncores number of cores used by pulsar
#'
#' @return list with the infinite output. Contains a data frame with the interactions inferred by InfIntE, the pulsar output and the abduction output.
#' @export
#'
#' @examples
#' infinte(otu_tb, hypothesis, thresh = 0.01, exclusion = TRUE, nperms = 50, search.depth = 2, depth = NULL, absolute_abundance = NULL)
infinte <- function(otu_tb, thresh = 0.01, exclusion = TRUE, nperms = 50, search.depth = 2, hypothesis=NULL, depth = NULL, absolute_abundance = NULL, ncores=1) {

  if(is.null(hypothesis)){
    if (exclusion){
      hypothesis<- c("abundance(C1,C2,S1,up):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_up(S2,S1)",
                     "abundance(C1,C2,S1,app):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_up(S2,S1)",
                     "abundance(C1,C2,S1,down):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_down(S2,S1)",
                     "abundance(C1,C2,S1,dis):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_down(S2,S1)")
    }else{
      hypothesis<- c("abundance(C1,C2,S1,up):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_up(S2,S1)",
                     "abundance(C1,C2,S1,down):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_down(S2,S1)")

    }
  }
  # Join absolute and compositional data in a table
  otu_data <- join_abundances(otu_tb, absolute_abundance, depth)

  # All possible pairs of samples
  comparisons <- get_comparsions(length(otu_data$samp_names))


  # Get head logic clauses
  head_clauses <- lapply(rownames(otu_data$otu_tb), function(otu) {
    pos <- which(rownames(otu_data$otu_tb) == otu)
    abundances <- do.call(
      what = otu_data$abundance_function[pos],
      args = list(
        "otu_abundance" = otu_data$otu_tb[pos, , drop = FALSE],
        "comparisons" = comparisons, "depth" = otu_data$depth, "exclusion" = exclusion
      )
    )
    return(abundances)
  })
  head_clauses <- unlist(head_clauses)

  # Get Body logic clauses
  body_clauses <- get_presence(otu_data)

  # Produce bottom clause
  bottom_clauses <- get_bottom_clause(otu_data = otu_data,
                                      head_clauses = head_clauses, body_clauses = body_clauses)

  # Abduce effects
  abduced_effects <- abduce(bottom = bottom_clauses, hypothesis = hypothesis, head_clauses = head_clauses)

  # Get I values
  abduced_effects <- get_I_values(abduced_effects)

  # Length observations
  mx <- length(bottom_clauses$head)

  # Lambda distribution
  lambda <- pulsar::getLamPath(max = mx, min = 0, len = nperms)

  # Pulsar execution
  pulsar_output <- pulsar::pulsar(t(otu_data$otu_tb),
    fun = pulsar_infinte,
    fargs = list(lambda = lambda, bottom_clauses = bottom_clauses,
                 hypothesis = hypothesis, exclusion = exclusion, otu_data = otu_data),
    rep.num = nperms, lb.stars = FALSE, ub.stars = FALSE, thresh = thresh, ncores = ncores
  )

  # Format output to dataframe
  fitted_model <- pulsar::refit(pulsar_output, criterion = "stars")
  interactions <- data.frame(igraph::get.edgelist(igraph::graph_from_adjacency_matrix(fitted_model$refit$stars)))

  # Take values from abduced effects dataframe
  interactions <- abduced_effects[paste0(abduced_effects[, 1], abduced_effects[, 2])
                                  %in% paste0(interactions[, 1], interactions[, 2]), ]

  # Classify and give back original names
  interactions <- classify_interactions(interactions)
  interactions <- return_names(interactions, otu_data$otu_names)

  # Prepare output object
  infinte_output <- list(selected_interactions = interactions,
                         pulsar_result = pulsar_output,
                         abduced_table = return_names(abduced_effects, otu_data$otu_names))

  return(infinte_output)
}
