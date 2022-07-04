# Function to apply chisq test to all replicate comparison for one species and obtain the abundance progol input
do_chi_abundance <- function(comparison, otu_abundance, depth, exclusion) {

  # Detect the pair of values and depth to compare
  pair <- otu_abundance[comparison]
  samps<- colnames(otu_abundance)
  spec<- rownames(otu_abundance)
  pair_depth <- depth[comparison]

  # If both pair values no zero and the samples are not equal to depth
  if (any(pair != pair_depth)) {
    if (all(pair != 0)) {

      # Make table and add margins
      contingency_tab <- as.matrix(rbind(pair, pair_depth - pair))
      contingency_tab <- addmargins(A = contingency_tab, margin = seq_along(dim(contingency_tab)))

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
        abundance_unit <- paste0("abundance(", samps[comparison[1]], ",", samps[comparison[2]], ",", spec, ",", up.down, ").")
      } else {
        abundance_unit <- as.character(NA)
      }
    } else {
      if (any(pair != 0) & exclusion) {
        # Make contingency_table and add margins
        contingency_tab <- as.matrix(rbind(pair, pair_depth - pair))
        contingency_tab <- addmargins(A = contingency_tab, margin = seq_along(dim(contingency_tab)))

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
          abundance_unit <- paste0("abundance(", samps[comparison[1]], ",", samps[comparison[2]], ",",  spec, ",", up.down, ").")
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


get_compositional_observations<- function(otu_abundance, comparisons, depth, exclusion) {

  # Run pairwise
  abundance_clauses <- vapply(comparisons, do_chi_abundance, otu_abundance, depth, exclusion, FUN.VALUE = character(1))
  # Delete NAs
  abundance_clauses <- abundance_clauses[!is.na(abundance_clauses)]

  return(abundance_clauses)
}

get_absolute_observations <- function(otu_abundance, comparisons, depth, exclusion) {

  samps<- colnames(otu_abundance)
  spec<- rownames(otu_abundance)

  log_abundance <- log(otu_abundance)
  log_abundance[log_abundance == -Inf] <- 0

  # For each comparison
  abundance_clauses <- vapply(comparisons, function(x) {
    # Subset pair of samples
    pair_samps <- log_abundance[,x]
    abundance_unit <- as.character(NA)
    if (all(pair_samps != 0)) {
      # If log samp1 is 0.05 smaller than samp2
      if (abs(pair_samps[1] - pair_samps[2]) > 0.05 & pair_samps[1] < pair_samps[2]) {
        abundance_unit <- paste0("abundance(", samps[x[1]], ",", samps[x[2]], ",", spec, ",up).")
      }
      # If log samp1 is 0.5 bigger than samp2
      if (abs(pair_samps[1] - pair_samps[2]) > 0.05 & pair_samps[1] > pair_samps[2]) {
        abundance_unit <- paste0("abundance(", samps[x[1]], ",", samps[x[2]], ",", spec, ",down).")
      }
    } else {
      if (exclusion) {
        if (pair_samps[1] == 0 & pair_samps[2] > 0) {
          abundance_unit <- paste0("abundance(", samps[x[1]], ",", samps[x[2]], ",", spec, ",app).")
        }
        if (pair_samps[1] > 0 & pair_samps[2] == 0) {
          abundance_unit <- paste0("abundance(", samps[x[1]], ",", samps[x[2]], ",", spec, ",dis).")
        }
      }
    }
    return(abundance_unit)
  }, FUN.VALUE = character(1))
  # Delete NAs
  abundance_clauses <- abundance_clauses[!is.na(abundance_clauses)]
  return(abundance_clauses)
}

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
