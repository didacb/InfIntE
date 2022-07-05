#' Join abundances
#' Joins in a single data frame the compositional and absolute OTU counts and returns a list with the new data frame and all the information of samples and OTUs
#' @param otu_tb data frame with the OTUs counts with samples on columns and OTUs in rows.
#' @param absolute_abundance data frame with the absolute counts with samples on columns and OTUs in rows.
#' @param depth numeric vector alternative sequencing depth in case OTU table is a subset.
#'
#' @return list joining all the OTU information including an OTU table with absolute and/or compositional sequence counts.
#' @export
#'
#' @examples
#' join_abundances(otu_tb, absolute_abundance, depth = NULL)
join_abundances <- function(otu_tb, absolute_abundance, depth = NULL) {
  change_position <- nrow(otu_tb)

  if (!is.null(absolute_abundance) & !is.null(otu_tb)) {

    # Join tables
    colnames(absolute_abundance) <- colnames(otu_tb)
    otu_tb <- rbind(otu_tb, absolute_abundance)

    end_position <- nrow(otu_tb) - change_position
  }

  if (is.null(otu_tb)) {
    change_position <- 0
    otu_tb <- absolute_abundance
    end_postion <- nrow(otu_tb)
  }

  # Save given OTU and sample names
  samp_names <- colnames(otu_tb)
  otu_names <- rownames(otu_tb)

  # Set simplified names
  rownames(otu_tb) <- paste0("s", seq_len(nrow(otu_tb)))
  colnames(otu_tb) <- paste0("c", seq_len(ncol(otu_tb)))

  # Check table
  otu_tb <- check_OTU_table(otu_tb)

  if (is.null(depth) & !is.null(otu_tb)) {
    depth <- colSums(otu_tb[1:change_position, ])
  }

  abundance_function <- c(
    rep("get_compositional_observations", change_position),
    rep("get_absolute_observations", end_position)
  )

  names(abundance_function) <- rownames(otu_tb)
  names(depth) <- colnames(otu_tb)

  otu_data <- list(
    "otu_tb" = otu_tb, "otu_names" = otu_names, "samp_names" = samp_names,
    "depth" = depth, "abundance_function" = abundance_function
  )
  return(otu_data)
}
