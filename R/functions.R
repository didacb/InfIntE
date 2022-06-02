
#' Title
#'
#' @return
#' @export
#'
#' @examples
loadInfIntE <- function() {
  # package.location<-  system.file(package = "InfIntE-nets")
  package.location <- "/home/didac/Desktop/InfIntE/str"
  current_wd <- getwd()
  setwd(package.location)

  # reticulate::use_python(system2("which", "python3", stdout = TRUE))

  reticulate::use_python("/usr/bin/python3")

  reticulate::source_python(file.path(package.location, "pygolm_abduce.py"), envir = globalenv())
  setwd(current_wd)
  on.exit(setwd( current_wd))
}


#' Title
#'
#' @param tb
#'
#' @return
#' @export
#'
#' @examples
check_OTU_table <- function(otu_tab) {
  # Check if all cells are numeric
  otu_tab <- apply(otu_tab, c(1, 2), as.numeric)
  if (any(!is.numeric(otu_tab))) {
    return(stop("Please provide an ASV table with numeric values"))
  }
  # Check if there are NA
  if (any(is.na(otu_tab))) {
    return(stop("Please provide an ASV table without NA"))
  }
  return(otu_tab)
}

#' Title
#'
#' @param tb
#' @param asv.names
#'
#' @return
#' @export
#'
#' @examples
return_names <- function(infered, asv.names) {
  infered[, 1:2] <- apply(infered[, 1:2], c(1, 2), function(x) ifelse(grepl("^s", x), asv.names[x], x))
  return(infered)
}

# Function to apply chisq test to all replicate comparison for one species and obtain the abundance progol input
abundanceChi <- function(comparisons, asv.table, spec, samps, read.depth, exclusion) {
  # Detect the pair of values and depth to compare
  pair <- asv.table[spec, comparisons]
  depth <- read.depth[comparisons]

  # If both pair values no zero and the samples are not equal to depth
  if (any(pair != depth)) {
    if (all(pair != 0)) {

      # Make table and add margins
      tab <- as.matrix(rbind(pair, depth - pair))
      tab <- addmargins(A = tab, margin = seq_along(dim(tab)))

      # Do the cishq test
      tst <- suppressWarnings(chisq.test(tab))
      if (tst[3] < 0.05) {
        # If it is significant establish the up and down
        if (pair[1] / depth[1] < pair[2] / depth[2]) {
          up.down <- "up"
        } else {
          up.down <- "down"
        }
        # Build logical abundance input
        abundance.unit <- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",", up.down, ").")
      } else {
        # Build logical abundance when no significant
        # abundance.unit<- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",zero).")
        abundance.unit <- as.character(NA)
      }
    } else {
      if (any(pair != 0) & exclusion) {
        # Make table and add margins
        tab <- as.matrix(rbind(pair, depth - pair))
        tab <- addmargins(A = tab, margin = seq_along(dim(tab)))

        # Do the chisq test
        tst <- suppressWarnings(chisq.test(tab))
        if (tst[3] < 0.05) {
          # If it is significant establish the up and down
          if (pair[1] / depth[1] < pair[2] / depth[2]) {
            up.down <- "app"
          } else {
            up.down <- "dis"
          }
          # Build logic abundance input
          abundance.unit <- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",", up.down, ").")
        } else {
          # Build logic abundance when no significant
          # abundance.unit<- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",zero).")
          abundance.unit <- as.character(NA)
        }
      } else {
        abundance.unit <- as.character(NA)
      }
    }
  } else {
    abundance.unit <- as.character(NA)
  }
  return(abundance.unit)
}


getObservations <- function(tb, depth = NULL, exclusion, cores = 1) {
  if (is.null(depth)) {
    # Obtain read depth
    read.depth <- colSums(tb)
  } else {
    read.depth <- depth
  }
  pathogen <- NULL
  if (any(rownames(tb) == "pathogen")) {
    pathogen <- tb["pathogen", , drop = FALSE]
    tb <- tb[rownames(tb) != "pathogen", ]
    tb <- data.frame(apply(tb, c(1, 2), as.numeric))
  }

  # All samples to compare
  comparisons <- combn(seq_len(ncol(tb)), 2, NULL, FALSE)

  # Add the opposite direction
  other.direction <- lapply(comparisons, function(x) {
    c(x[2], x[1])
  })
  comparisons <- c(comparisons, other.direction)

  # Run pairwise test function for all species
  abundance <- parallel::mclapply(seq_len(nrow(tb)), function(i) {
    # Run pairwise
    abundance.specie <- vapply(comparisons, abundanceChi, tb, i, colnames(tb), as.numeric(read.depth), exclusion, FUN.VALUE = character(1))
    # Delete NAs
    abundance.specie <- abundance.specie[!is.na(abundance.specie)]
    return(abundance.specie)
    # Number of cores
  }, mc.cores = cores)

  # Unlist and return
  abundance <- unlist(abundance)

  if (!is.null(pathogen)) {
    pathogen <- apply(pathogen, c(1, 2), function(x) {
      x <- as.numeric(x)
      x <- ifelse(is.na(x), 0, x)
      return(x)
    })

    abundance <- c(abundance, getObservationsQpcr(pathogen, exclusion))
  }

  return(abundance)
}

getObservationsQpcr <- function(qpcr.mat, exclusion = FALSE) {

  # Compare all qpcr samples
  comparisons <- combn(seq_len(ncol(qpcr.mat)), 2, NULL, FALSE)
  # Add the opposite direction
  other.direction <- lapply(comparisons, function(x) {
    c(x[2], x[1])
  })
  comparisons <- c(comparisons, other.direction)

  # Sample names
  sn <- colnames(qpcr.mat)

  log.qpcr <- log(qpcr.mat)
  log.qpcr[log.qpcr == -Inf] <- 0

  # For each comparison
  qpcr.abundance <- vapply(comparisons, function(x) {
    # Subset pair of samples
    psamp <- log.qpcr[, x]
    abundance.unit <- as.character(NA)
    if (all(psamp != 0)) {
      # If log samp1 is 0.05 smaller than samp2
      if (abs(psamp[1] - psamp[2]) > 0.05 & psamp[1] < psamp[2]) {
        abundance.unit <- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",up).")
      }
      # If log samp1 is 0.5 bigger than samp2
      if (abs(psamp[1] - psamp[2]) > 0.05 & psamp[1] > psamp[2]) {
        abundance.unit <- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",down).")
      }
    } else {
      if (exclusion) {
        if (psamp[1] == 0 & psamp[2] > 0) {
          abundance.unit <- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",app).")
        }
        if (psamp[1] > 0 & psamp[2] == 0) {
          abundance.unit <- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",dis).")
        }
      }
    }

    return(abundance.unit)
  }, FUN.VALUE = character(1))
  # Delete NAs
  qpcr.abundance <- qpcr.abundance[!is.na(qpcr.abundance)]
  return(qpcr.abundance)
}

getPresence <- function(tb) {
  # Iterate by ASV table rows
  presence <- lapply(seq_len(nrow(tb)), function(x) {

    # Check presence in each sample
    yes.no <- ifelse(as.numeric(tb[x, ]) > 0, "yes", "no")

    # Build character vector
    presence <- paste0("presence(", colnames(tb), ",", rownames(tb)[x], ",", unname(yes.no), ").")
    return(presence)
  })
  return(unlist(presence))
}


formatOutput <- function(abduced.tables) {

  # Split effect and species
  eff <- vapply(strsplit(abduced.tables[, 1], split = "\\("), function(x) x[1], FUN.VALUE = character(1))
  s1 <- vapply(strsplit(abduced.tables[, 1], split = "\\("), function(x) x[2], character(1))
  s2 <- vapply(strsplit(s1, split = ","), function(x) x[2], character(1))

  # Format species
  s2 <- gsub(")", "", s2)
  s1 <- vapply(strsplit(s1, split = ","), function(x) x[1], character(1))

  # Create data.frame
  df <- data.frame(s1, s2, eff, abduced.tables[, 2], stringsAsFactors = FALSE)
  colnames(df) <- c("sp1", "sp2", "lnk", "comp")
  df$comp <- as.numeric(df$comp)
  return(df)
}

normaliseCompression <- function(abduced.table, presence) {
  presence <- strsplit(presence, ",")
  yes.no <- gsub(").", "", vapply(presence, function(x) {
    x[3]
  }, FUN.VALUE = character(1)))
  samp <- gsub("presence\\(", "", vapply(presence, function(x) {
    x[1]
  }, FUN.VALUE = character(1)))
  spec <- vapply(presence, function(x) {
    x[2]
  }, character(1))
  no.table <- table(paste0(spec, yes.no))

  spec <- spec[yes.no == "yes"]
  samp <- samp[yes.no == "yes"]

  cooc <- character()
  for (z in seq_len(length(unique(samp)))) {
    sm <- samp[samp == unique(samp)[z]]
    if (length(sm) > 1) {
      sm.spec <- spec[samp == sm[1]]
      cmb <- t(combn(sm.spec, 2))
      cooc <- c(cooc, paste(cmb[, 1], cmb[, 2], sep = ""), paste(cmb[, 2], cmb[, 1], sep = ""))
    }
  }
  coo.tab <- table(cooc)
  abduced.table[, 4] <- apply(abduced.table, 1, function(x) {
    as.numeric(x[4]) * abs(log(coo.tab[paste0(x[1], x[2])] / no.table[paste0(x[1], "no")]))
  })
  abduced.table[is.na(abduced.table[, 4]), 4] <- 0
  return(abduced.table)
}

Abduce <- function(bottom, hypothesis, abundance = NULL) {
  if (is.null(abundance)) {
    abundance <- bottom$abundance
  }

  # Define abducibles
  abducible <- c("effect_up", "effect_down")

  # Source python function
  loadInfIntE()

  # Execute abduction using InfIntE
  coverage <- abduction(bottom$clauses, abducible, positive_example_list = abundance, constant_set = bottom$const, meta_rule = hypothesis, metric = "predictive_power")

  # Extract compression values from python object
  compressions <- vapply(names(coverage), function(x) {
    coverage[[x]]
  }, FUN.VALUE = numeric(1))

  # Join compression with interaction names
  coverage <- cbind(names(coverage), unname(compressions))

  # Format the output into a table
  coverage <- formatOutput(coverage)

  return(coverage)
}

getBottom <- function(tb, depth = NULL, exclusion = FALSE, cores = 1, qpcr = NULL, search.depth = 2) {
  if (!is.null(qpcr)) {
    qpcr <- matrix(qpcr, nrow = 1)
    rownames(qpcr) <- "pathogen"
    tb <- rbind(tb, qpcr)
  }

  # Get observations and presence
  abundance <- getObservations(tb, depth, exclusion, cores)
  presence <- getPresence(tb)
  # Get constants
  const <- c("yes", "no", "up", "down", "app", "dis")
  const <- c(const, rownames(tb))

  presence1 <- gsub("presence", "presence1", presence)

  # Source InfIntE
  loadInfIntE()

  # Generate bottom clauses
  P <- generate_bottom_clause(c(presence, presence1), const, abundance, NULL, container = "memory", depth = search.depth)

  # Format bottom clause
  P <- reticulate::dict(P, convert = TRUE)

  # Create and object with all the elements necessaty for abduction
  tb[is.na(tb)] <- 0
  bottom <- list(P, abundance, presence, presence1, const, tb)
  names(bottom) <- c("clauses", "abundance", "presence", "presence1", "const", "table")
  return(bottom)
}

#############################################################################################################
############################# Bootstrap ##################################################################
################################################################################################ ################

InfIntEPulsar <- function(tb, lambda, bot, hypothesis, exclusion, depth = NULL) {
  # Subset depth
  if (!is.null(depth)) {
    depth <- depth[rownames(tb)]
  } else {
    depth <- rowSums(tb)
  }
  snames <- row.names(bot$table)
  # Get observations from subsampled table
  bot$abundance <- getObservations(t(data.frame(tb)), depth = depth, exclusion = exclusion)
  # Abduce
  ab <- Abduce(bottom = bot, hypothesis = hypothesis)

  # Obtain final value
  abduced.table <- plyr::ddply(ab, .(sp1, sp2, lnk), summarise, comp = max(comp))
  abduced.table <- plyr::ddply(abduced.table, .(sp1, sp2), summarise, lnk = lnk[comp == max(comp)][1], comp = if (length(comp) > 1) {
    max(comp) - min(comp)
  } else {
    comp
  })

  # Add non interacting asvs
  noi.asvs <- snames[!snames %in% unique(c(abduced.table$sp1, abduced.table$sp2))]
  noi.tab <- data.frame(noi.asvs, noi.asvs, rep("no_effect", length(noi.asvs)), rep(0, length(noi.asvs)))
  colnames(noi.tab) <- colnames(abduced.table)
  abduced.table <- rbind(abduced.table, noi.tab)

  # Construct adjacency matrix
  g <- igraph::graph_from_data_frame(abduced.table[, c(1, 2, 4)])
  ad <- igraph::get.adjacency(g, attr = "comp", sparse = TRUE)

  # For each lambda
  pt <- lapply(lambda, function(lam) {
    # Compression bigger than lambda
    tmp <- ad > lam
    ga <- igraph::graph_from_adjacency_matrix(tmp)
    tmp <- igraph::get.adjacency(ga, sparse = TRUE)

    diag(tmp) <- FALSE
    return(tmp)
  })
  # Return the path over lambda
  return(list(path = pt))
}



InfIntE <- function(otu.table, hypothesis, thresh = 0.01, exclusion = FALSE, qpcr = NULL, depth = NULL, nperms = 50, search.depth = 2) {
  tb <- otu.table

  # Check ASV tables
  tb <- checkASVtable(tb)

  # Save given ASV and sample names
  asv.names <- rownames(tb)
  sample.names <- colnames(tb)

  # Set simplified names
  rownames(tb) <- paste0("s", seq_len(nrow(tb)))
  names(asv.names) <- paste0("s", seq_len(nrow(tb)))

  colnames(tb) <- paste0("c", seq_len(ncol(tb)))

  names(depth) <- colnames(tb)

  # Produce bottom clause
  bot <- getBottom(tb, exclusion = exclusion, qpcr = qpcr, depth = depth, search.depth = search.depth)
  tb <- bot$table

  # Get final values
  ab <- Abduce(bottom = bot, hypothesis = hypothesis)

  # Obtain final value
  ab <- plyr::ddply(ab, .(sp1, sp2, lnk), summarise, comp = max(comp))
  ab <- plyr::ddply(ab, .(sp1, sp2), summarise, lnk = lnk[comp == max(comp)][1], comp = if (length(comp) > 1) {
    max(comp) - min(comp)
  } else {
    comp
  })

  # Length observations
  mx <- length(bot$abundance)

  # Lambda distribution
  lambda <- pulsar::getLamPath(max = mx, min = 0, 50, FALSE)

  # Pulsar execution
  pr <- pulsar::pulsar(t(tb), InfIntEPulsar, fargs = list(lambda = lambda, bot = bot, hypothesis = hypothesis, exclusion = exclusion, depth = depth), rep.num = nperms, lb.stars = TRUE, ub.stars = TRUE, thresh = thresh)

  # Format output to dataframe
  fit.p <- pulsar::refit(pr, criterion = "stars")

  df <- data.frame(igraph::get.edgelist(igraph::graph_from_adjacency_matrix(fit.p$refit$stars)))
  df <- ab[paste0(ab[, 1], ab[, 2]) %in% paste0(df[, 1], df[, 2]), ]

  df <- classifyInteraction(df)
  df<- return_names(df, asv.names)
  resul <- list(selected_interactions = returnNames(df, asv.names), pulsar_result = pr, abduced_table = returnNames(ab, asv.names))

  return(resul)
}

##############################################################################
############################# Classify interactions #########################
#############################################################################
classifyInteraction <- function(tb) {
  right <- paste0(tb[, 1], tb[, 2])
  oposite <- paste0(tb[, 2], tb[, 1])

  interactions <- vapply(right, function(x) {
    if (any(oposite == x)) {
      if (which(oposite == x) < which(right == x)) {
        int <- paste0(tb[which(right == x), 3], "/", tb[which(oposite == x), 3])
      } else {
        int <- "take_out"
      }
    } else {
      int <- tb[which(right == x), 3]
    }
  }, FUN.VALUE = character(1))

  interactions <- gsub("effect_down/effect_down", "competition", interactions)
  interactions <- gsub("effect_up/effect_up", "mutualism", interactions)
  interactions <- gsub("effect_up/effect_down", "predation", interactions)
  interactions <- gsub("effect_down/effect_up", "take_out", interactions)
  interactions <- gsub("^effect_down$", "amensalism", interactions)
  interactions <- gsub("^effect_up$", "commensalism", interactions)

  tb[, 3] <- interactions
  tb <- tb[tb[, 3] != "take_out", ]

  return(tb)
}


assignTaxonomy <- function(tb, phylos) {
  tx <- data.frame(tax_table(phylos))
  tx$Genus <- gsub("g__", "", tx$Genus)
  tx$Species <- gsub("s__", "", tx$Species)

  ge <- paste0(tx$Genus, " ", tx$Species)
  names(ge) <- rownames(tx)

  tb[, 1:2] <- apply(tb[, 1:2], c(1, 2), function(x) {

  })
}
