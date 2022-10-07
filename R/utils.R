
#' Check OTU table
#' Checks the presence of NA and non-numeric cells
#' @param otu_tab data frame with the OTUs counts with samples on columns and OTUs in rows.
#'
#' @return revised data frame
#' @export
#'
#' @examples
#' check_OTU_table(otu_tab)
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

#' Return names
#' Returns de original names to the OTUs
#' @param abduced a data frame, containing the significant effects where each row contains a pair of OTU, the effect on the abundance or interaction and the I value supporting the effect.
#' @param otu_names named character vector with the original otu names in the same order as provided in the original input.
#'
#' @return the same input data frame substituting the executing names by the original ones.
#' @export
#'
#' @examples
#' return_names(infered, otu_names)
return_names <- function(infered, otu_names) {
  names(otu_names) <- paste0("s", seq_len(length(otu_names)))
  infered[, 1:2] <- apply(infered[, 1:2], c(1, 2), function(x) {
    ifelse(grepl("^s", x), otu_names[x], x)
  })
  return(infered)
}


#' Format abduction output
#' Gives data frame format to the raw abduction output produced by PyGol
#' @param abduced_raw data frame containing the raw output of PyGol abduction and the compression values.
#'
#' @return data frame where each row contains a pair of OTU, the effect on the abundance and the compression supporting the effect.
#' @export
#'
#' @examples
#' formatOutput(abduced_raw)
formatOutput <- function(abduced_raw) {

  # Split effect and species
  eff <- vapply(strsplit(abduced_raw[, 1], split = "\\("), function(x) {
    x[1]
  }, FUN.VALUE = character(1))
  s1 <- vapply(strsplit(abduced_raw[, 1], split = "\\("), function(x) {
    x[2]
  }, character(1))
  s2 <- vapply(strsplit(s1, split = ","), function(x) x[2], character(1))

  # Format species
  s2 <- gsub(")", "", s2)
  s1 <- vapply(strsplit(s1, split = ","), function(x) x[1], character(1))

  # Create data.frame
  df <- data.frame(s1, s2, eff, abduced_raw[, 2], stringsAsFactors = FALSE)
  colnames(df) <- c("sp1", "sp2", "lnk", "comp")
  df$comp <- as.numeric(df$comp)
  return(df)
}

#' Get comparisons
#' Produces a list with all the possible comparisons between samples in both directions.
#' @param n_samples numeric vector with the position of tha samples to compare.
#'
#' @return a list with the pairs of samples to compare.
#' @export
#'
#' @examples
#' get_comparsions(n_samples)
get_comparsions <- function(n_samples) {
  # All samples to compare
  comparisons <- combn(n_samples, 2, NULL, FALSE)
  # Add the opposite direction
  other.direction <- lapply(comparisons, function(x) {
    c(x[2], x[1])
  })
  comparisons <- c(comparisons, other.direction)
  return(comparisons)
}

#' Load PyGol
#' Loads PyGol function for abduction using reticulate.
#' @return
#' @export
#'
#' @examples
#' load_PyGol()
load_PyGol <- function() {
  package.location<-  system.file(package = "InfIntE")
  python.file<- file.path(package.location, "python", "pygol_abduce.py")
  current_wd <- getwd()

  #install python3-dev
  #install cython
  setwd(file.path(package.location, "python"))

  python.ver<- system2("python3", "--version", stdout = TRUE)
  python.ver<- gsub(" ", "", python.ver)
  python.ver<- gsub("P", "p", python.ver)
  python.ver<- gsub("\\.[1-9]$", "", python.ver)

  if(!file.exists(file.path(package.location, "python", "pygolm_V1.so"))){

    system2("gcc", args = paste("-I", paste0("/usr/include/", python.ver), "-c", "-fPIC", "pygolm_V1.c", "-o" ,"pygolm_V1.o", sep = " "))
    system2("gcc", args = c( "pygolm_V1.o", "-shared", "-o", "pygolm_V1.so"))
  }

  reticulate::use_python(system2("which", python.ver, stdout = TRUE))

  reticulate::source_python(python.file, envir = globalenv())
  setwd(current_wd)
  on.exit(setwd(current_wd))
}


