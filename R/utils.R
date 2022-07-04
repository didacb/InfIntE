
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

return_names <- function(infered, asv.names) {
  names(asv.names)<- paste0("s", seq_len(length(asv.names)))
  infered[, 1:2] <- apply(infered[, 1:2], c(1, 2), function(x) ifelse(grepl("^s", x), asv.names[x], x))
  return(infered)
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

get_comparsions<- function(n_samples){
  # All samples to compare
  comparisons <- combn(n_samples, 2, NULL, FALSE)
  # Add the opposite direction
  other.direction <- lapply(comparisons, function(x) {
    c(x[2], x[1])
  })
  comparisons <- c(comparisons, other.direction)
  return(comparisons)
}

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
