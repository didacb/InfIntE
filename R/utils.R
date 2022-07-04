
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


join_abundances<- function(otu_tb, absolute_abundance, depth=NULL){

  change_position<- nrow(otu_tb)

  if(!is.null(absolute_abundance) & !is.null(otu_tb)){

    #Join tables
    colnames(absolute_abundance)<- colnames(otu_tb)
    otu_tb<- rbind(otu_tb, absolute_abundance)

    end_position<- nrow(otu_tb) - change_position
  }

  if(is.null(otu_tb)){
    change_position<- 0
    otu_tb<- absolute_abundance
    end_postion<- nrow(otu_tb)
  }

  # Save given OTU and sample names
  samp_names <- colnames(otu_tb)
  otu_names <- rownames(otu_tb)

  # Set simplified names
  rownames(otu_tb) <- paste0("s", seq_len(nrow(otu_tb)))
  colnames(otu_tb) <- paste0("c", seq_len(ncol(otu_tb)))

  #Check table
  otu_tb <- check_OTU_table(otu_tb)

  if(is.null(depth) & !is.null(otu_tb)){
    depth<- colSums(otu_tb[1:change_position,])
  }

  abundance_function<- c(rep("get_compositional_observations", change_position),
                         rep("get_absolute_observations", end_position))

  names(abundance_function)<- rownames(otu_tb)
  names(depth)<- colnames(otu_tb)

  otu_data<- list("otu_tb"=otu_tb, "otu_names"=otu_names, "samp_names"=samp_names,
                  "depth"=depth,"abundance_function"=abundance_function)
  return(otu_data)
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

get_I_values<- function(abduced){
  # Obtain final value
  abduced <- plyr::ddply(abduced, .(sp1, sp2, lnk), summarise, comp = max(comp))
  abduced <- plyr::ddply(abduced, .(sp1, sp2), summarise,
                               lnk = lnk[comp == max(comp)][1], comp = if (length(comp) > 1) {
                                 max(comp) - min(comp)
                               } else {
                                 comp
                               })
  return(abduced)
}
