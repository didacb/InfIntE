
loadPyGolM<- function(){
  #package.location<-  system.file(package = "pyGolm-nets")
  package.location<- "/home/didac/Desktop/PyGolM-nets/str"
  ex_wd <- getwd()
  on.exit(setwd(ex_wd))
  setwd(package.location)
  reticulate::source_python(file.path(package.location, "pygolm_abduce.py"), envir = globalenv())
  setwd(ex_wd)
  #return("PyGolM loaded")
}


checkASVtable<- function(tb){
  #Check if all cells are numeric
  tb<- apply(tb, c(1,2), as.numeric)
  if(any(!is.numeric(tb))){
    return(stop("Please provide an ASV table with numeric values"))  
  }
  #Check if there are NA
  if(any(is.na(tb))){
    return(stop("Please provide an ASV table without NA"))
  }
  return(tb)
}

#Function to apply chisq test to all replicate comparison for one species and obtain the abundance progol input
abundanceChi<- function(comparisons, asv.table, spec, samps, read.depth, exclusion){
  #Detect the pair of values and depth to compare
  pair<- asv.table[spec, comparisons]
  depth<- read.depth[comparisons]
  
  #If both pair values no zero and the samples are not equal to depth 
  if(any(pair != depth)){
    if (all(pair != 0)){
      
      #Make table and add margins
      tab<- as.matrix(rbind(pair, depth - pair))
      tab<- addmargins(A=tab, margin= seq_along(dim(tab)))
      
      #Do the cishq test  
      tst<- suppressWarnings(chisq.test(tab))
      if (tst[3] < 0.05){
        #If it is significant establish the up and down 
        if (pair[1]/depth[1] < pair[2]/depth[2]){
          up.down<- "up"
        } else {
          up.down<- "down"
        }
        #Build logical abundance input
        abundance.unit<- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",", up.down, ").")
      } else {
        #Build logical abundance when no significant
        #abundance.unit<- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",zero).")
        abundance.unit<- as.character(NA)
      }
    }else{
      if(any(pair != 0) & exclusion){
        #Make table and add margins
        tab<- as.matrix(rbind(pair, depth - pair))
        tab<- addmargins(A=tab, margin= seq_along(dim(tab)))
        
        #Do the chisq test  
        tst<- suppressWarnings(chisq.test(tab))
        if (tst[3] < 0.05){
          #If it is significant establish the up and down 
          if (pair[1]/depth[1] < pair[2]/depth[2]){
            up.down<- "app"
          } else {
            up.down<- "dis"
          }
          #Build logic abundance input
          abundance.unit<- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",", up.down, ").")
        }else{
          #Build logic abundance when no significant
          #abundance.unit<- paste0("abundance(", samps[comparisons[1]], ",", samps[comparisons[2]], ",", paste0("s", spec), ",zero).")
          abundance.unit<- as.character(NA)
        } 
      }else{
        abundance.unit<- as.character(NA)  
    }  
  } 
  }else{
    abundance.unit<- as.character(NA)  
  }  
  return(abundance.unit)
}


getObservations<- function(tb, depth=NULL, exclusion, cores=1){
  
  if(is.null(depth)){
    #Obtain read depth
    read.depth<- colSums(tb)
  }else{
    read.depth<- depth
  }
  pathogen<-NULL
  if(any(rownames(tb)== "pathogen")){
    pathogen<- tb["pathogen",,drop=FALSE]
    tb<- tb[rownames(tb) != "pathogen",]
    tb<- data.frame(apply(tb,c(1,2),as.numeric))
  }
    
  #All samples to compare
  comparisons<- combn(seq_len(ncol(tb)),2,NULL,FALSE)
  
  #Add the opposite direction
  other.direction<- lapply(comparisons, function(x){c(x[2], x[1])})
  comparisons<- c(comparisons, other.direction)
  
  #Run pairwise test function for all species
  abundance<- parallel::mclapply(seq_len(nrow(tb)), function(i){
    #Run pairwise
    abundance.specie<- vapply(comparisons, abundanceChi, tb, i, colnames(tb), as.numeric(read.depth), exclusion, FUN.VALUE = character(1))
    #Delete NAs
    abundance.specie<- abundance.specie[!is.na(abundance.specie)]
    return(abundance.specie)
    #Number of cores
  },mc.cores = cores)
  
  #Unlist and return
  abundance<- unlist(abundance)
  
  if(!is.null(pathogen)){
    pathogen<- apply(pathogen, c(1,2), function(x){
      x <-as.numeric(x)
      x<- ifelse(is.na(x), 0, x)
      return(x)
    })
    
    abundance<- c(abundance, getObservationsQpcr(pathogen, exclusion))
  }
  
  return(abundance)
}

getObservationsQpcr<- function(qpcr.mat, exclusion=FALSE){
  
  #Compare all qpcr samples
  comparisons<- combn(seq_len(ncol(qpcr.mat)),2,NULL,FALSE)
  #Add the opposite direction
  other.direction<- lapply(comparisons, function(x){c(x[2], x[1])})
  comparisons<- c(comparisons, other.direction)
  
  #Sample names
  sn<- colnames(qpcr.mat)
  
  log.qpcr<- log(qpcr.mat) 
  log.qpcr[log.qpcr == -Inf]<- 0
  
  #For each comparison
  qpcr.abundance<- vapply(comparisons, function(x){
    #Subset pair of samples
    psamp<- log.qpcr[,x]
    abundance.unit<- as.character(NA)
    if(all(psamp != 0)){
      #If log samp1 is 0.05 smaller than samp2
      if(abs(psamp[1]- psamp[2]) > 0.05 & psamp[1] < psamp[2]){
        abundance.unit<- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",",rownames(qpcr.mat), ",up).")
      }
      #If log samp1 is 0.5 bigger than samp2
      if(abs(psamp[1]- psamp[2]) > 0.05 & psamp[1] > psamp[2]){
        abundance.unit<- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",down).")
      }  
    }else{
      if(exclusion){
        if(psamp[1] == 0 & psamp[2] > 0){
          abundance.unit<- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",app).")
        }
        if(psamp[1] > 0 & psamp[2] == 0){
          abundance.unit<- paste0("abundance(", sn[x[1]], ",", sn[x[2]], ",", rownames(qpcr.mat), ",dis).")
        }
      }
    }
    
    return(abundance.unit)
  }, FUN.VALUE = character(1))  
  #Delete NAs
  qpcr.abundance<- qpcr.abundance[!is.na(qpcr.abundance)]
  return(qpcr.abundance)
}

getPresence<- function(tb){
  #Iterate by ASV table rows
  presence<- lapply(seq_len(nrow(tb)), function(x){
    
    #Check presence in each sample
    yes.no<- ifelse(as.numeric(tb[x,])>0, "yes", "no")
    
    #Build character vector
    presence<- paste0("presence(", colnames(tb), ",", rownames(tb)[x],",", unname(yes.no), ").")
    return(presence)
  })
  return(unlist(presence))
}

getRelativePresence<- function(tb, depth=NULL){
  
  if(!is.null(depth)){
    depth<- colSums(tb)
  }
  #Iterate by ASV table rows
  rel.presence<- lapply(seq_len(nrow(tb)), function(x){
    
    #Check presence in each sample
    hilow<- ifelse(as.numeric(tb[x,])>(depth/2), "high", "low")
    
    #Build character vector
    relpres<- paste0("relative(", colnames(tb), ",", rownames(tb)[x],",", unname(hilow), ").")
    relpres<- relpres[tb[x,]!=0]
    return(relpres)
  })
  return(unlist(rel.presence))
}


formatOutput<- function(abduced.tables){
  
  #Split effect and species
  eff<- vapply(strsplit(abduced.tables[,1], split = "\\("), function(x)x[1],FUN.VALUE = character(1))
  s1<- vapply(strsplit(abduced.tables[,1], split = "\\("), function(x)x[2],character(1))
  s2<- vapply(strsplit(s1, split = ","), function(x)x[2],character(1))
  
  #Format species
  s2<- gsub(")", "", s2)
  s1<- vapply(strsplit(s1, split = ","), function(x)x[1], character(1))
  
  #Create data.frame
  df<- data.frame(s1,s2,eff,abduced.tables[,2], stringsAsFactors = FALSE)
  colnames(df)<- c("sp1", "sp2", "lnk", "comp")
  df$comp<- as.numeric(df$comp)
  return(df)
}

normaliseCompression<- function(abduced.table, presence){
  presence<- strsplit(presence, ",")
  yes.no<- gsub(").", "", vapply(presence, function(x){x[3]}, FUN.VALUE = character(1)))
  samp<- gsub("presence\\(", "", vapply(presence, function(x){x[1]}, FUN.VALUE = character(1)))
  spec<- vapply(presence, function(x){x[2]}, character(1))
  no.table<- table(paste0(spec, yes.no))
  
  spec<- spec[yes.no == "yes"]
  samp<- samp[yes.no == "yes"]
  
  cooc<- character()
  for (z in seq_len(length(unique(samp)))){
    
    sm<- samp[samp==unique(samp)[z]]
    if (length(sm)>1){
      sm.spec<- spec[samp ==sm[1]]
      cmb<-t(combn(sm.spec, 2))
      cooc<- c(cooc, paste(cmb[,1], cmb[,2], sep = ""), paste(cmb[,2], cmb[,1], sep = ""))
    }
  }
  coo.tab<- table(cooc)
  abduced.table[,4]<- apply(abduced.table, 1, function(x){ as.numeric(x[4])*abs(log(coo.tab[paste0(x[1],x[2])] / no.table[paste0(x[1], "no")]))})
  abduced.table[is.na(abduced.table[,4]),4]<- 0
  return(abduced.table)
}

Abduce<- function(bottom, hypothesis, abundance=NULL){
  
  if(is.null(abundance)){
    abundance<-bottom$abundance
  }
  
  #Define abducibles
  abducible<- c('effect_up','effect_down')
  
  #Source python function
  loadPyGolM()
  
  #Execute abduction using PyGolM
  coverage<- abduction(bottom$clauses, abducible,  positive_example_list=abundance, constant_set=bottom$const, meta_rule=hypothesis, metric="predictive_power")
  
  #Extract compression values from python object
  compressions<- vapply(names(coverage), function(x){coverage[[x]]}, FUN.VALUE = numeric(1))
  
  #Join compression with interaction names 
  coverage<- cbind(names(coverage), unname(compressions))
  
  #Format the output into a table
  coverage<- formatOutput(coverage)
  
  return(coverage)
  
}

getBottom<- function(tb, depth=NULL, exclusion=FALSE, cores=1, qpcr=NULL,search.depth=2){
  
  if(!is.null(qpcr)){
    qpcr<- matrix(qpcr, nrow=1)
    rownames(qpcr)<- "pathogen"
    tb<- rbind(tb, qpcr)
  }
  
  #Get observations and presence
  abundance<- getObservations(tb, depth, exclusion, cores)
  presence<- getPresence(tb)
  #Get constants
  const<- c("yes", "no", "up", "down", "app", "dis")
  const<- c(const, rownames(tb))
  
  presence1<- gsub("presence", "presence1", presence)
  
  #Source PyGolM
  loadPyGolM()
  
  #Generate bottom clauses
  P<-  generate_bottom_clause(c(presence, presence1), const,abundance, NULL,  container="memory", depth=search.depth)
  
  #Format bottom clause
  P<- reticulate::dict(P, convert = TRUE)
  
  #Create and object with all the elements necessaty for abduction
  tb[is.na(tb)]<- 0
  bottom<- list(P, abundance, presence, presence1, const, tb)
  names(bottom)<- c("clauses", "abundance", "presence", "presence1", "const", "table")
  return(bottom) 
}

#############################################################################################################
############################# Bootstrap ##################################################################
################################################################################################ ################

pygolmPulsar<- function(tb, lambda, bot, hypothesis, exclusion, depth=NULL){
  #Subset depth
  if(!is.null(depth)){
    depth<- depth[rownames(tb)]
  }else{
    depth<-rowSums(tb)
  }
  snames<- row.names(bot$table)
  #Get observations from subsampled table
  bot$abundance<- getObservations(t(data.frame(tb)), depth = depth, exclusion = exclusion)
  #Abduce
  ab<- Abduce(bottom = bot, hypothesis = hypothesis)
  
  #Obtain final value
  abduced.table<- plyr::ddply(ab, .(sp1, sp2, lnk), summarise, comp=max(comp))   
  abduced.table<- plyr::ddply(abduced.table, .(sp1, sp2), summarise, lnk = lnk[comp == max(comp)][1], comp = if(length(comp)>1){max(comp) - min(comp)}else{comp})
  
  #Add non interacting asvs
  noi.asvs<- snames[!snames %in% unique(c(abduced.table$sp1, abduced.table$sp2))]
  noi.tab<- data.frame(noi.asvs, noi.asvs, rep("no_effect", length(noi.asvs)), rep(0, length(noi.asvs)))
  colnames(noi.tab)<- colnames(abduced.table)
  abduced.table<- rbind(abduced.table, noi.tab)
   
  #Construct adjacency matrix    
  g<- igraph::graph_from_data_frame(abduced.table[,c(1,2,4)])
  ad<- igraph::get.adjacency(g, attr = "comp",sparse = TRUE)
  
  #For each lambda
  pt<- lapply(lambda, function(lam){
    #Compression bigger than lambda
    tmp <- ad > lam
    ga<- igraph::graph_from_adjacency_matrix(tmp)
    tmp<-igraph::get.adjacency(ga, sparse = TRUE)
    
    diag(tmp) <- FALSE
    return(tmp)
  })
  #Return the path over lambda
  return(list(path=pt))
}

getPvals<-function(abduction.result){
  
  final.tabs<- lapply(abduction.result, function(x){
    x[,4]<- as.numeric(x[,4])
    x<- x[x[,1] != x[,2],]
    df<- ddply(x, .(sp1, sp2), summarise, comp = if(length(comp)>1){
      comp[grepl("_up", lnk)] - comp[grepl("_down", lnk)]}else{comp * ifelse(grepl("_up", lnk), 1, -1)})
    return(df)
  })
  
  final.tab<- Reduce(function(x, y) merge(x, y, by=c(1,2), all=TRUE), final.tabs)
  
  final.tab<- apply(final.tab, c(1,2), function(x){ifelse(is.na(x), 0, x)})
  
  colnames(final.tab)<- c("sp1", "sp2", "correct", paste0("p", seq_len(ncol(final.tab)-3)))
  
  pvals<- apply(final.tab, 1, function(r){
    correct<- as.numeric(r[3])
    vals<- as.numeric(r[4:length(r)])
    if (correct != 0){
      if(correct < 0){
        pval<- mean(vals < correct)
      }else{
        pval<- mean(vals > correct)
      }
    }else{
      pval<-1
    }
    return(pval)
  })
  names(pvals)<- paste0(final.tab[,1], "-", final.tab[,2])
  return(pvals)
}



boot.fun<- function(df, i){
  df1<- df
  
  #Randomize interactions 
  df1[,5]<- df1[i,5]
  
  #Keep new vallues
  df1<- df1[df1[,5]=="yes",] 
  
  #Take only ups or downs
  up<- df1[df1[,3]=="effect_up",4]
  dn<- df1[df1[,3]=="effect_down",4]
  
  #Obtain I statistic value
  if(length(up)==0){
    up<-0
  }else{
    up<- max(up)
  }
  if(length(dn)==0){
    dn<-0
  }else{
    dn<- max(dn)
  }
  return(up-dn)
}


bootstrapCompression<- function(rel.table, cores=1){
  
  #Make a tag for each interaction
  boot.tag<- paste(rel.table[,1], rel.table[,2], sep = "-")
  rel.table$comp<- as.numeric(rel.table$comp)
  #Obtain pvalue for each possible interaction
  pvals<- parallel::mclapply(sort(unique(boot.tag)), function(x){
    
    #Obtain both interacting ASVs
    sp1<- vapply(strsplit(x, split = "-"), function(y){y[1]}, FUN.VALUE = character(1)) 
    sp2<- vapply(strsplit(x, split = "-"), function(y){y[2]}, FUN.VALUE = character(1)) 
    
    #Subset data frame to keep interactions involving at least one of the species
    rl<- rel.table[rel.table[,1] %in% c(sp1,sp2) | rel.table[,2] %in% c(sp2,sp1),] 
    
    #Make interaction tag for the subsetetd dataframe
    tag<- paste(rl[,1], rl[,2], sep = "-")                  
    
    #Mark interactions involving both species
    rl<- cbind(rl, ifelse(tag==x, "yes", "no"))
    #Set colnames
    colnames(rl)<- c("sp1","sp2", "lnk", "comp", "correct")
    
    #Bootstrap
    bt<-boot::boot(data = rl, statistic = boot.fun, R = 999, sim = "ordinary",parallel = "no")
    
    #Calculate the p-value
    if(bt$t0 < 0){
      pval<- mean(bt$t < bt$t0)  
    }else{
      pval<- mean(bt$t > bt$t0)
    }
    names(pval)<- x
    return(pval)
  }, mc.cores = cores)
  pvals<- unlist(pvals)
  return(pvals)
}

checkInteractions<- function(dtf, ns){
  sp1<-character()
  sp2<- character()
  for (i in seq_len(ns/2)){
    sp1<- c(sp1, paste0("s", i*2-1))
    sp2<- c(sp2, paste0("s", i*2))
  }
  ct<- paste0(sp2, sp1)
  
  sp1<-character()
  sp2<- character()
  for (i in c(6:10, 16:20)){
    sp1<- c(sp1, paste0("s", i*2-1))
    sp2<- c(sp2, paste0("s", i*2))
  }
  
  ct<- c(ct, paste0(sp1,sp2))
  dt<- paste0(dtf[,1], dtf[,2])
  return(paste0("TP: ", sum(dt %in% ct), " (negative: ", sum(grepl("down", dtf[dt %in% ct,3])),"), FP: ", sum(! dt %in% ct), ", total: ", length(dt)))
}

checkSparsity<- function(tb, plotg=FALSE){
  sps<- apply(tb, 1, function(lin){
    sum(lin != 0)/length(lin)
  })
  if(plotg){
    return(plot(sort(sps)))
  }else{
    return(sps)  
  }
}

PyGolMnets<- function(otu.table, hypothesis, thresh=0.01, exclusion=FALSE, qpcr=NULL, depth=NULL, nperms=50, search.depth=2){
  tb<- otu.table
  
  #Check ASV tables
  tb<- checkASVtable(tb)
  
  #Save given ASV and sample names
  asv.names<- rownames(tb)
  sample.names<- colnames(tb)
  
  #Set simplified names
  rownames(tb)<- paste0("s", seq_len(nrow(tb)))
  colnames(tb)<- paste0("c", seq_len(ncol(tb)))
  
  #Produce bottom clause
  bot<- getBottom(tb, exclusion = exclusion, qpcr = qpcr,depth = depth, search.depth = search.depth)
  tb<- bot$table
  
  #Get final values
  ab<- Abduce(bottom = bot,hypothesis =  hypothesis)
  #Obtain final value
  ab<- plyr::ddply(ab, .(sp1, sp2, lnk), summarise, comp=max(comp))
  ab<- plyr::ddply(ab, .(sp1, sp2), summarise, lnk = lnk[comp == max(comp)][1], comp = if(length(comp)>1){max(comp) - min(comp)}else{comp})
  
  #Length observations
  mx<- length(bot$abundance)
  
  #Lambda distribution
  lambda<-pulsar::getLamPath(max = mx, min = 0, 50, FALSE)
  
  #Pulsar execution
  pr<- pulsar::pulsar(t(tb), pygolmPulsar, fargs = list(lambda=lambda, bot=bot, hypothesis=hypothesis, exclusion=exclusion), rep.num = nperms, lb.stars = TRUE,ub.stars = TRUE, thresh = thresh)
  
  #Format output to dataframe
  fit.p <- pulsar::refit(pr, criterion = "stars")
  
  df<- data.frame(igraph::get.edgelist(igraph::graph_from_adjacency_matrix(fit.p$refit$stars)))
  df<- ab[paste0(ab[,1], ab[,2]) %in% paste0(df[,1], df[,2]),]
  
  df<- classifyInteraction(df)
  
  resul<-list(selected_interactions=df, pulsar_result=pr, abduced_table=ab)

  return(resul)
}

##############################################################################
############################# Classify interactions #########################
#############################################################################
classifyInteraction<- function(tb){
  right<- paste0(tb[,1], tb[,2])
  oposite<- paste0(tb[,2], tb[,1])
  
  interactions<- vapply(right, function(x){
    if(any(oposite == x)){
      if(which(oposite == x) < which(right == x)){
        int<- paste0(tb[which(right == x), 3], "/", tb[which(oposite == x), 3])
      }else{
        int<- "take_out"
      }
    }else{
      int<- tb[which(right == x), 3]
    }
  }, FUN.VALUE = character(1))
  
  interactions<- gsub("effect_down/effect_down", "competition", interactions)
  interactions<- gsub("effect_up/effect_up", "mutualism", interactions)
  interactions<- gsub("effect_up/effect_down", "predation", interactions)
  interactions<- gsub("effect_down/effect_up", "take_out", interactions)
  interactions<- gsub("^effect_down$", "amensalism", interactions)
  interactions<- gsub("^effect_up$", "commensalism", interactions)
  
  tb[,3]<- interactions
  tb<- tb[tb[,3] != "take_out",]
  
  return(tb)
}

