# InfintE
Microbial interaction inference using logic machine learning

#Dependencies
##R packages
pulsar
igraph
plyr

##Python
pygolm_V1

#Basic usage in R

Define PyGolM location
pyGolM.location<- "~/Desktop/pygolm_V1.py"

Define interaction hypothesis
hypothesis<- c("abundance(C1,C2,S1,up):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_up(S2,S1)",                
                 "abundance(C1,C2,S1,down):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_down(S2,S1)")
  
OTU/ASV table with OTUs/ASVs in rows and samples in columns
otu.table<- matrix(trunc(rlnorm(20*20, meanlog = 0.5, sdlog = 3)),nrow=20)

Iteraction inference using PyGolM-nets
inference<- PyGolMnets(otu.table, hypothesis, pyGolM.location)

inference object contains the table with the interactions classified, de results of the StARS selection, the result of the abduction.


gcc -I/usr/include/python3.9 I/usr/lib/R  -c -fPIC pygolm_V1.c -o pygolm_V1.o

gcc pygolm_V1.o -shared -o pygolm_V1.so
