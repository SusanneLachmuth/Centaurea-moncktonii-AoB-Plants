# Code for Figure 2 and (if applicable) corresponding statistical analyses for publication
# Lachmuth S, Molofsky J, Suda J, Milbrath L, Keller SR (accepted for AoB Plants on 17 June 2019) 
# Associations between genomic ancestry, genome size and capitula morphology in the invasive meadow knapweed hybrid complex 
# (Centaurea ×moncktonii C.E. Britton) in eastern North America

### Clean working directory:
rm(list=ls())
gc()


### Install and load packages:

# install dependencies 
# install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)
# library(devtools)
# install pophelper package from GitHub
# devtools::install_github('royfrancis/pophelper')
library(pophelper)


### Read data:

setwd("D:/data")
# Load individual and reformat meta-data 
metadata<-read.delim("metadata.txt", header=T,stringsAsFactors=F)
metadata2<-metadata[,c(3,1)]
names(metadata2)<-c("Taxon","Pop")
# Rename taxa
metadata2$Taxon<-ifelse(metadata2$Taxon=="nigra","C.cf.nigra",ifelse(metadata2$Taxon=="jacea","C.cf.jacea",ifelse(metadata2$Taxon=="hybridNY","NY hybrids","VT hybrids")))

# Load "Admixture" assignment probabilities  
afiles <- list.files(path="./distruct",full.names=TRUE)
# create a qlist
admixK2_K3<-readQ(afiles)

# Load "NewHybrids" hybrids class assignment probabilities (C. jacea panel)
NH_LDjacea <- list.files(path="./NH/jacea",full.names=TRUE)
# create a qlist
NH_LDjacea <-readQ(NH_LDjacea)

# Load hybrids class assignment probabilities (C. nigra panel)
NH_LDnigra <- list.files(path="./NH/nigra",full.names=TRUE)
# create a qlist
NH_LDnigra <-readQ(NH_LDnigra)

   
##########################################################
### Admixture plot (Fig. 2A to copy paste into graphics program):

# Data formating
str(admixK2_K3)
names(admixK2_K3$K2.txt)<-c("Cluster1", "Cluster2")
names(admixK2_K3$K3.txt)<-c("Cluster1", "Cluster2", "Cluster3")

# Colors
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"
rgb(red=223,green=194,blue=125,names="Lbeige",maxColorValue=255)
# light beige: "#DFC27D" 
col_admix<-(c("#0000E9","#640000","#BF812D"))

# Plot
plotQ(qlist=admixK2_K3[1:2],grplab=metadata2,ordergrp=TRUE,imgoutput="join",subsetgrp=c("C.cf.nigra","C.cf.jacea","NY hybrids","VT hybrids"),outputfilename="Admixture_barplot_190521"
      ,clustercol = col_admix, splab=c("K2","K3"),splabsize=3,barbordersize=0.0,grplabsize=1.1) 


##########################################################
### NewHybrids plot (C. jacea panel (Fig. 2B upper panel) to copy paste into graphics program)

# Colors
col_NH <- c("#0000E9","#640000","#DFC27D" ,"#BF812D","#A6C0FF","#FF9DA3","#6682EC","#AA343B")

# Plot
plotQ(qlist=NH_LDjacea,grplab=metadata2,ordergrp=TRUE,subsetgrp=c("C.cf.nigra","C.cf.jacea","NY hybrids","VT hybrids"),outputfilename="NH_barplot_LDjacea_190521"
      ,clustercol = col_NH, splabsize=3,barbordersize=0.0,grplabsize=1.1)


##########################################################
# NewHybrids plot (C. nigra panel (Fig. 2B lower panel) to copy paste into graphics program)

col_NH <- c("#0000E9","#640000","#DFC27D" ,"#BF812D","#A6C0FF","#FF9DA3","#6682EC","#AA343B")

plotQ(qlist=NH_LDnigra,grplab=metadata2,ordergrp=TRUE,subsetgrp=c("C.cf.nigra","C.cf.jacea","NY hybrids","VT hybrids"),outputfilename="NH_barplot_LDnigra_190521"
      ,clustercol = col_NH, splabsize=3,barbordersize=0.0,grplabsize=1.1)


##########################################################
# Legend (to copy paste into graphics program)

col_legend <- c("#0000E9","#6682EC","#A6C0FF","#DFC27D" ,"#BF812D","#FF9DA3","#AA343B","#640000")
plot(1:8,pch=19,col=col_legend)
legend("topleft",col=col_legend, pch=19,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F1","F2","BC1nigra","BC2nigra","C.cf.nigra"),bty="n",cex=0.8)

