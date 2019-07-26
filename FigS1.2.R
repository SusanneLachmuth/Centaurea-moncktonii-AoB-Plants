# Code for Figure S1.2 and (if applicable) corresponding statistical analyses for publication
# Lachmuth S, Molofsky J, Suda J, Milbrath L, Keller SR (accepted for AoB Plants on 17 June 2019) 
# Associations between genomic ancestry, genome size and capitula morphology in the invasive meadow knapweed hybrid complex 
# (Centaurea ×moncktonii C.E. Britton) in eastern North America

### Clean working directory:

rm(list=ls())
gc()


### Load packages:
library(ggplot2)


### Read and format data:
setwd("D:/data")

# Individuals' hybrid classes (C. jacea panel)
HybClass_jacea<-read.table("ind_HybClasses_Cjacea_Panel.txt", header=T,sep="\t") 
summary(HybClass_jacea)
levels(HybClass_jacea$HybClass)
levels(HybClass_jacea$HybClass)<-c("BC1 jacea","BC2 jacea","BC2 nigra","F2","C. cf. jacea","C. cf. nigra")
# Individuals' hybrid classes (C. nigra panel)
HybClass_nigra<-read.table("ind_HybClasses_Cnigra_Panel.txt", header=T,sep="\t") 
levels(HybClass_nigra$HybClass)
levels(HybClass_nigra$HybClass)<-c("BC1 jacea","BC1 nigra","BC2 jacea","BC2 nigra","F2","C. cf. jacea","C. cf. nigra")
summary(HybClass_nigra)


# Merge data
dat<-data.frame(HybClass_jacea=factor(HybClass_jacea$HybClass,levels=c("C. cf. jacea","BC2 jacea","BC1 jacea","F2","BC2 nigra","C. cf. nigra"),ordered=F),HybClass_nigra=factor(HybClass_nigra$HybClass,levels=c("C. cf. jacea","BC2 jacea","BC1 jacea","F2","BC1 nigra","BC2 nigra","C. cf. nigra"),ordered=F))
str(dat)

# Plot
p<-ggplot(dat, aes(dat$HybClass_jacea,dat$HybClass_nigra)) + geom_count()
p+labs(x = "C. cf. jacea data set",y="C. cf. nigra data set")

