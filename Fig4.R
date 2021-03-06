# Code for Figure 4 and (if applicable) corresponding statistical analyses for publication
# Lachmuth S, Molofsky J, Suda J, Milbrath L, Keller SR (accepted for AoB Plants on 17 June 2019) 
# Associations between genomic ancestry, genome size and capitula morphology in the invasive meadow knapweed hybrid complex 
# (Centaurea ×moncktonii C.E. Britton) in eastern North America

### Clean working directory:

rm(list=ls())
gc()


### Install and load packages:

library(stats)
library(MASS)
library(lme4)
library(car)
#library(devtools)
#install_github("ggbiplot", "vqv",force=T)
library(ggbiplot)
library(grid)
library(ggplot2)


### Read data:

setwd("D:/data")
cap <- read.table("capitula_data_PCA.txt", header=T, sep = "\t") 
summary(cap)
str(cap)


###########################################################
### Perform PCA for morphometric trait data: 

pca_dat <- cap[,5:7]
names(pca_dat)<- c("Pectinate bract rows (%)","Appendage center width/length","Rel. width appendage center" ) #
summary(pca_dat)

pca1 <- prcomp(pca_dat,
               center = TRUE,
               scale. = TRUE) 

print(pca1)
summary(pca1)


###########################################################
### Morphometric PCA plot (Figure 4):

# Colors according to  NewHybrids hybrid classes
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
#dark beige: "#BF812D"
colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",2),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#FF9DA3","#AA343B","#640000")    
symbols<- c(15,23,22,19,24,25,17)                                          
HybClass<- factor(cap$HybClass_LDnigra,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC1_nigra","BC2_nigra","nigra"))

#Plot
PCAplot<-ggbiplot(pca1, obs.scale = 1, var.scale = 1, group=HybClass,
         varname.size = 3.7, labels.size=1.5, ellipse = TRUE, circle = TRUE) +
  scale_color_manual(name="Hybrid class", values=colors, labels=c("C. cf. jacea","BC2 jacea","BC1 jacea","F2","BC1 nigra","BC2 nigra","C. cf. nigra")) +  
  scale_shape_manual(name="Hybrid class", values=symbols, labels=c("C. cf. jacea","BC2 jacea","BC1 jacea","F2","BC1 nigra","BC2 nigra","C. cf. nigra")) +
  scale_fill_manual(name="Hybrid class", values=bg, labels=c("C. cf. jacea","BC2 jacea","BC1 jacea","F2","BC1 nigra","BC2 nigra","C. cf. nigra")) +
  geom_point(aes(colour=HybClass, shape=HybClass,fill=HybClass), size = 2) +
  theme(legend.direction ="horizontal", legend.title = element_text(size = 10),
        legend.position = "top",legend.text = element_text(size = 10))+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))

# Save Figure 4
ggsave("Fig4.jpg", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 17, height = NA, units =  "cm",
       dpi = 600, limitsize = TRUE)

dev.off()



