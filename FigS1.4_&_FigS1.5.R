# Code for Figures S1.4, S1.5 and (if applicable) corresponding statistical analyses for publication
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
library(ggbiplot)
library(multcomp)


### Read and format data:
setwd("D:/data")
gen_dat<-read.table("pop_genetic_data.txt",header = T, sep="\t")
summary(gen_dat)
# Order hybrid class levels
gen_dat$HybClass_LDnigra<-factor(gen_dat$HybClass_LDnigra,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC1_nigra","BC2_nigra","nigra"))

morpho_dat <- read.table("PCA_data_CnigraPanel.txt", header=T, sep = "\t") 
summary(morpho_dat)
# Order hybrid class levels:
morpho_dat$HybClass_LDnigra<- factor(morpho_dat$HybClass_LDnigra,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC1_nigra","BC2_nigra","nigra"))

###########################################################
### Analyze relationships between genome size and genomic as well as morphometric PC-scores:

# genome size ~ genomic PC-scores (axis 1)
lmer1<-lmer(GenSize~PC1+(1|pop),data=gen_dat,REML=T,na.action=na.omit) 
dropterm(lmer1,test="Chisq") 
summary(lmer1)
#plot(lmer1) 
#hist(resid(lmer1)) 

# genome size ~ genomic PC-scores (axis 2)
lmer2<-lmer(GenSize~PC2+(1|pop),data=gen_dat,REML=F) 
dropterm(lmer2,test="Chisq")
summary(lmer2)
#plot(lmer2) 
#hist(resid(lmer2)) 


# Different data set: Mean values (of morphometric PCA scores) per indivdual:

# genome size ~ morphometric PC-scores (axis 1)
lmer3<-lmer(GenSize~PC1_morpho+(1|pop),data=morpho_dat,REML=F) 
dropterm(lmer3,test="Chisq") 
#plot(lmer3) 
#hist(resid(lmer3)) 

# genome size ~ morphometric PC-scores (axis 2)
lmer4<-lmer(GenSize~PC2_morpho+(1|pop),data=morpho_dat,REML=F) 
dropterm(lmer4,test="Chisq")
#plot(lmer4) 
#hist(resid(lmer4)) 


###########################################################
### Figures:

# Reorder factor levels:
gen_dat$HybClass_LDnigra<- factor(gen_dat$HybClass_LDnigra,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC1_nigra","BC2_nigra","nigra"))

# Colors and symbols
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"
colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",2),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#FF9DA3","#AA343B","#640000")    
symbols<- c(15,23,22,19,24,25,17)                                          


###########################################################
### Figure S1.4:

tiff("FigS1.4.tiff",width=3.5,height=3.5,units="in",res=800)

par(mar=c(5,4,1,2),mgp=c(2.5,0.5,0),cex.axis=0.8, cex.lab=1.1,cex=0.6,lwd=1) 

with(gen_dat,plot(y=GenSize,x=PC2,xlab="Genetic principal component 2",ylab="1C genome size (pg)",pch=symbols[HybClass_LDnigra],col=colors[HybClass_LDnigra],bg=bg[HybClass_LDnigra])) 

legend(3.8,2.17,col=colors, pch=symbols,pt.bg=bg,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")


dev.off()

###########################################################
### Figure S1.5

tiff("FigS1.5.tiff",width=3.5,height=7,units="in",res=800)

par(mfrow=c(2,1),mar=c(5,4,1,2),cex.axis=0.8, cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1) 

# Different data set: Mean values (of morphometric PCA scores) per indivdual:
with(morpho_dat,plot(y=GenSize,x=PC1_morpho,xlab="Morphometric principal component 1",ylab="1C genome size (pg)",pch=symbols[HybClass_LDnigra],col=colors[HybClass_LDnigra],bg=bg[HybClass_LDnigra])) 
text(-2.8,2.18,"(A)")

legend(-3,2.16,col=colors, pch=symbols,pt.bg=bg,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

with(morpho_dat,plot(y=GenSize,x=PC2_morpho,xlab="Morphometric principal component 2",ylab="1C genome size (pg)",pch=symbols[HybClass_LDnigra],col=colors[HybClass_LDnigra],bg=bg[HybClass_LDnigra])) 
text(-2.2,2.18,"(B)")


dev.off()
