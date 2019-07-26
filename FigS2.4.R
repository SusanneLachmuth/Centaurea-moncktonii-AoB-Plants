# Code for Figure S2.4 and (if applicable) corresponding statistical analyses for publication
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
library(multcomp)
library(multcompView)


### Read data:

setwd("D:/data")
dat <- read.table("PCA_data_CjaceaPanel.txt", header=T, sep = "\t") 
summary(dat)
# Order hybrid class levels:
dat$HybClass_LDjacea<- factor(dat$HybClass_LDjacea,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC2_nigra","nigra"))


###########################################################
### Analyze relationships between morophometric PC-scores and genomic PC-scores:


# morophometric PC-scores (axis 1) ~ genomic PC-scores (axis 1)
lmer1<-lmer(PC1_morpho~PC1+(1|pop),data=dat,REML=F) 
dropterm(lmer1,test="Chisq") 
#plot(lmer1) 

# morophometric PC-scores (axis 1) ~ genomic PC-scores (axis 2)
lmer2<-lmer(PC1_morpho~PC2+(1|pop),data=dat,REML=F) 
dropterm(lmer2,test="Chisq") 
#plot(lmer2) 

# morophometric PC-scores (axis 2) ~ genomic PC-scores (axis 1)
lmer3<-lmer(PC2_morpho~PC1+(1|pop),data=dat,REML=F) 
dropterm(lmer3,test="Chisq")
#(lmer3) 

# morophometric PC-scores (axis 2) ~ genomic PC-scores (axis 2)
lmer4<-lmer(PC2_morpho~PC2+(1|pop),data=dat,REML=F) 
dropterm(lmer4,test="Chisq") 
#plot(lmer4) 


###########################################################
### Figure S2.4:

# Colors and symbols
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D" 
colors<-c("#0000E9",rep("grey20",2),"#BF812D" ,rep("grey20",1),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D" ,"#AA343B","#640000")    
symbols<- c(15,23,22,19,25,17)                                        

# Plot:
tiff("FigS2.4.tiff",width=7,height=7,units="in",res=800)

par(mfrow=c(2,2), cex.axis=0.8, cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1, mar=c(5,4,1,2))


###########################################################
### Scatter plots

# morophometric PC-scores (axis 1) ~ genomic PC-scores (axis 1)
with(dat,plot(y=PC1_morpho,x=PC1,xlab="Genetic principal component 1",ylab="Morphometric principal component 1",pch=symbols[HybClass_LDjacea],bg=bg[HybClass_LDjacea],col=colors[HybClass_LDjacea],ylim=c(-3,3.5),yaxp=c(-3,3,6)))
lines(seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1),predict(lmer1,newdata=data.frame(PC1=seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1)),re.form=NA))
text(-6,3.5,"(A)")
text(6.5,3.5,"chisq (1) = 10.6, p < 0.01")  # korrekt (180622)
legend(-6.5,3.2,col=colors,pt.bg =bg, pch=symbols,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

# morophometric PC-scores (axis 1) ~ genomic PC-scores (axis 2)
with(dat,plot(y=PC1_morpho,x=PC2,xlab="Genetic principal component 2",ylab="Morphometric principal component 1",pch=symbols[HybClass_LDjacea],bg=bg[HybClass_LDjacea],col=colors[HybClass_LDjacea],ylim=c(-3,3.5),yaxp=c(-3,3,6)))
text(-4.5,3.5,"(B)")

# morophometric PC-scores (axis 2) ~ genomic PC-scores (axis 1)
with(dat,plot(y=PC2_morpho,x=PC1,xlab="Genetic principal component 1",ylab="Morphometric principal component 2",pch=symbols[HybClass_LDjacea],bg=bg[HybClass_LDjacea],col=colors[HybClass_LDjacea]))
text(-6,2.15,"(C)")
text(6.5,2.15,"chisq (1) = 5.5, p < 0.05")  # korrekt (180622)
lines(seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1),predict(lmer3,newdata=data.frame(PC1=seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1)),re.form=NA))

# morophometric PC-scores (axis 2) ~ genomic PC-scores (axis 2)
with(dat,plot(y=PC2_morpho,x=PC2,xlab="Genetic principal component 2",ylab="Morphometric principal component 2",pch=symbols[HybClass_LDjacea],bg=bg[HybClass_LDjacea],col=colors[HybClass_LDjacea]))
text(-4.5,2.2,"(D)")


dev.off()




