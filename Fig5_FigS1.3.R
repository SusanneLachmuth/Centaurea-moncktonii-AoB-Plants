# Code for Figures 5, Figure S1.3 and (if applicable) corresponding statistical analyses for publication
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
dat <- read.table("PCA_data_CnigraPanel.txt", header=T, sep = "\t") 
summary(dat)
# Order hybrid class levels:
dat$HybClass_LDnigra<- factor(dat$HybClass_LDnigra,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC1_nigra","BC2_nigra","nigra"))


###########################################################
### Analyze relationships between morophometric PC-scores and hybrid classes as well as genomic PC-scores:

# morophometric PC-scores (axis 1) ~ Hybrid Class (based on C. nigra panel)
lmer1<-lmer(PC1_morpho~HybClass_LDnigra+(1|pop),data=dat,REML=F) 
dropterm(lmer1,test="Chisq") 
#plot(lmer1) 
#hist(resid(lmer1))
posthoc1<-glht(lmer1, linfct = mcp(HybClass_LDnigra = "Tukey"))
cld(posthoc1)

# morophometric PC-scores (axis 2) ~ Hybrid Class (based on C. nigra panel)
lmer2<-lmer(PC2_morpho~HybClass_LDnigra+(1|pop),data=dat,REML=F) 
dropterm(lmer2,test="Chisq")
#plot(lmer2) 
#hist(resid(lmer2)) 
posthoc2<-glht(lmer2, linfct = mcp(HybClass_LDnigra = "Tukey"))
cld(posthoc2)

# morophometric PC-scores (axis 1) ~ genomic PC-scores (axis 1)
lmer3<-lmer(PC1_morpho~PC1+(1|pop),data=dat,REML=F) 
dropterm(lmer3,test="Chisq") 
#plot(lmer3) 

# morophometric PC-scores (axis 2) ~ genomic PC-scores (axis 1)
lmer3a<-lmer(PC2_morpho~PC1+(1|pop),data=dat,REML=F) 
dropterm(lmer3a,test="Chisq")
#plot(lmer3a) 

# morophometric PC-scores (axis 1) ~ genomic PC-scores (axis 2)
lmer4<-lmer(PC1_morpho~PC2+(1|pop),data=dat,REML=F) 
dropterm(lmer4,test="Chisq") 
#plot(lmer4)

# morophometric PC-scores (axis 2) ~ genomic PC-scores (axis 2)
lmer4a<-lmer(PC2_morpho~PC2+(1|pop),data=dat,REML=F) 
dropterm(lmer4a,test="Chisq")
#plot(lmer4a) 


###########################################################
### Figure 5:

# Colors and symbols
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"
colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",2),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#FF9DA3","#AA343B","#640000")    
symbols<- c(15,23,22,19,24,25,17)                                          


# Plot
jpeg("Fig5.jpg",width=17,height=17,units="cm",res= 600)

par(mfrow=c(2,2),mar=c(5,4,1,2),cex.axis=0.9,cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1)

####################################################################
#Boxplots (.~ Hybrid Class)

# morophometric PC-scores (axis 1)
with(dat,plot(y=PC1_morpho,x=HybClass_LDnigra,col=bg,xlab="Hybrid Class (NewHybrids)",ylab="Morphometric principal component 1",names=c("C.cf.jac","BC2jac","BC1jac","F2","BC1nig","BC2 nig","C.cf.nig"),ylim=c(-3,4),yaxp=c(-3,3,6)))
means1<-tapply(dat$PC1_morpho,dat$HybClass_LDnigra,mean)
points(means1,pch=21,bg="grey20",col="white",cex=1.5)
text(0.5,4,"(A)",cex=1.3)
text(seq(1,7,1),rep(3.5,7),c("b","ab","ab","ab","ab","ab","a"))   
text(6,4,"chisq (1) = 13.1, p < 0.05")   
legend(x=0.25,y=3,col=bg, pch=19,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

# morophometric PC-scores (axis 2)
with(dat,plot(y=PC2_morpho,x=HybClass_LDnigra,xlab="Hybrid Class (NewHybrids)",ylab="Morphometric principal component 2",pch=19,col=bg,names=c("C.cf.jac","BC2jac","BC1jac","F2","BC1nig","BC2nig","C.cf.nig"),ylim=c(-3,4),yaxp=c(-3,2,5)))
means2<-tapply(dat$PC2_morpho,dat$HybClass_LDnigra,mean)
points(means2,pch=21,bg="grey20",col="white",cex=1.5)
text(0.5,4,"(B)",cex=1.3)   
text(seq(1,7,1),rep(3,7),c("a","ab","ab","ab","b","ab","ab"))    
text(6,4,"chisq (1) = 14.8, p < 0.05")  

####################################################################
# Sactterplots (.~ genetic PC1)

# Colors and symbols
colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",2),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#FF9DA3","#AA343B","#640000")    
symbols<- c(15,23,22,19,24,25,17)                                          


# morophometric PC-scores (axis 1)
with(dat,plot(y=PC1_morpho,x=PC1,xlab="Genetic principal component 1",ylab="Morphometric principal component 1",pch=symbols[dat$HybClass_LDnigra],bg=bg[dat$HybClass_LDnigra],col=colors[dat$HybClass_LDnigra],ylim=c(-3,3.5),yaxp=c(-3,3,6)))
lines(seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1),predict(lmer3,newdata=data.frame(PC1=seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1)),re.form=NA))
text(-6.1,3.5,"(C)",cex=1.3)
text(6,3.5,"chisq (1) = 10.6, p < 0.01")  #korrekt (180622)
legend(-6.5,3.2,col=colors,pt.bg =bg, pch=symbols,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), #
       cex=0.75, bty = "n", title = "Hybrid Class")

# morophometric PC-scores (axis 2)
with(dat,plot(y=PC2_morpho,x=PC1,xlab="Genetic principal component 1",ylab="Morphometric principal component 2",pch=symbols[dat$HybClass_LDnigra],bg=bg[dat$HybClass_LDnigra],col=colors[dat$HybClass_LDnigra]))
text(-6.1,2.15,"(D)",cex=1.3)
text(6,2.15,"chisq (1) = 5.5, p < 0.05")  #korrekt (180622)
lines(seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1),predict(lmer3a,newdata=data.frame(PC1=seq(min(dat$PC1,na.rm=T),max(dat$PC1,na.rm=T),0.1)),re.form=NA))


dev.off()


#######################################################################################################################################
#######################################################################################################################################

####################################################################
### FigS1.3


tiff("FigS1.3.tiff",width=3.5,height=7,units="in",res=800)

par(mfrow=c(2,1),mar=c(5,4,1,2),cex.axis=1, cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1) 

# morophometric PC-scores (axis 1)
with(dat,plot(y=PC1_morpho,x=PC2,xlab="Genetic principal component 2",ylab="Morphometric principal component 1",pch=symbols[HybClass_LDnigra],bg=bg[HybClass_LDnigra],col=colors[HybClass_LDnigra],ylim=c(-3,4.5),yaxp=c(-3,3,6)))
text(-4.7,4.5,"(A)")
legend(-4.8,4.3,col=colors,pt.bg =bg, pch=symbols,legend=c("C.cfjacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

# morophometric PC-scores (axis 2)
with(dat,plot(y=PC2_morpho,x=PC2,xlab="Genetic principal component 2",ylab="Morphometric principal component 2",pch=symbols[HybClass_LDnigra],bg=bg[HybClass_LDnigra],col=colors[HybClass_LDnigra]))
text(-4.7,2.2,"(B)")

dev.off()

