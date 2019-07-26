# Code for Figure S2.3 and (if applicable) corresponding statistical analyses for publication
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
dat$HybClass_LDjacea<-factor(dat$HybClass_LDjacea,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC2_nigra","nigra"))


###########################################################
### Analyze relationships between morophometric PC-scores and hybrid classes:

# morophometric PC-scores (axis 1) ~ Hybrid Class (based on C. jacea panel)
lmer1<-lmer(PC1_morpho~HybClass_LDjacea+(1|pop),data=dat,REML=F) 
dropterm(lmer1,test="Chisq") 
#plot(lmer1) 
#hist(resid(lmer1)) 
posthoc1<-glht(lmer1, linfct = mcp(HybClass_LDjacea = "Tukey"))
cld(posthoc1)  

# morophometric PC-scores (axis 2) ~ Hybrid Class (based on C. jacea panel)
lmer2<-lmer(PC2_morpho~HybClass_LDjacea+(1|pop),data=dat,REML=F) 
dropterm(lmer2,test="Chisq")
#plot(lmer2) 
#hist(resid(lmer2))
posthoc2<-glht(lmer2, linfct = mcp(HybClass_LDjacea = "Tukey"))
cld(posthoc2)


###########################################################
### Figure 5:

# Colors and symbols
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"

colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",1),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#AA343B","#640000")    
symbols<- c(15,23,22,19,25,17)     

# File
tiff("FigS2.3.tiff",width=3.5,height=7,units="in",res=800)

# Plot
par(mfrow=c(2,1),mar=c(5,4,1,2),cex.axis=0.9,cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1)

#####################################################################
#Boxplots (.~ Hybrid Class)

# morophometric PC-scores (axis 1)
with(dat,plot(y=PC1_morpho,x=HybClass_LDjacea,col=bg,ylab="Morphometric principal component 1",names=c("C.cf.jac","BC2jac","BC1jac","F2","BC2nig","C.cf.nig"),ylim=c(-3,6),yaxp=c(-3,3,6)))
means1<-tapply(dat$PC1_morpho,dat$HybClass_LDjacea,mean)
points(means1,pch=21,bg="grey20",col="white",cex=1.5)
text(0.5,6.1,"(A)")
text(seq(1,7,1),rep(3.3,7),c("b","a","ab","a","ab","a")) 
text(5.3,6.1,"chisq (1) = 20.8, p < 0.001")
legend(x=0.3,y=5.9,col=bg, pch=19,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")


# morophometric PC-scores (axis 2)
with(dat,plot(y=PC2_morpho,x=HybClass_LDjacea, xlab= "Hybrid Class (NewHybrids)",ylab="Morphometric principal component 2",pch=19,col=bg,names=c("C.cf.jac","BC2jac","BC1jac","F2","BC2nig","C.cf.nig"),ylim=c(-3,3),yaxp=c(-3,2,5)))
means2<-tapply(dat$PC2_morpho,dat$HybClass_LDjacea,mean)
points(means2,pch=21,bg="grey20",col="white",cex=1.5)
text(0.5,3,"(B)")
text(5.3,3,"chisq (1) = 3.3, p > 0.05")


dev.off()




