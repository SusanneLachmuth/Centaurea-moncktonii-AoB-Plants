# Code for Figure S1.1 and (if applicable) corresponding statistical analyses for publication
# Lachmuth S, Molofsky J, Suda J, Milbrath L, Keller SR (accepted for AoB Plants on 17 June 2019) 
# Associations between genomic ancestry, genome size and capitula morphology in the invasive meadow knapweed hybrid complex 
# (Centaurea ×moncktonii C.E. Britton) in eastern North America

### Clean working directory:

rm(list=ls())
gc()


### Load packages:
library(stats)
library(MASS)
library(lme4)
library(car)


### Read and format data:
setwd("D:/data")

### Read and check data:
cap <- read.table("capitula_data.txt", header=T, sep = "\t") 
summary


###########################################################
### Select morphometric traits that differ significantly between C. cf. jacea and C. cf. nigra

# Data subset
cap_parent<-droplevels(cap[cap$taxon!="hybrid1"&cap$taxon!="hybrid2",])
summary(cap_parent)

# Nr. of bract rows
glmer1<-glmer(Brac_RowNr~taxon+(1|Pop/IndID),data=cap_parent,family = "poisson")
dropterm(glmer1,test="Chisq") # Using lmer (without any transformation) the effect is significant

# % pectinate bract rows
lmer2<-lmer(Perc_SerBrac~taxon+(1|Pop/IndID),data=cap_parent,REML=F)
dropterm(lmer2,test="Chisq")
#plot(lmer2)

# % pectinate bract rows (logit transformed)
lmer2a<-lmer(logit(Perc_SerBrac/100)~taxon+(1|Pop/IndID),data=cap_parent,REML=F)
dropterm(lmer2a,test="Chisq")
#plot(lmer2a)

# bract color
lmer3<-lmer(Brac_Col~taxon+(1|Pop/IndID),data=cap_parent,REML=F)
dropterm(lmer3,test="Chisq")
#plot(lmer3) 

# Ratio capitulum width to length
lmer4<-lmer(Cap_Ratio~taxon+(1|Pop/IndID),data=cap_parent,REML=F)
dropterm(lmer4,test="Chisq")
#plot(lmer4) 

# Ratio appendage center width to length
lmer5<-lmer(ApCen_Ratio~taxon+(1|Pop/IndID),data=cap_parent,REML=F)
dropterm(lmer5,test="Chisq")
#plot(lmer5)

# Rel. width of bract appendage center
lmer6<-lmer(relWid_BracApCen~taxon+(1|Pop/IndID),data=cap_parent,REML=F)
dropterm(lmer6,test="Chisq")
#plot(lmer6) 


###########################################################
### Figure S1.1:

# Colors
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"
rgb(red=223,green=194,blue=125,names="Lbeige",maxColorValue=255)
# light beige: "#DFC27D" 

par(mfrow=c(3,2),cex.main=1.2,mar=c(3.5,4.3,2.2,0.5),xpd=NA,cex.axis=1.2,cex.lab=1.2)

with(cap,boxplot(Brac_RowNr~taxon,main="Nr. of bract rows",names=c ("C. cf. jacea","C. cf. nigra","Hybrids NY", "Hybrids VT"), col= c("#0000E9","#640000","#BF812D","#DFC27D"),cex=1.5))
text(1.5,11.5,"n.sig")
lines(c(1.0,1.25,2.0),rep(11,3),lwd=1.8)
means<-tapply(cap$Brac_RowNr,cap$taxon,mean)
points(means,pch=21,bg="grey20",col="white",cex=1.5)

with(cap,boxplot(Perc_SerBrac~taxon,main="Pectinate bract rows (%)",names=c ("C. cf. jacea","C. cf. nigra","Hybrids NY", "Hybrids VT"), col= c("#0000E9","#640000","#BF812D","#DFC27D"),ylim=c(0,130),yaxp=c(0,100,5),cex=1.5))
text(1.5,120,"***",cex=1.5)
lines(c(1.0,1.25,2.0),rep(110,3),lwd=1.8)
means<-tapply(cap$Perc_SerBrac,cap$taxon,mean)
points(means,pch=21,bg="grey20",col="white",cex=1.5)

with(cap,boxplot(Brac_Col~taxon,main="Bract color",names=c ("C. cf. jacea","C. cf. nigra","Hybrids NY", "Hybrids VT"), col= c("#0000E9","#640000","#BF812D","#DFC27D"),ylim=c(1,3.8),yaxp=c(1,3,2),cex=1.5))
text(1.5,3.6,"n.sig")
lines(c(1.0,1.25,2.0),rep(3.4,3),lwd=1.8)
means<-tapply(cap$Brac_Col,cap$taxon,mean)
points(means,pch=21,bg="grey20",col="white",cex=1.5)

with(cap,boxplot(Cap_Ratio~taxon,main="Ratio capitulum width to length",names=c ("C. cf. jacea","C. cf. nigra","Hybrids NY", "Hybrids VT"), col= c("#0000E9","#640000","#BF812D","#DFC27D"),ylim=c(0.5,1.5),cex=1.5))
text(1.5,1.45,"n.sig")
lines(c(1.0,1.25,2.0),rep(1.35,3),lwd=1.8)
means<-tapply(cap$Cap_Ratio,cap$taxon,mean)
points(means,pch=21,bg="grey20",col="white",cex=1.5)

with(cap,boxplot(ApCen_Ratio~taxon,main="Ratio appendage center width to length",names=c ("C. cf. jacea","C. cf. nigra","Hybrids NY", "Hybrids VT"), col= c("#0000E9","#640000","#BF812D","#DFC27D"),cex=1.5))
text(1.5,1.6,"*",cex=1.5)
lines(c(1.0,1.25,2.0),rep(1.45,3),lwd=1.8)
means<-tapply(cap$ApCen_Ratio,cap$taxon,mean)
points(means,pch=21,bg="grey20",col="white",cex=1.5)

with(cap,boxplot(relWid_BracApCen~taxon,main="Rel. width of bract appendage center",names=c ("C. cf. jacea","C. cf. nigra","Hybrids NY", "Hybrids VT"), col= c("#0000E9","#640000","#BF812D","#DFC27D"),cex=1.5))
text(1.5,1,"**",cex=1.5)
lines(c(1.0,1.25,2.0),rep(0.93,3),lwd=1.8)
means<-tapply(cap$relWid_BracApCen,cap$taxon,mean)
points(means,pch=21,bg="grey20",col="white",cex=1.5)



