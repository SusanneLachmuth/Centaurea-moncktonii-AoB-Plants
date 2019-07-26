# Code for Figures S2.5, S2.6 and (if applicable) corresponding statistical analyses for publication
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

gen_dat <- read.table("pop_genetic_data.txt", header=T, sep = "\t") 
summary(gen_dat)
# Order hybrid class levels:
gen_dat$HybClass_LDjacea<- factor(gen_dat$HybClass_LDjacea,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC2_nigra","nigra"))


###########################################################
### Analyze relationships between genome size, genomic as well as morphometric PC-scores and hybrid class:



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


# genome size ~ morphometric PC-scores (axis 1)
lmer3<-lmer(GenSize~PC1_morpho+(1|pop),data=dat,REML=F)  #[gen_cap_pca_mean$GenSize<2.05,]
dropterm(lmer3,test="Chisq")
#plot(lmer3) 
#hist(resid(lmer3)) 

# genome size ~ morphometric PC-scores (axis 2)
lmer4<-lmer(GenSize~PC2_morpho+(1|pop),data=dat,REML=F) 
dropterm(lmer4,test="Chisq")
#plot(lmer4) 
#hist(resid(lmer4)) 

# genome size ~ Hybrid Class (based on C. jacea panel)
lmer5<-lmer(GenSize~HybClass_LDjacea+(1|pop),data=gen_dat,REML=F) 
dropterm(lmer5,test="Chisq") 
#plot(lmer5) 
#hist(resid(lmer5)) 
posthoc5<-glht(lmer5, linfct = mcp(HybClass_LDjacea = "Tukey"))
cld(posthoc5)


###########################################################
### Figures:

# Colors and symbols
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"
colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",1),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#AA343B","#640000")    #,"#FF9DA3"
symbols<- c(15,23,22,19,25,17)      # ,24                              

###########################################################
### Figure S2.5


tiff("FigS2.5.tiff",width=3.5,height=10.5,units="in",res=800)

par(mfrow=c(3,1), mar=c(5,4,1,2), cex.axis=0.8, cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1)

gen_dat$HybClass_LDjacea<-factor(gen_dat$HybClass_LDjacea,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC2_nigra","nigra"))
with(gen_dat,plot(y=GenSize,x=HybClass_LDjacea,ylab="1C genome size (pg)",xlab="Hybrid Class (NewHybrids)",names=c("C.cf.jac","BC2jac","BC1jac","F2","BC2nig","C.cf.nig"),pch=19,col=bg,ylim=c(1.85,2.22),yaxp=c(1.85,2.15,6))) 
means3<-tapply(gen_dat$GenSize,gen_dat$HybClass_LDjacea,mean,na.rm=T)
points(means3,pch=21,bg="grey20",col="white",cex=1.5)
text(0.55,2.22,"(A)",cex=1.3)
text(seq(1,7,1),rep(2.2,7),c("b","b","ab","b","ab","a")) #
text(5.3,2.22,"chisq (1) = 19.5, p < 0.01")
legend(0.5,2.16,col=colors, pch=symbols,pt.bg=bg,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

with(gen_dat,plot(y=GenSize,x=PC1,xlab="Genetic principal component 1",ylab="1C genome size (pg)",pch=symbols[HybClass_LDjacea],col=colors[HybClass_LDjacea],bg=bg[HybClass_LDjacea],ylim=c(1.85,2.22),yaxp=c(1.85,2.15,6)))
text(-6.6,2.22,"(B)",cex=1.3)
lines(seq(min(gen_dat$PC1,na.rm=T),max(gen_dat$PC1,na.rm=T),0.1),predict(lmer1,newdata=data.frame(PC1=seq(min(gen_dat$PC1,na.rm=T),max(gen_dat$PC1,na.rm=T),0.1)),re.form=NA))
text(7.5,2.22,"chisq (1) = 5.0, p < 0.05")

with(gen_dat,plot(y=GenSize,x=PC2,xlab="Genetic principal component 2",ylab="1C genome size (pg)",pch=symbols[HybClass_LDjacea],col=colors[HybClass_LDjacea],bg=bg[HybClass_LDjacea],ylim=c(1.85,2.22),yaxp=c(1.85,2.15,6)))
text(-4.7,2.22,"(C)",cex=1.3)

dev.off()


###########################################################
### Figure S2.6

tiff("FigS2.6.tiff",width=3.5,height=7,units="in",res=800)

par(mfrow=c(2,1), mar=c(5,4,1,2), cex.axis=0.8, cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1)

with(dat,plot(y=GenSize,x=PC1_morpho,xlab="Morphometric principal component 1",ylab="1C genome size (pg)",pch=symbols[HybClass_LDjacea],col=colors[HybClass_LDjacea],bg=bg[HybClass_LDjacea],ylim=c(1.85,2.22),yaxp=c(1.85,2.15,6))) 
text(-2.65,2.22,"(A)",cex=1.3)
legend(-2.75,2.19,col=colors, pch=symbols,pt.bg=bg,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

with(dat,plot(y=GenSize,x=PC2_morpho,xlab="Morphometric principal component 2",ylab="1C genome size (pg)",pch=symbols[HybClass_LDjacea],col=colors[HybClass_LDjacea],bg=bg[HybClass_LDjacea],ylim=c(1.85,2.22),yaxp=c(1.85,2.15,6))) 
text(-2.1,2.22,"(B)",cex=1.3)

dev.off()
