# Code for Figure 6 and (if applicable) corresponding statistical analyses for publication
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


### Read data:
setwd("D:/data")
gen_dat<-read.table("pop_genetic_data.txt",header = T, sep="\t")
summary(gen_dat)
# Order hybrid class levels
gen_dat$HybClass_LDnigra<-factor(gen_dat$HybClass_LDnigra,levels=c("jacea","BC2_jacea","BC1_jacea","F2","BC1_nigra","BC2_nigra","nigra"))


###########################################################
### Analyze relationships between genome size and hybrid classes as well as genomic PC-scores:

# genome size ~ Hybrid Class (based on C. nigra panel)
summary(gen_dat)
lmer1<-lmer(GenSize~HybClass_LDnigra+(1|pop),data=gen_dat,REML=F) 
dropterm(lmer1,test="Chisq") 
#plot(lmer1)
#hist(resid(lmer1))
posthoc1<-glht(lmer1, linfct = mcp(HybClass_LDnigra = "Tukey"))
cld(posthoc1)

# genome size ~ genomic PC-scores (axis 1)
lmer2<-lmer(GenSize~PC1+(1|pop),data=gen_dat,REML=T,na.action=na.omit) #[gen_dat$GenSize<2.05,]
dropterm(lmer2,test="Chisq") #  0.02487 * (180625)
summary(lmer2)
#plot(lmer2)
#hist(resid(lmer2)) 


###########################################################
### Figure 6:

# Colors and symbols
rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"
colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",2),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#FF9DA3","#AA343B","#640000")    
symbols<- c(15,23,22,19,24,25,17)                                          


### Plot:

jpeg("Fig6.jpg",width=8,height=16,units="cm",res= 600)

par(mfrow=c(2,1),mar=c(5,4,1,2),cex.axis=0.8, cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1) 

with(gen_dat,plot(y=GenSize,x=HybClass_LDnigra,ylab="1C genome size (pg)",xlab="Hybrid Class (NewHybrids)",pch=19,col=bg,names=c("C.cf.jac","BC2jac","BC1jac","F2","BC1nig","BC2nig","C.cf.nig"),ylim=c(1.85,2.22),yaxp=c(1.85,2.15,6))) 
means3<-tapply(gen_dat$GenSize,gen_dat$HybClass_LDnigra,mean,na.rm=T)
points(means3,pch=21,bg="grey20",col="white",cex=1.5)
text(0.5,2.22,"(A)",cex=1.3)
text(seq(1,7,1),rep(2.2,7),c("b","b","b","b","ab","ab","a")) # 
text(6,2.22,"chisq (1) = 18.26, p < 0.01")
legend(x=0.5,y=2.18,col=bg, pch=19,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")


with(gen_dat,plot(y=GenSize,x=PC1,xlab="Genetic principal component 1",ylab="1C genome size (pg)",pch=symbols[HybClass_LDnigra],col=colors[HybClass_LDnigra],bg=bg[HybClass_LDnigra],ylim=c(1.85,2.2),yaxp=c(1.85,2.15,6))) 
text(-6.5,2.2,"(B)",cex=1.3)
lines(seq(min(gen_dat$PC1,na.rm=T),max(gen_dat$PC1,na.rm=T),0.1),predict(lmer2,newdata=data.frame(PC1=seq(min(gen_dat$PC1,na.rm=T),max(gen_dat$PC1,na.rm=T),0.1)),re.form=NA))
text(7,2.2,"chisq (1) = 5.03, p < 0.05")
legend(-6.7,2.18,col=colors, pch=symbols,pt.bg=bg,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC1nigra","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")

dev.off()



