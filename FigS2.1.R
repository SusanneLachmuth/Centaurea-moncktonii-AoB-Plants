# Code for Figure S2.1 and (if applicable) corresponding statistical analyses for publication
# Lachmuth S, Molofsky J, Suda J, Milbrath L, Keller SR (accepted for AoB Plants on 17 June 2019) 
# Associations between genomic ancestry, genome size and capitula morphology in the invasive meadow knapweed hybrid complex 
# (Centaurea ×moncktonii C.E. Britton) in eastern North America

### Clean working directory:

rm(list=ls())
gc()


### Load packages:

library(data.table)
library(vcfR)
library(adegenet)


### Read genomic data and connect them with population information:

setwd("D:/data")

# SNP data
vcf1 <- read.vcfR("ind_SNP_data.vcf")
gl1 <- vcfR2genlight(vcf1)
str(gl1)

dat <- read.table("pop_genetic_data.txt",header=T,sep="\t") 
summary(dat)

popID<-dat[,1:2]
summary(popID)

popID1<-droplevels(cbind(gl1$ind.names,popID[match(gl1$ind.names,popID$Ind),]))
summary(popID1)
gl1$pop <- popID1$pop # assign locality info



##########################################################
# Perform PCA:

pca1 <- glPca(gl1, nf=4, parallel=T) 
print(pca1)


##########################################################
# Plot (run code directly in R, e.g. RStudio confuses the population symbols):

tiff("FigS2.1.tiff",width=7,height=7,units="in",res=800)

par(mfrow=c(2,2),cex.axis=1,cex.lab=1.1,mgp=c(2.5,0.5,0),cex=0.6,lwd=1,mar=c(5,4,1,2))

### Plot 3A: PCA plot with grey symbols labeling populations

popcode1 <-popID1$pop

plot(pca1$scores[,1], pca1$scores[,2], 
     pch=c(1:18,21,22)[gl1$pop], col="grey20",bg="grey"
     ,ann=F 
)
mtext(side = 1, text = "Genetic principal component 1", line = 2.5,cex=0.7)
mtext(side = 2, text = "Genetic principal component 2", line = 2.5,cex=0.7) 

legend(8,6.6, cex=0.75, legend=levels(as.factor(popcode1)), pch=c(1:18,21,22), 
       col="black", pt.bg = "grey", pt.cex=1,bty = "n", title = "Population")
text(-6.5,6.4,"(A)")


### Plot S2.1B: PCA plot with colored symbols labeling hybrid classes

# Define hybrid classes
head(dat)
levels(dat$HybClass_LDjacea)
dat$HybClass_LDjacea<-factor(dat$HybClass_LDjacea,levels = c("jacea","BC2_jacea","BC1_jacea","F2","BC2_nigra","nigra"))
levels(dat$HybClass_LDjacea)
names(dat)
HybClass<-dat[,c(2,11)]

hybID1<-droplevels(cbind(gl1$ind.names,HybClass[match(gl1$ind.names,HybClass$Ind),]))
head(hybID1)
# assign HybClass info (has to be named pop to be compatible with genlight)
gl1$pop <- factor(hybID1$HybClass_LDjacea,levels=c("jacea","BC2_jacea","BC1_jacea","F2", "BC2_nigra", "nigra")) 
gl1$pop[1:3]
hybCode1 <-factor(hybID1$HybClass_LDjacea,levels=c("jacea","BC2_jacea","BC1_jacea","F2", "BC2_nigra", "nigra"))
levels(gl1$pop)

rgb(red=191,green=129,blue=45,names="Dbeige",maxColorValue=255)
# dark beige: "#BF812D"

colors<-c("#0000E9",rep("grey20",2),"#BF812D",rep("grey20",1),"#640000") 
bg <- c("#0000E9","#6682EC","#A6C0FF","#BF812D","#AA343B","#640000")     
symbols<- c(15,23,22,19,25,17)                                          



plot(pca1$scores[,1], pca1$scores[,2], 
     pch=symbols[gl1$pop],col=colors[gl1$pop], bg=bg[gl1$pop]
     ,ann=F 
     #,ylab="" 
     #,main="",xlim=c(-10,15)
)
mtext(side = 1, text = "Genetic principal component 1", line = 2.5,cex=0.7)
mtext(side = 2, text = "Genetic principal component 2", line = 2.5,cex=0.7)


legend(6,6.6,col=colors,pt.bg =bg, pch=symbols,legend=c("C.cf.jacea","BC2jacea","BC1jacea","F2","BC2nigra","C.cf.nigra"), 
       cex=0.75, bty = "n", title = "Hybrid Class")
text(-6.5,6.4,"(B)")


### Plot S2.1C-D: Relate Admixture and hybrid classes to PC scores

# (C) PC1 and C. nigra cluster (K=2)

plot(dat$PC1,dat$AnK2_K1, ylab= "",xlab=""
     , ann=F
     ,pch=symbols[dat$HybClass_LDjacea],col=colors[dat$HybClass_LDjacea], bg=bg[dat$HybClass_LDjacea]
) 
mtext(side = 1, text = "Genetic principal component 1", line = 2.5,cex=0.7)
mtext(side = 2, text = "Assignment probability to C.cf. nigra cluster", line = 2.5,cex=0.7)

text(-6.5,1,"(C)")


# (D) PC1 and hybrid cluster (K=3)

plot(dat$PC2,dat$AnK3_K3, ylab= "",xlab=""
     , ann=F
     ,pch=symbols[dat$HybClass_LDjacea],col=colors[dat$HybClass_LDjacea], bg=bg[dat$HybClass_LDjacea]
) 
mtext(side = 1, text = "Genetic principal component 2", line = 2.5,cex=0.7)
mtext(side = 2, text = "Assignment probability to hybrid cluster", line = 2.5,cex=0.7)

text(-4.75,1,"(D)")


dev.off()
