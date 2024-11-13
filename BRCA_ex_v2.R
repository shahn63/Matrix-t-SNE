## Example: BRCA data
## by Hedenfalk et al. (2001). 'Gene-expression profiles in hereditary breast caner. The New Engliand Journal od Medicine.
## In this study, we used most important 176 genes in distinguishing BRCA1 and BRCA2
## The 176 genes' list are available in Hedenfalk et al. (2001)

## X1: squared distance between patients (mutation types)
## X2: squared distance between genes
## Load 'BRCA176.RData' which contains squared distance matrices X1 and X2

load("BRCA176.RData")

source("Matrix_tSNE_v2.R")
library(Rtsne)

##################
## row-wise t-SNE (Matrix t-SNE with alpha = 1) using Rtsne
##################

set.seed(901)
row.res<-Rtsne(sqrt(X1),pca=FALSE,perplexity=3,is_distance=T,theta=0,max_iter=1000)
ratio=NULL
for(i in 1:8){
	km<-kmeans(row.res$Y,centers = i, iter.max = 10000)
	ratio[i]<-1-km$tot.withinss/km$totss
}

## select the number of subject clusters (mutation)

G1<-which(ratio>0.9)[which.min(ratio>0.9)]

## generate and re-arrange group memberships vector
set.seed(901)
km<-kmeans(row.res$Y,centers = G1, iter.max = 10000)
grp1=km$cl
idx<-unique(grp1)
for(i in 1:G1){
	grp1[km$cluster==idx[i]]<-i
}

##################
## column-wise t-SNE (Matrix t-SNE with alpha = 0)
##################

set.seed(901)
col.res<-Rtsne(sqrt(X2),pca=FALSE,perplexity=30,is_distance=T,theta=0,max_iter=1000)
ratio=NULL
for(i in 1:8){
	set.seed(901)
	km<-kmeans(col.res$Y,centers = i, iter.max = 10000)
	ratio[i]<-1-km$tot.withinss/km$totss
}

## select the number of subject clusters (mutation)

G0<-which(ratio>0.9)[which.min(ratio>0.9)]

## generate and re-arrange group memberships vector

set.seed(901)
km<-kmeans(col.res$Y,centers = G0, iter.max = 10000)
grp0=km$cl
idx<-unique(grp0)
for(i in 1:G0){
	grp0[km$cluster==idx[i]]<-i
}

## set the number of clusters for Matrix t-SNE

G=G0*G1

## search the optimal alpha using equation (7)

ratio=NULL
aa=seq(0.02,0.98,0.02)
for(ii in 1:length(aa)){
	set.seed(901)
	Score<-Matrix_tsne(X1,X2, alpha=aa[ii], k=2, perplexity1=3, perplexity2=30)
	km<-kmeans(cbind(c(Score[,,1]),c(Score[,,2])),centers = G, iter.max = 10000)
	ratio[ii]<-1-km$tot.withinss/km$totss
}
opt.alpha=aa[which.max(ratio)]

##################################
# Matrix t-SNE with optimal alpha
##################################

set.seed(901)
Score<-Matrix_tsne(X1,X2, alpha=opt.alpha, k=2,perplexity1=3, perplexity2=30)

# BRCA1: 1/ BRCA2: 2
col.type=rep(c(rep(1,7),rep(2,8)),dim(Score)[2])
pch.type=rep(c(rep(1,7),rep(4,8)),dim(Score)[2])
cex.type=rep(c(rep(1.8,7),rep(1,8)),dim(Score)[2])
vec.grp0=rep(grp0,each=15)

par(mfrow=c(3,1))

## Matrix t-SNE result
## 1. embedding plot
plot(c(Score[,,1]),c(Score[,,2]),col=1,pch=1,cex=1,lwd=1, ylim=c(-80,80),
     main=paste("Matrix t-SNE (alpha=0.88)"),xlab="Score 1",ylab="Score 2")

## 2. plot with coloring with mutation type (alpha=0)
plot(c(Score[,,1]),c(Score[,,2]),col=col.type,pch=pch.type,cex=1,lwd=1, ylim=c(-80,80),
     main=paste("Mutation"),xlab="Score1",ylab="Score 2")
legend("bottomright",c("BRCA1","BRCA2"),col=1:2,lty=0,cex=1,lwd=1,pch=c(1,4),box.col=0)

## 3. plot with coloring with gene group (alpha=1)
plot(c(Score[,,1]),c(Score[,,2]),col=vec.grp0+3,cex=0.8, pch=vec.grp0+3,lwd=1, ylim=c(-80,80),
     main="Gene",xlab="Score 1",ylab="Score 2")
legend("bottomright",c("Gene group 1","Gene group 2","Gene group 3","Gene group 4","Gene group 5"),
       col=c(4:8),lty=0,cex=1.05,lwd=1,pch=4:8,box.col=0)






