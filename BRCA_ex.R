## Example: BRCA data
## by Hedenfalk et al. (2001). 'Gene-expression profiles in hereditary breast caner. The New Engliand Journal od Medicine.
## In this study, we used most important 176 genes in distinguishing BRCA1 and BRCA2
## The 176 genes' list are available in Hedenfalk et al. (2001)

## X1: squared distance between patients (mutation types)
## X2: squared distance between genes

source("Matrix_tSNE.R")

# Matrix t-SNE with alpha=0
set.seed(20230701)
Score0<-Matrix_tsne(X1,X2, alpha=0, k=2,perplexity1=5,perplexity2=15)
plot(Score0[,,1],Score0[,,2])

# kmeans clustering for column-wise group membership
Score<-Score0

set.seed(20230701)
wss<-ratio<-NULL
for(i in 1:10){
  km<-kmeans(cbind(c(Score[,,1]),c(Score[,,2])), centers = i, iter.max = 10000)
  wss[i]<-km$tot.withinss
  ratio[i]<-1-km$tot.withinss/km$totss
}
plot(ratio)

num.cl<-4

km<-kmeans(cbind(c(Score[,,1]),c(Score[,,2])), centers = num.cl, iter.max = 10000)

plot(c(Score[,,1]),c(Score[,,2]),col=km$cluster,cex=0.9, pch=km$cluster,lwd=1,
     main="Genes (alpha=0)",xlab="Score 1",ylab="Score 2")

grp0<-km$cluster
idx<-unique(grp0)

grp0[km$cluster==idx[1]]<-1
grp0[km$cluster==idx[2]]<-2
grp0[km$cluster==idx[3]]<-3
grp0[km$cluster==idx[4]]<-4


# BRCA1: 1/ BRCA2: 2
col.type=rep(c(rep(1,7),rep(2,8)),dim(Score)[2])
pch.type=rep(c(rep(1,7),rep(4,8)),dim(Score)[2])
cex.type=rep(c(rep(1.8,7),rep(1,8)),dim(Score)[2])

plot(c(Score[,,1]),c(Score[,,2]),col=col.type,pch=pch.type,cex=cex.type,lwd=1,
     main=paste("Mutation (alpha=0)"),xlab="Score1",ylab="Score 2")


# Matrix t-SNE with alpha=1
set.seed(20230701)
Score1<-Matrix_tsne(X1,X2, alpha=1, k=2,perplexity1=5,perplexity2=15,max_iter=1500)
plot(Score1[,,1],Score1[,,2])

Score<-Score1

plot(c(Score[,,1]),c(Score[,,2]),col=col.type,pch=pch.type,cex=1,lwd=1, #xlim=c(-40,80),
     main=paste("Mutation (alpha=0.86)"),xlab="Score1",ylab="Score 2")
legend("topright",c("BRCA1","BRCA2"),col=1:2,lty=0,cex=1,lwd=2,pch=1,box.col=0)

plot(c(Score[,,1]),c(Score[,,2]),col=grp0+4,cex=(grp0+6)/5, pch=grp0+4,lwd=1, #xlim=c(-40,80),
     main="Genes (alpha=0.86)",xlab="Score 1",ylab="Score 2")
legend("topright",c("Gene group 1","Gene group 2","Gene group 3","Gene group 4"),
       col=5:8,lty=0,cex=1.1,lwd=2,pch=5:8,box.col=0)

# Visualization with selected alpha=0.86

set.seed(20230701)
Score86<-Matrix_tsne(X1,X2, alpha=0.86, k=2,perplexity1=5,perplexity2=15)

Score<-Score86

plot(c(Score[,,1]),c(Score[,,2]),col=col.type,pch=pch.type,cex=1,lwd=1, #xlim=c(-40,80),
     main=paste("Mutation (alpha=0.86)"),xlab="Score1",ylab="Score 2")
legend("topright",c("BRCA1","BRCA2"),col=1:2,lty=0,cex=1,lwd=2,pch=1,box.col=0)

plot(c(Score[,,1]),c(Score[,,2]),col=grp0+4,cex=(grp0+6)/10, pch=grp0+4,lwd=1, #xlim=c(-40,80),
     main="Genes (alpha=0.86)",xlab="Score 1",ylab="Score 2")
legend("topright",c("Gene group 1","Gene group 2","Gene group 3","Gene group 4"),
       col=5:8,lty=0,cex=1.1,lwd=2,pch=5:8,box.col=0)

