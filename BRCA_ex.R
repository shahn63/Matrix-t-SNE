## Example: BRCA data
## X1: distance between patients (mutation types)
## X2: distance between genes

# Matrix t-SNE with alpha=0
ydata0<-tsne(X1,X2, alpha=0, k=2,perplexity1=5,perplexity2=5)

# kmeans clustering for column-wise group membership
ydata<-ydata0

wss<-ratio<-NULL
for(i in 1:10){
  km<-kmeans(cbind(c(ydata[,,1]),c(ydata[,,2])), centers = i, iter.max = 10000)
  wss[i]<-km$tot.withinss
  ratio[i]<-1-km$tot.withinss/km$totss
}

num.cl<-which(ratio>0.8)[1]

km<-kmeans(cbind(c(ydata[,,1]),c(ydata[,,2])), centers = num.cl, iter.max = 10000)

plot(c(ydata[,,1]),c(ydata[,,2]),col=km$cluster,cex=0.9, pch=km$cluster,lwd=1,
     main="Genes (alpha=0)",xlab="Score 1",ylab="Score 2")

grp0<-km$cluster
idx<-unique(grp0)

grp0[km$cluster==idx[1]]<-1
grp0[km$cluster==idx[2]]<-2
grp0[km$cluster==idx[3]]<-3


col.type=rep(c(rep(1,7),rep(2,8)),dim(ydata)[2])
pch.type=rep(c(rep(1,7),rep(4,8)),dim(ydata)[2])
cex.type=rep(c(rep(1.8,7),rep(1,8)),dim(ydata)[2])

plot(c(ydata[,,1]),c(ydata[,,2]),col=col.type,pch=pch.type,cex=cex.type,lwd=1,
     main=paste("Mutation (alpha=0)"),xlab="Score1",ylab="Score 2")


# Matrix t-SNE with alpha=1

ydata1<-tsne(X1,X2, alpha=1, k=2,perplexity1=5,perplexity2=5)

ydata<-ydata1

wss<-ratio<-NULL
for(i in 1:10){
  km<-kmeans(cbind(c(ydata[,,1]),c(ydata[,,2])), centers = i, iter.max = 10000)
  wss[i]<-km$tot.withinss
  ratio[i]<-1-km$tot.withinss/km$totss
}
plot(wss,type="b")
plot(ratio,type="b")

num.cl<-which(ratio>0.8)[1]

km<-kmeans(cbind(c(ydata[,,1]),c(ydata[,,2])), centers = num.cl, iter.max = 10000)

grp1<-km$cluster

plot(c(ydata[,,1]),c(ydata[,,2]),col=grp1,cex=0.9, pch=grp1,lwd=1,
     main="Genes (alpha=0)",xlab="Score 1",ylab="Score 2")


# Visualization with selected alpha=0.83

ydata83<-tsne(X1,X2, alpha=0.83, k=2,perplexity1=5,perplexity2=5)

ydata<-ydata83

plot(c(ydata[,,1]),c(ydata[,,2]),col=col.type,pch=pch.type,cex=1,lwd=1,xlim=c(-40,80),
     main=paste("Mutation (alpha=0.83)"),xlab="Score1",ylab="Score 2")
legend("topright",c("BRCA1","BRCA2"),col=1:2,lty=0,cex=1,lwd=2,pch=1,box.col=0)

plot(c(ydata[,,1]),c(ydata[,,2]),col=grp0+4,cex=(grp0+6)/5, pch=grp0+4,lwd=1,xlim=c(-40,80),
     main="Genes (alpha=0.83)",xlab="Score 1",ylab="Score 2")
legend("topright",c("Gene group 1","Gene group 2","Gene group 3"),
       col=5:7,lty=0,cex=1.1,lwd=2,pch=5:7,box.col=0)

