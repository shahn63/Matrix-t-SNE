dim(hedenfalk.expr)
head(hedenfalk.expr)

# re-arrange the data according to mutation type
hedenfalk.expr<-hedenfalk.expr[,c(1:6,18,7:10,19:22,11:17)]
dat<-t(hedenfalk.expr)

dim(dat)
dat<-dat[-c(16:22),] #remove the Sporadic type
grp<-c(rep(1,7),rep(2,8)) # membership

# t-test
zz<-pval<-NULL
for(i in 1:dim(dat)[2]){
  pval[i]<-anova(lm(dat[,i]~as.factor(grp)))[1,5]
  zz[i]<-summary(lm(dat[,i]~as.factor(grp)))$coef[2,3]
}

library(locfdr)
lfdr<-locfdr(zz)

sel.dat<-dat[,which(lfdr$fdr<0.2)]

## patient(row)-wise distance

n<-dim(sel.dat)[1]
distE<-matrix(0,n,n)

for(i in 1:(n-1)){
  for(j in (i+1):n){
    distE[i,j]<-stats::dist(rbind(sel.dat[i,],sel.dat[j,]), method = "euclidean")
  }
}
distE

distE[lower.tri(distE)]<-t(distE)[lower.tri(distE)]

X1<-distE

## gene(column)-wise distance

sel.dat<-t(sel.dat)
n<-dim(sel.dat)[1]
distE<-matrix(0,n,n)

for(i in 1:(n-1)){
  for(j in (i+1):n){
    distE[i,j]<-stats::dist(rbind(sel.dat[i,],sel.dat[j,]), method = "euclidean")
  }
}
dim(distE)

distE[lower.tri(distE)]<-t(distE)[lower.tri(distE)]

X2<-distE