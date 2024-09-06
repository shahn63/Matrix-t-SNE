## this is the code for two given distances as inputs
## two distances for a matrix data: row-wise squared distance and column-wise squared distance
## e.g.) Euclidean distance for general data/ DTW for trajectory data
## alpha is a weight
## perplexity1 for row (perp^r)
## perplexity2 for column (perp^c)

.libPaths("C:/R-4.1.1/library")
setwd("D:/matrix_tSNE/")
require(Rcpp)
require(RcppArmadillo)
sourceCpp("jointprob.cpp")
sourceCpp("gradient.cpp")
source("GaussianK.R")
source("GaussianK_sub.R")


Matrix_tsne<- function(X1, X2, alpha, initial_config = NULL, k=2, perplexity1=5, perplexity2=5, max_iter = 1000, min_cost=0,epoch=100){
  n1<-nrow(X1) # X1: the row-wise squared distance
  n2<-nrow(X2) # X2: the column-wise squared distance
  
  momentum = .5 # initial momentum
  final_momentum = .8 # value to which momentum is changed
  mom_switch_iter = 250 # iteration at which momentum is changed
  
  epsilon = 1000 # initial learning rate
  min_gain = .01 # minimum gain for delta-bar-delta
  
  eps = 2^(-52) # typical machine precision
  Score<-array(0,dim=c(n1,n2,k)) # low-dimensional scores
  
  # set initial scores
  # set the same random initial score over columns if alpha=1
  # set the same random initial score over rows if alpha=0
  # set the random initial score, otherwise
  if(alpha==1){
    new<-rnorm(k*n1)
    for(i in 1:k){
      Score[,,i]<-rep(new[(n1*(i-1)+1):(n1*i)],n2)
    }
  }else if(alpha==0){
    new<-rnorm(k*n2)
    for(i in 1:k){
      Score[,,i]<-rep(new[(n2*(i-1)+1):(n2*i)],each=n1)
    }
  }else{
    Score1=Score2=array(0,dim=c(n1,n2,k))
    new<-rnorm(k*n1)
    for(i in 1:k){
      Score1[,,i]<-rep(new[(n1*(i-1)+1):(n1*i)],n2)
    }
    new<-rnorm(k*n2)
    for(i in 1:k){
      Score2[,,i]<-rep(new[(n2*(i-1)+1):(n2*i)],each=n1)
    }
    Score=Score1+Score2
  }
  
  # compute joint probabilities
  # make sure p-values are set properly
  
  P1=d2p(X1, perplexity1, 1e-5)$P
  P2=d2p(X2, perplexity2, 1e-5)$P
  
  P1 = .5 * (P1 + t(P1))
  P2 = .5 * (P2 + t(P2))
  
  P1[P1 < eps]<-eps; 	P2[P2 < eps]<-eps
  P1 = P1/sum(P1); P2 = P2/sum(P2)
  P1 = P1 * 4; P2 = P2 * 4
  
  grads =  array(0,dim=c(n1,n2,k))
  grads1 =  array(0,dim=c(n1,n2,k)); gradsc1 = matrix(0, nrow=n1,ncol=n1+n2)
  grads2 =  array(0,dim=c(n1,n2,k)); gradsc2 = matrix(0, nrow=n1,ncol=n1+n2)
  incs =   array(1,dim=c(n1,n2,k))
  gains = array(0,dim=c(n1,n2,k))
  sum_Score1=matrix(0,nrow=k,ncol=n1)
  sum_Score2=matrix(0,nrow=k,ncol=n2)
  
  num1=matrix(0,n1,n1); num2=matrix(0,n2,n2)
  
  for (iter in 1:max_iter){
    if (iter %% epoch == 0) { # epoch
      cost =  sum(P1 * log((P1+eps)/(Q1+eps)))+sum(P2 * log((P2+eps)/(Q2+eps)))
      message("Epoch: Iteration #",iter," error is: ",cost)
      if (cost < min_cost) break
    }
    
    num1=matrix(0,n1,n1); num2=matrix(0,n2,n2)
    
    # compute joint probablity that point i and j are neighbors
    
    num1<-jointprob(Score[,,1],Score[,,2], n1, n2)
    num1 = 1/(1+num1) 

    num2<-jointprob(t(Score[,,1]),t(Score[,,2]), n2, n1)
    num2 = 1/(1+num2) 
   
    diag(num1)=0; diag(num2)=0 # set diagonal to zero
    Q1 = num1 / sum(num1); Q2 = num2 / sum(num2) # normalize to get probabilities
    if (any(is.nan(num1))) message ('NaN in grad. descent')
    if (any(is.nan(num2))) message ('NaN in grad. descent')
    Q1[Q1 < eps] = eps; Q2[Q2 < eps] = eps
    
    # compute the gradients
    stiffnesses1 = 4 * alpha^2 * (P1-Q1) * num1 
    stiffnesses2 = 4 * (1-alpha)^2 * (P2-Q2) * num2 
    
    gradsc1<-gradient(Score[,,1],Score[,,2], stiffnesses1, n1, n2)
    grads1[,,1]<-gradsc1[,1:n2]
    grads1[,,2]<-gradsc1[,-c(1:n2)]

    gradsc2<-gradient(t(Score[,,1]),t(Score[,,2]), stiffnesses2, n2, n1)
    grads2[,,1]<-t(gradsc2[,1:n1])
    grads2[,,2]<-t(gradsc2[,-c(1:n1)])

    grads= grads1+grads2
    
    # update the solution
    gains = ((gains + .2) * abs(sign(grads) != sign(incs)) +
               gains * .8 * abs(sign(grads) == sign(incs)))
    
    gains[gains < min_gain] = min_gain
    
    incs =   momentum * incs - epsilon * (gains * grads)
    
    Score = Score + incs
    
    if (iter == mom_switch_iter) momentum = final_momentum
    
    if (iter == 100 && is.null(initial_config)){P1 = P1/4; P2 = P2/4}
  }
  Score # output (low-dimensional scores)
}
