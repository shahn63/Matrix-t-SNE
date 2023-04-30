## this is the code for two given distances as inputs
## two distances for a matrix-variate data: row-wise distance and column-wise distance
## e.g.) Euclidean distance for general data/ DTW for trajectory data
## alpha is a weight
## perplexity1 for row
## perplexity2 for column

Matrix_tsne<-
  function(X1, X2, alpha, initial_config = NULL, k=2, perplexity1=5, perplexity2=5, max_iter = 1000, min_cost=0,epoch=100){
    n1<-nrow(X1) # X1: the row-wise distance
    n2<-nrow(X2) # X2: the column-wise distance
    
    momentum = .5 # initial momentum
    final_momentum = .8 # value to which momentum is changed
    mom_switch_iter = 250 # iteration at which momentum is changed
    
    epsilon = 1000 # initial learning rate
    min_gain = .01 # minimum gain for delta-bar-delta
    
    eps = 2^(-52) # typical machine precision
    ydata<-array(0,dim=c(n1,n2,k)) # low-dimensional scores
    
    # set initial scores
    # set the same random initial score over columns if alpha=1
    # set the same random initial score over rows if alpha=0
    # set the random initial score, otherwise
    if(alpha==1){
      new<-rnorm(k*n1)
      for(i in 1:k){
        ydata[,,i]<-rep(new[(n1*(i-1)+1):(n1*i)],n2)
      }
    }else if(alpha==0){
      new<-rnorm(k*n2)
      for(i in 1:k){
        ydata[,,i]<-rep(new[(n2*(i-1)+1):(n2*i)],each=n1)
      }
    }else{
      ydata1=ydata2=array(0,dim=c(n1,n2,k))
      new<-rnorm(k*n1)
      for(i in 1:k){
        ydata1[,,i]<-rep(new[(n1*(i-1)+1):(n1*i)],n2)
      }
      new<-rnorm(k*n2)
      for(i in 1:k){
        ydata2[,,i]<-rep(new[(n2*(i-1)+1):(n2*i)],each=n1)
      }
      ydata=ydata1+ydata2
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
    grads1 =  array(0,dim=c(n1,n2,k))
    grads2 =  array(0,dim=c(n1,n2,k))
    incs =   array(1,dim=c(n1,n2,k))
    gains = array(0,dim=c(n1,n2,k))
    sum_ydata1=matrix(0,nrow=k,ncol=n1)
    sum_ydata2=matrix(0,nrow=k,ncol=n2)
    
    num1=matrix(0,n1,n1); num2=matrix(0,n2,n2)
    
    for (iter in 1:max_iter){
      if (iter %% epoch == 0) { # epoch
        cost =  sum(P1 * log((P1+eps)/(Q1+eps)))
        message("Epoch: Iteration #",iter," error is: ",cost)
        if (cost < min_cost) break
      }
      
      num1=matrix(0,n1,n1); num2=matrix(0,n2,n2)
      
      # compute joint probablity that point i and j are neighbors
      
      for(jj in 1:n2){
        sk2 = apply(ydata[,jj,]^2, 1, sum)
        num1 = num1+sk2+sweep(-2 * ydata[,jj,] %*% t(ydata[,jj,]),2, -t(sk2))
      }
      num1 = 1/(1+num1) 
      
      for(ii in 1:n1){
        sk2 = apply(ydata[ii,,]^2, 1, sum)
        num2 = num2+sk2+sweep(-2 * ydata[ii,,] %*% t(ydata[ii,,]),2, -t(sk2))
      }
      num2 = 1/(1+num2) 
      
      diag(num1)=0; diag(num2)=0 # set diagonal to zero
      Q1 = num1 / sum(num1); Q2 = num2 / sum(num2) # normalize to get probabilities
      if (any(is.nan(num1))) message ('NaN in grad. descent')
      if (any(is.nan(num2))) message ('NaN in grad. descent')
      Q1[Q1 < eps] = eps; Q2[Q2 < eps] = eps
      
      # compute the gradients
      stiffnesses1 = 4 * alpha^2 * (P1-Q1) * num1 
      stiffnesses2 = 4 * (1-alpha)^2 * (P2-Q2) * num2 
      
      for (jj in 1:n2){
        for(ii in 1:n1){
          grads1[ii,jj,] = apply(sweep(-ydata[,jj,], 2, -ydata[ii,jj,]) * stiffnesses1[,ii],2,sum)
        }
      }
      for (ii in 1:n1){
        for(jj in 1:n2){
          grads2[ii,jj,] = apply(sweep(-ydata[ii,,], 2, -ydata[ii,jj,]) * stiffnesses2[jj,],2,sum)
        }
      }
      grads= grads1+grads2
      
      # update the solution
      gains = ((gains + .2) * abs(sign(grads) != sign(incs)) +
                 gains * .8 * abs(sign(grads) == sign(incs)))
      
      gains[gains < min_gain] = min_gain
      
      incs =   momentum * incs - epsilon * (gains * grads)
      
      ydata = ydata + incs

      if (iter == mom_switch_iter) momentum = final_momentum
      
      if (iter == 100 && is.null(initial_config)){P1 = P1/4; P2 = P2/4}
    }
    ydata # output (low-dimensional scores)
  }