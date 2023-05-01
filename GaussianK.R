## subfunction that computes the Gaussian kernel
## X is a distance data

d2p<- function(X,perplexity = 15,tol = 1e-5){
    D=as.matrix(X)
    n=dim(D)[1]
    P = matrix(0,n,n)
    beta= rep(1,n)
    logU=log(perplexity)
    
    for (i in 1:n){
      # set minimum and maximum values for precision
      betamin = -Inf
      betamax = Inf
      Di = D[i, -i]
      
      # Compute the Gaussian kernel and entropy
      hbeta = .Hbeta(Di, beta[i])
      H = hbeta$H; 
      thisP = hbeta$P
      Hdiff = H - logU;
      tries = 0;
      
      # binary search for sigma_i^2 for each i
      while(abs(Hdiff) > tol && tries < 2000){
        if (Hdiff > 0){
          betamin = beta[i]
          if (is.infinite(betamax)) beta[i] = beta[i] * 2
          else beta[i] = (beta[i] + betamax)/2
        } else{
          betamax = beta[i]
          if (is.infinite(betamin))  beta[i] = beta[i]/ 2
          else beta[i] = ( beta[i] + betamin) / 2
        }
        
        hbeta = .Hbeta(Di, beta[i])
        H = hbeta$H
        thisP = hbeta$P
        Hdiff = H - logU
        tries = tries + 1
      }	
      P[i,-i]  = thisP	
    }	
    
    r = NULL
    r$P = P
    
    r
}

