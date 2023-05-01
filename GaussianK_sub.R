## subfunction that computes the Gaussian kernel 
## beta = 1/ 2*sigma^2
## D: distance matrix
## H = H(perplexity)
## P = P_j|i

.Hbeta <-
  function(D, beta){
    P = exp(-D * beta)
    sumP = sum(P)
    if (sumP == 0){
      H = 0
      P = D * 0
    } else {
      H = log(sumP) + beta * sum(D %*% P) /sumP
      P = P/sumP
    }
    r = NULL
    r$H = H
    r$P = P
    r
  }
