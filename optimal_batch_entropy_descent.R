optimal_batch_entropy_descent <- function(I, Sig, nsko){
  lenI <- length(I)
  p <- ncol(Sig)
  I_IP <- matrix(0, nrow = lenI, ncol = p)
  I_IP[(I - 1)*lenI + (1:lenI)] <- 1
  I_PI <- t(I_IP)
  
  Df <- function(s0){
    if (length(s0) == 1){
      Df <- diag(solve(((nsko + 1)/nsko)*Sig - s0*I_PI%*%I_IP))
    } else{
      Df <- diag(solve(((nsko + 1)/nsko)*Sig - I_PI%*%diag(s0)%*%I_IP))
    }
    return(Df[I] - nsko/s0)
  }
  

  mat <- function(s0){
    if (length(s0) == 1){
      return(((nsko + 1)/nsko)*Sig - s0*I_PI%*%I_IP)
    } else{
      return(((nsko + 1)/nsko)*Sig - I_PI%*%diag(s0)%*%I_IP)
    }
  }
  

  s0 <- rep(0.000001, lenI)
  count <- 0
  if (lenI < 20){
    factor <- 1
  } else {
    factor <- lenI/20
  }
  
  start = 0.0001*factor

  while (count < 120000){
    asdf = Df(s0)
    length = norm(as.matrix(asdf), type = "1")
    s1 <- s0 - start*asdf/length
    if( all(s1 < 1) & all(s1 > 0) & min(eigen(mat(s1))$values) > 0){
      s0 <- s1
      count <- count + 1
    } else{
      break
    }
  }
  print(Df(s0))
  return(s0)
}
