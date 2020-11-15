create_candidate_batch_G <- function(I, Sig, nsko, s0){
  p <- ncol(Sig)
  lenI <- length(I)
  lenG <- nsko*lenI + p
  G <- matrix(0, lenG, lenG)
  P <- 1:p
  if (length(s0) == 1){
    diag_s <- diag(rep(s0, p))
    Sig0 <- Sig - diag_s
  } else if (length(s0) > 1){
    diag_s <- diag(c(s0))
    I_IP <- diag(1, nrow = lenI, ncol = p)
    I_IP[(I - 1)*lenI + (1:lenI)] <- 1
    I_PI <- t(I_IP)
    Sig0 <- Sig - I_PI%*%diag_s%*%I_IP
  }
  
  SigII <- Sig[I,I]
  Sig0II <- Sig0[I,I]
  Sig0IP <- Sig0[I,P]
  Sig0PI <- Sig0[P,I]
  
  G[P, P] <- Sig
  for (i in 1:nsko){
    G[P , (p+(i-1)*lenI+1) : (p+i*lenI)] <- Sig0PI
    G[(p+(i-1)*lenI+1) : (p+i*lenI), P] <- Sig0IP
  }
  for (i in 1:nsko){
    G[(p+(i-1)*lenI+1) : (p+i*lenI), (p+(i-1)*lenI+1) : (p+i*lenI)] <- SigII
    for (j in 1:nsko){
      if (j != i){
        G[(p+(i-1)*lenI+1) : (p+i*lenI) , (p+(j-1)*lenI+1) : (p+j*lenI)] <- Sig0II
      }
    }
  }
  return(G)
}