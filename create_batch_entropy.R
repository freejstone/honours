create_batch_entropy <- function(X, I, nsko, Qall, diagRall){
  n <- nrow(X)
  p <- ncol(X)
  Sig <- t(X)%*%X
  s0 <- optimal_batch_entropy_G(I, Sig, nsko)
  G <- create_candidate_batch_G(I, Sig, nsko, s0)
  lenI <- length(I)
  dimG <- nsko*lenI + p #p + d|I|
  Decom <- svd(G, nu = n, nv = 0)
  Decom$d[abs(Decom$d) < 10**(-12)] <- 0
  X0 <- Decom$u %*% diag(sqrt(Decom$d)) %*% t(Decom$u) #construction of X0 through SVD
  qr_X0 <- qr(X0)
  R0 <- qr.R(qr(X0))
  R0[, qr_X0$pivot] <- R0
  
  inds <- kronecker(matrix(1:nsko, nrow = 1), repmat(p, 1, lenI)) + kronecker(ones(1, nsko), matrix(I, nrow = 1))
  QI = cbind(Qall[,(1:p)], Qall[,c(inds)])
  diagRI = diagRall[1:p]
  
  #makes R0 agree with R so that first p of X1 is X!
  signs = cbind(matrix(sign(diagRI/diag(R0[1:p,1:p])), nrow = 1), ones(1, nsko * lenI))
  Qs = QI %*% diag(c(signs))
  X1 = Qs %*% R0
  X_ko = X1[, (p+1):dimG]
  if (any( abs(t(Xall)%*%Xall - G) > 10^(-7))){
    warning('OH NO: Gram matrix is different to G')
  }
  return(list(X_ko, s0))
}