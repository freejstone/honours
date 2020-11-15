# n must be greater than p to begin with for row extension to take place
row_extension = function(X, y, n_simul_ko, est_sigma){
  n <- dim(X)[1]
  p <- dim(X)[2]
  k <- n_simul_ko + 1
  if (n < k*p){
    if (est_sigma <= 0){
      Decom <- svd(X, nu = n, nv = 0)
      U_2 <- Decom$u[, c((p + 1): n)]
      sigma <- sqrt(mean((t(U_2)%*%y)^2)) #sqrt(RSS/(n-p)), (fourier expansion of colspace(U) + pythag will get you there)
    } else {
      sigma <- est_sigma #use est_sigma
    }
    y_extra <- randn(k*p - n, 1)*sigma
    y <- rbind(y, y_extra)
    X <- rbind(X, zeros(k*p-n, p))
  }
  return(list(X, y))
}