optimal_batch_SDP_G <- function(I, Sig, nsko){
  lenI <- length(I)
  p <- ncol(Sig)
  s1 <- Variable(lenI)
  I_IP <- matrix(0, nrow = lenI, ncol = p)
  I_IP[(I - 1)*lenI + (1:lenI)] <- 1
  I_PI <- t(I_IP)
  
  objective <- Minimize(sum(abs(1- s1)))
  
  
  # #Setting up constraints
  constraints = list(s1 >= 0, lambda_min(((nsko + 1)/nsko)*Sig - I_PI%*%diag(s1)%*%I_IP) >= 0)
  
  problem <- Problem(objective, constraints = constraints)
  result <- solve(problem, eps = 1e-08, max_iters = 20000)
  s0_optimal <- result$getValue(s1)
  return(s0_optimal)
}