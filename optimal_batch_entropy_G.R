optimal_batch_entropy_G <- function(I, Sig, nsko){
  lenI <- length(I)
  p <- ncol(Sig)
  s1 <- Variable(lenI)
  I_IP <- matrix(0, nrow = lenI, ncol = p)
  I_IP[(I - 1)*lenI + (1:lenI)] <- 1
  I_PI <- t(I_IP)
  
  #Objective function
  #Had to present the objective function in a peculiar way so that CVXR accepts this as a convex function
  if (lenI == 1){
    objective <- Minimize(-log_det( ((nsko + 1)/nsko)*Sig - s1*I_PI%*%I_IP) - nsko*log(s1))
  } else {
    objective <- Minimize(-log_det( ((nsko + 1)/nsko)*Sig - I_PI%*%diag(s1)%*%I_IP) - nsko*sum(log(s1)))
  }
  
  
  #Setting up constraints
  constraints = list(s1 >= 0, lambda_min(((nsko + 1)/nsko)*Sig - I_PI%*%diag(s1)%*%I_IP) >= 0)
  
  problem <- Problem(objective, constraints = constraints, eps = 1e-08, max_iters = 20000)
  result <- solve(problem, solver = "SCS")
  s0_optimal <- result$getValue(s1)
  return(s0_optimal)
}