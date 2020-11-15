#packages
library('knockoff', lib = '/users/stfour/jack/Documents')
library('glmnet', lib = '/users/stfour/jack/Documents')
library('pracma', lib = '/users/stfour/jack/Documents')
library('SuperPCA', lib = '/users/stfour/jack/Documents')
library('MASS', lib = '/users/stfour/jack/Documents')
library('scs', lib = '/users/stfour/jack/Documents')
library('clusterGeneration', lib = '/users/stfour/jack/Documents')
suppressWarnings(library('CVXR', lib = '/users/stfour/jack/Documents', warn.conflicts=FALSE))

#functions
scores_to_ranks = function(scores){
  n_scores <- numel(scores)
  rperm <- randperm(n_scores)
  scores1 <- scores[rperm]
  sort_inds <- sort(scores1, index.return = T)$ix
  ranks <- 1:n_scores
  ranks[sort_inds] <- 1:n_scores
  ranks[rperm] <- ranks
  ranks <- matrix(ranks, nrow = nrow(scores), ncol = ncol(scores))
  return(ranks)
}
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
    I_IP <- matrix(0, nrow = lenI, ncol = p)
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
    G[P, (p+(i-1)*lenI+1) : (p+i*lenI)] <- Sig0PI
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
  
  # p <- ncol(Sig)
  # lenI <- length(I)
  # lenG <- (nsko + 1)*p
  # G <- matrix(0, lenG, lenG)
  # P <- 1:p
  # J <- setdiff(P, I)
  # if (length(s0) == 1){
  #   diag_s <- diag(rep(s0, p))
  #   Sig0 <- Sig - diag_s
  # } else if (length(s0) > 1){
  #   diag_s <- diag(c(s0))
  #   I_IP <- matrix(0, nrow = lenI, ncol = p)
  #   I_IP[(I - 1)*lenI + (1:lenI)] <- 1
  #   I_PI <- t(I_IP)
  #   Sig0 <- Sig - I_PI%*%diag_s%*%I_IP
  # }
  # 
  # all_mats <- cbind(Sig, matrix(Sig0, nrow(Sig0), ncol(Sig0)*nsko))
  # G[P, ] <- all_mats
  # inds <- 1:ncol(G)
  # for (i in 1:nsko){
  #   inds <- c(inds[(nsko*p + 1): (p*(nsko + 1))], inds[1:(nsko*p)])
  #   G[(i*p + 1):((i + 1)*p), ] <- all_mats[, inds]
  # }
  # 
  # I = mod(I, p)
  # 
  # inds <- union(P, p + which((mod((p + 1):ncol(G), p) %in% I) == TRUE))
  return(G)
}
optimal_batch_G <- function(I, Sig, nsko){
  minEig <- function(s){
    mineig <- min(eigen(create_candidate_batch_G(I, Sig, nsko, s))$values)
    if (abs(mineig) < 1e-13){
      mineig <- 0
    }
    return(mineig)
  }
  sigeig <- eigen(Sig)$values
  s0 <- min(sigeig)
  tol <- 1e-5
  if (any(sigeig <= tol*max(sigeig))){
    warning('Model is not identifiable, but proceeding with knockoffs')
  }
  # analytical result says that s_equi for NO batching is sigeig*(nsko + 1)/nsko
  # so, batching has to be at least as good as this
  if (minEig(s0*(nsko + 1)/nsko) == 0){
    s_equi <- s0*(nsko + 1)/nsko
  } else {
    s_equi <- as.numeric(fzero(minEig, c(s0*(nsko + 1)/nsko, 2))$x)
  }
  G <- create_candidate_batch_G(I, Sig, nsko, min(1, s_equi))
  return(list(G, s_equi))
}
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
  
  
  # #Setting up constraints
  constraints = list(s1 >= 0, lambda_min(((nsko + 1)/nsko)*Sig - I_PI%*%diag(s1)%*%I_IP) >= 0)
  
  
  # #this takes a realllllly long time!!
  # #3.7 seconds for a p = 20, lenI = 2 problem
  problem <- Problem(objective, constraints = constraints)
  result <- solve(problem, solver = "SCS")
  s0_optimal <- result$getValue(s1)
  # s0_optimal[s0_optimal < 0] <- 0
  return(s0_optimal)
}
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
  
  
  # #this takes a realllllly long time!!
  # #3.7 seconds for a p = 20, lenI = 2 problem
  problem <- Problem(objective, constraints = constraints)
  result <- solve(problem, eps = 1e-08, max_iters = 20000)
  s0_optimal <- result$getValue(s1)
  # s0_optimal[s0_optimal < 0] <- 0
  return(s0_optimal)
}
create_batch_ko <- function(X, I, nsko, Qall, diagRall){
  n <- nrow(X)
  p <- ncol(X)
  Sig <- t(X)%*%X
  Results <- optimal_batch_G(I, Sig, nsko)
  G <- Results[[1]] #batch specific Gram
  s0 <- Results[[2]] #corresponding s0
  lenI <- length(I)
  dimG <- nsko*lenI + p #p + d|I|
  Decom <- svd(G, nu = nrow(G), nv = 0)
  X0 <- Decom$u %*% diag(sqrt(Decom$d)) %*% t(Decom$u) #construction of X0 through SVD
  qr_X0 <- qr(X0)
  Q0 <- qr.Q(qr_X0)
  R0 <- qr.R(qr_X0)
  R0[, qr_X0$pivot] <- R0
  
  
  inds <- kronecker(matrix(1:nsko, nrow = 1), repmat(p, 1, lenI)) + kronecker(ones(1, nsko), matrix(I, nrow = 1))
  QI = cbind(Qall[,(1:p)], Qall[,c(inds)])
  diagRI = diagRall[1:p]
  
  #makes R0 agree with R so that first p of X1 is X!
  signs = cbind(matrix(sign(diagRI/diag(R0[1:p,1:p])), nrow = 1), ones(1, nsko * lenI))
  Qs = QI %*% diag(c(signs))
  X1 = Qs %*% R0
  
  X_ko = X1[, (p+1):dimG]
  Xall <- cbind(X, X_ko)
  if (any( abs(t(Xall)%*%Xall - G) > 10^(-7))){
    print('OH NO')
  }
  return(list(X_ko, s0))
}
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
  return(list(X_ko, s0))
}
create_batch_SDP <- function(X, I, nsko, Qall, diagRall){
  n <- nrow(X)
  p <- ncol(X)
  Sig <- t(X)%*%X
  s0 <- optimal_batch_SDP_G(I, Sig, nsko)
  s0[s0 < 0] <- 0
  G <- create_candidate_batch_G(I, Sig, nsko, s0)
  lenI <- length(I)
  dimG <- nsko*lenI + p #p + d|I|
  Decom <- svd(G, nu = nrow(G), nv = 0)
  D <- Decom$d
  D[D < 0] <- 0
  X0 <- Decom$u %*% diag(sqrt(Decom$d)) %*% t(Decom$u) #construction of X0 through SVD
  qr_X0 <- qr(X0)
  Q0 <- qr.Q(qr_X0)
  R0 <- qr.R(qr_X0)
  R0[, qr_X0$pivot] <- R0
  
  
  inds <- kronecker(matrix(1:nsko, nrow = 1), repmat(p, 1, lenI)) + kronecker(ones(1, nsko), matrix(I, nrow = 1))
  QI = cbind(Qall[,(1:p)], Qall[,c(inds)])
  diagRI = diagRall[1:p]
  
  #makes R0 agree with R so that first p of X1 is X!
  signs = cbind(matrix(sign(diagRI/diag(R0[1:p,1:p])), nrow = 1), ones(1, nsko * lenI))
  Qs = QI %*% diag(c(signs))
  X1 = Qs %*% R0
  
  X_ko = X1[, (p+1):dimG]
  Xall <- cbind(X, X_ko)
  if (any( abs(t(Xall)%*%Xall - G) > 10^(-4))){
    print('OH NO')
  }
  return(list(X_ko, s0))
}
TDC_flex_c = function(decoy_wins, target_wins, BC1, thresholds){
  c <- thresholds@c
  lambda <- thresholds@lambda
  nTD <- as.matrix(cumsum(target_wins) + .Machine$double.xmin)
  nDD <- as.matrix(cumsum(decoy_wins))
  fdps <- as.matrix(pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals <- as.matrix(flipud(as.matrix(cummin(flipud(fdps)))))
  return(list(fdps, 'qvals', nTD, 'nDD'))
}
number_of_discoveries = function(obs_scores, decoys, method, alpha, printing, seed){
  set.seed(seed)
  obs_scores <- as.matrix(obs_scores)
  decoys <- as.matrix(decoys)
  m <- dim(decoys)[1]
  n_decoys <- dim(decoys)[2]
  if (method@type == 'random' || method@type == 'mirandom'){
    c <- method@threshold
  } else {
    c <- 0.5
  }
  if (method@lambda < c){
    if (method@lambda_lt_c){
      lambda <- c
    }
  } else {
    lambda <- method@lambda
  }
  mixed_scores <- -.Machine$double.xmax * ones(m, 1)
  
  sorted_decoys <- t(apply(decoys,1,sort))
  n_decoys_LT_obs <- as.matrix(rowSums(sweep(decoys, 1, obs_scores) < 0, 2))
  n_decoys_EQ_obs <- as.matrix(rowSums(sweep(decoys, 1, obs_scores) == 0, 2))
  n_decoys_GT_obs <- n_decoys - (n_decoys_LT_obs + n_decoys_EQ_obs)
  n_decoys_NE_obs <- n_decoys - n_decoys_EQ_obs
  p_obs <- (n_decoys_NE_obs - n_decoys_LT_obs + 1)/ (n_decoys_NE_obs + 1)
  obs_wins <- as.matrix(p_obs <= c)
  mixed_scores[which(obs_wins == 1)] <- obs_scores[which(obs_wins == 1)]
  
  if (method@type == 'mirandom'){
    #	Define the mapping of the losing observation ranks to the winning decoys ranks
    n_decoys_in_draw <- floor(c*(n_decoys + 1)) # the number of decoys the losing observations will be mapped to
    n_obs_ranks <- n_decoys - ceil((n_decoys + 1) * lambda) + 1	# the number of ranks in which the observation loses
    max_mapped_decoy_rank <- zeros(n_obs_ranks, 1) # the highest decoy rank each losing observation is mapped to
    max_mapped_decoy_prob <- zeros(n_obs_ranks, 1) # the probability of mapping to that maximal decoy rank
    current_decoy_rank <- n_decoys # the decoy rank currently being mapped into, starting from the highest rank
    current_decoy_coverage <- 0		# how much of that decoy rank did we cover so far
    uni_decoy_coverage <- n_obs_ranks / n_decoys_in_draw	# how much coverage should each mapped-to-decoy get
    for (i in 1 : n_obs_ranks){ # loop defining the map, one observed rank at a time
      max_mapped_decoy_rank[i] <- current_decoy_rank		# assign the current highest decoy rank available
      if (current_decoy_coverage + 1 > uni_decoy_coverage){ # mapping the current observation (i) to this decoy overfills its quota
        max_mapped_decoy_prob[i] <- uni_decoy_coverage - current_decoy_coverage	# so the probability is (quota) - (coverage so far)
        remainder_obs_coverage_prob <- 1 - max_mapped_decoy_prob[i]		# portion of the current observation we still need to map
        current_decoy_coverage <- uni_decoy_coverage;				# current decoy is saturated
      } else { # current observation can be fully mapped to current decoy
        max_mapped_decoy_prob[i] <- 1					# so the probability is 1
        current_decoy_coverage <- 1 + current_decoy_coverage	# update current decoy coverage
        remainder_obs_coverage_prob <- 0				# there is no "change", the current observation was fully mapped
      }
      if (current_decoy_coverage >= uni_decoy_coverage - 1e-12){ # is current decoy saturated (up to a roundoff error)?
        # floor(...) below is the number of additional decoy ranks the current observation can saturate
        current_decoy_rank <- current_decoy_rank - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) - 1
        current_decoy_coverage <- remainder_obs_coverage_prob - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) * uni_decoy_coverage	# and the remainder of the observation is allocated to the new decoy
      }
    }
    if (n_decoys_in_draw > 0 && (current_decoy_coverage > 1e-12 || n_decoys - current_decoy_rank != n_decoys_in_draw)){	# sanity check - increase 1e-12 for larger m
      print('Tell Uri the mapping doesnt work...')
    }
    # next, apply the mapping to every observation that lost
    decoy_wins <- (n_decoys_in_draw > 0) & (n_decoys_GT_obs > 0) & !obs_wins
    counted_decoy_wins <- decoy_wins & (p_obs > lambda)
    obs_ranks <- n_decoys_LT_obs[counted_decoy_wins] + floor((n_decoys_EQ_obs[counted_decoy_wins]+1)*rand(sum(counted_decoy_wins), 1))	# randomly broken ties
    obs_ranks <- pmin(obs_ranks, n_obs_ranks-1)		# because ties are randomly broken the initial rank can exceed n_obs_ranks-1
    # 	if c <= 0.5 mapped_obs_ranks = max_mapped_decoy_rank(1+obs_ranks) - (rand(sum(decoy_wins), 1) > max_mapped_decoy_prob(1+obs_ranks));
    rands = rand(sum(counted_decoy_wins), 1)
    mapped_obs_ranks = max_mapped_decoy_rank[1+obs_ranks] - ceil((rands - max_mapped_decoy_prob[1+obs_ranks]) / uni_decoy_coverage - 1e-12)*(rands > max_mapped_decoy_prob[1+obs_ranks])		# this is a bit tricky but it works (look at the notes to see why)
    mixed_scores[which(counted_decoy_wins == T)] <- sorted_decoys[(c(mapped_obs_ranks) - 1)*nrow(sorted_decoys) + which(counted_decoy_wins == T)] # the selected scores
  }
  
  if (method@randomTieBreak){
    set.seed(seed)
    tie_perm <- randperm(m)	# row vector
    mixed_scores <- mixed_scores[tie_perm]
    decoy_wins <- decoy_wins[tie_perm]
    obs_wins <- obs_wins[tie_perm]
    counted_decoy_wins <- obs_wins[tie_perm]
    p_obs <- p_obs[tie_perm]
  } else {
    tie_perm <- (1:m)
  }
  sort_idxs0 <- sort(mixed_scores, decreasing = T, index.return = T)$ix
  setClass(
    "thresholds",
    slots = list(
      lambda = "numeric",
      c = "numeric"
    )
  )
  if (method@lambda_lt_c){
    thresholds <- new("thresholds", lambda = method@lambda)
    counted_decoy_wins <- p_obs > method@lambda
  } else {
    thresholds <- new("thresholds", lambda = lambda)
    counted_decoy_wins <- decoy_wins & (p_obs > lambda)
  }
  thresholds@c <- c
  results <- TDC_flex_c(as.matrix(counted_decoy_wins[sort_idxs0]), as.matrix(obs_wins[sort_idxs0]), method@BC1, thresholds)
  
  fdps <- results[[1]]
  #qvals <- results[[2]]
  nTD <- results[[3]]
  nDD <- results[[4]]
  
  if (length(which(fdps <= alpha)) > 0 ){
    if (printing){
      print('we have some POSSIBLE discoveries')
    }
    i_acp <- max(which(fdps <= alpha))
    total_no <- nTD[i_acp]
    discovery_inds <- rep(FALSE, m)
    discovery_inds[sort_idxs0[1:i_acp]] <- TRUE
    discovery_inds <- discovery_inds & obs_wins
    if (sum(discovery_inds) + .Machine$double.xmin != total_no){
      stop('oh shit')
    }
    discovery_inds[tie_perm] <- discovery_inds
  } else {
    total_no <- 0
    discovery_inds <- logical(0)
  }
  
  # fdps[sort_idxs0] <- fdps
  # fdps[tie_perm] <- fdps
  # qvals[sort_idxs0] <- qvals
  # qvals[tie_perm] <- qvals
  # obs_wins[tie_perm] <- obs_wins
  # counted_decoy_wins[tie_perm] <- counted_decoy_wins
  # sort_idxs <- tie_perm[sort_idxs0]
  # p_obs[tie_perm] <- p_obs
  
  target_dist <- cumsum(as.matrix(obs_wins[sort_idxs0] & ((beta[tie_perm])[sort_idxs0] != 0)))
  decoy_dist <- cumsum(as.matrix(obs_wins[sort_idxs0] & ((beta[tie_perm])[sort_idxs0] == 0)))
  asdf <- cumsum(as.matrix((beta[tie_perm])[sort_idxs0] == 0))
  return(list(total_no, discovery_inds, fdps, target_dist, decoy_dist, asdf, obs_wins[sort_idxs0]))
}
lasso_scores = function(Xe, ye, n_simul_ko, max_n_simul_ko, shared_max_lambda, n_ko_batches, clustered_batches, tilXs, lenAugI, n_lam_factor, n_lam_const){
  n <- dim(Xe)[1]
  p <- dim(Xe)[2]
  Zs = zeros(p, 1)
  tilZs <- zeros(p, n_simul_ko)
  tilXdpY <- zeros(p, n_simul_ko)
  if (shared_max_lambda == T){
    max_lam <- max(abs(t(Xe)%*%ye)) #the largest possible coefficient of ye's projection on the columns of Xe
  }
  maxI <- 0
  for (ib in 1:n_ko_batches){
    if (length(clustered_batches) != 0){
      I <- which(clustered_batches == ib)
    } else {
      minI <- maxI + 1
      maxI <- minI + floor(p / n_ko_batches) + as.numeric(ib <= mod(p, n_ko_batches)) - 1
      I <- minI:maxI
    }
    X_ko <- tilXs[[ib, n_simul_ko]]
    tilXs_dp_y <- t(X_ko)%*%ye
    if (shared_max_lambda == T && max(abs(tilXs_dp_y)) > max_lam){
      max_lam <- max(abs(tilXs_dp_y))
    }
    for (ik in 1:n_simul_ko){
      tilXdpY[I, ik] <- t(tilXs_dp_y[((ik-1)*lenAugI[ib]+1) : ((ik-1)*lenAugI[ib] + lenAugI[ib])])
    }
  }
  maxI <- 0
  for (ib in 1:n_ko_batches){
    if (length(clustered_batches) != 0){
      I <- which(clustered_batches == ib)
    } else {
      minI <- maxI + 1
      maxI <- minI + floor(p / n_ko_batches) + as.numeric(ib <= mod(p, n_ko_batches)) - 1
      I <- minI:maxI
    }
    X_ko <- tilXs[[ib, n_simul_ko]]
    X_hat <- cbind(Xe, X_ko)
    #inds <- 1 : (n_simul_ko * lenAugI[ib] + p)
    if (n_lam_const){
      n_lambdas <- (1 + max_n_simul_ko) * p * n_lam_factor
    } else {
      n_lambdas <- (1 + n_simul_ko) * p * n_lam_factor
    }
    
    
    if (shared_max_lambda){
      #Changed nlambda to 500 for comparison sakes
      Z <- lassoMaxLambdaU(X_hat, ye, max_lam/n, 500)
    } else {
      #removed nlambda = n_lambdas for comparison
      
      # swap = rbinom(ncol(Xe), 1, 0.5)
      # swap.M = matrix(swap, nrow = nrow(Xe), ncol = length(swap), byrow = TRUE)
      # Xe.swap = Xe * (1 - swap.M) + X_ko * swap.M
      # X_ko.swap = Xe * swap.M + X_ko * (1 - swap.M)
      #rperm <- randperm(n_simul_ko * lenAugI[ib] + p)
      
      Z <- matrix(knockoff:::lasso_max_lambda(X_hat, ye, method = 'glmnet'), nrow = 1)
      #Z[rperm] <- Z
      # Z.swap <- Z[1:ncol(Xe)]
      # Zk.swap <- Z[(ncol(Xe) + 1): length(Z)]
      # Z <- Z.swap*(1 - swap) + Zk.swap*swap
      # Zk <- Z.swap*swap + Zk.swap*(1-swap)
      # Z <- cbind(Z, Zk)
    }
    Zs[I] <- t(Z[I])
    for (ik in 1:n_simul_ko){
      tilZs[I, ik] <- t(Z[p+(((ik-1)*lenAugI[ib]+1) : ((ik-1)*lenAugI[ib]+numel(I)))])
    }
  }
  return(list(Zs, tilZs, tilXdpY))
}


setClass(
  "mdt",
  slots = list(
    type = "character",
    BC1 = "numeric",
    poolDecoys = "logical",
    randomTieBreak = "logical",
    threshold = "numeric",
    lambda = "numeric",
    lambda_lt_c = "logical",
    all_lambdas = "numeric"
  )
)

'-------------------'

set.seed(1234)
n = 200                # number of observations
p = 50                # number of variables
k = 20                 # number of variables with nonzero coefficients
amplitude = 6         # signal amplitude (for noise level = 1)
n_simul_ko = 3
alpha_vec <- c(1/9, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4)
rho <- 0.3
n_runs = 1000
mu <- rep(0,p)
batch_vec = c(30, 40)

seeds <- sample(.Machine$integer.max, n_runs)


est_FDR_entropy <- logical(0)
est_FDR_equi <- logical(0)
est_FDR_SDP <- logical(0)

est_power_entropy <- logical(0)
est_power_equi <- logical(0)
est_power_SDP <- logical(0)

est_FDR_var_entropy <- logical(0)
est_FDR_var_equi <- logical(0)
est_FDR_var_SDP <- logical(0)

est_power_var_entropy <- logical(0)
est_power_var_equi <- logical(0)
est_power_var_SDP <- logical(0)

alpha_overall <- logical(0)
batch_overall <- logical(0)
c_overall <- logical(0)
lambda_overall <- logical(0)
batch <- logical(0)
counts_overall <- logical(0)
rho_overall <- logical(0)
noise_overall <- logical(0)

seeds <- sample(.Machine$integer.max, n_runs)

Sigma <- rcorrmatrix(p)

prop_null_target_wins_entropy <- logical(0)
prop_null_target_wins_SDP <- logical(0)

Xs <- matrix(list(), nrow = n_runs)
ys <- matrix(list(), nrow = n_runs)
betas_all <- matrix(list(), nrow = n_runs)

for (i_run in 1:n_runs){
  Xs[[i_run]] <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  betas <- amplitude * (-1)^randi(2, k, 1)
  inds <- sample(p, k)
  beta <- zeros(p,1)
  beta[inds] <- betas
  betas_all[[i_run]] <- beta
  ys[[i_run]] <- Xs[[i_run]]%*%beta + rnorm(n, sd = 1)
}

for (batch in batch_vec){
  
  X_ko_SDP <-  matrix(list(), nrow = max(batch_vec), ncol = n_runs)
  X_ko_entropy <- matrix(list() , nrow = max(batch_vec), ncol = n_runs)
  X_ko_equi <- matrix(list() , nrow = max(batch_vec), ncol = n_runs)
  
  
  s_SDP <- matrix(list(), nrow = max(batch_vec))
  s_entropy <- matrix(list(), nrow = max(batch_vec))
  s_equi <- matrix(list(), nrow = max(batch_vec))
  
  ranks_entropy <- matrix(list(), nrow = n_runs)
  ranks_SDP <- matrix(list(), nrow = n_runs)
  ranks_equi <- matrix(list(), nrow = n_runs)
  
  max_n_ko_batches = batch
  
  
  SigmaInv <- solve(Sigma)
  Decom <- svd(Sigma, nu = nrow(Sigma))
  Decom$d[abs(Decom$d) < 10**(-12)] <- 0
  X0 <- Decom$u %*% diag(sqrt(Decom$d)) %*% t(Decom$u)
  X <- X0[, 1:p]
  zeros <- matrix(0, nrow = n - p, ncol = p)
  X <- rbind(X, zeros)
  Xnorm <- normc(X)
  tree_batches <- hclust(dist(t(Xnorm), method = "euclidean"), method = "average")
  clustered_batches <- cutree(tree_batches, k = max_n_ko_batches)
  
  for (ib in 1:max_n_ko_batches){
    I <- which(clustered_batches == ib)
    s_SDP[[ib]] <- optimal_batch_SDP_G(I, Sigma, n_simul_ko)
    s_entropy[[ib]] <- optimal_batch_entropy_G(I, Sigma, n_simul_ko)
  }
  
  for (i_run in 1:n_runs){
    X <- Xs[[i_run]]
    for (ib in 1:max_n_ko_batches){
      I <- which(clustered_batches == ib)
      lenI <- length(I)
      s0 <- s_SDP[[ib]]
      G <- create_candidate_batch_G(I, Sigma, n_simul_ko, s0)
      Sig_conditional <- G[(p + 1):(p + n_simul_ko*lenI), (p + 1):(p + n_simul_ko*lenI)] - G[(p + 1):(p + n_simul_ko*lenI), 1:p]%*%SigmaInv%*%G[1:p, (p + 1):(p + n_simul_ko*lenI)]
      
      s <- svd(Sig_conditional)
      V <- s$v
      if (length(Sig_conditional) == 1){
        D <- s$d
        D[D < 0] <- 0
        D <- sqrt(D)
        X_ko = as.numeric(Re(V * D)) * t(matrix(rnorm(n*(n_simul_ko*lenI)), (n_simul_ko*lenI)))
      } else{
        D <-diag(c(s$d))
        D[D < 0] <- 0
        D <- sqrt(D)
        X_ko = t(Re(V %*% D) %*% matrix(rnorm(n*(n_simul_ko*lenI)), (n_simul_ko*lenI)))
      }
      X_ko_SDP[[ib, i_run]] <- X_ko
      for (i in 1:n){
        mu_conditional <- matrix(rep(mu[I], n_simul_ko), ncol = 1) + G[(p + 1):(p + n_simul_ko*lenI), 1:p]%*%SigmaInv%*%(X[i,] - matrix(mu, ncol = 1))
        X_ko_SDP[[ib, i_run]][i, ] <- X_ko_SDP[[ib, i_run]][i, ] + matrix(mu_conditional, nrow = 1)
      }
    }
  }
  
  for (i_run in 1:n_runs){
    y <- ys[[i_run]]
    X <- Xs[[i_run]]
    Zs = zeros(p, 1)
    tilZs <- zeros(p, n_simul_ko)
    for (ib in 1:max_n_ko_batches){
      I <- which(clustered_batches == ib)
      lenI <- length(I)
      Xall <- scale(cbind(Xs[[i_run]], X_ko_SDP[[ib, i_run]]))
      
      nlambda = 5*(1 + n_simul_ko)*p
      lambda_max = max(abs(t(Xall) %*% y))/n
      lambda_min = lambda_max/2000
      k = (0:(nlambda - 1))/nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
      
      inds <- sample(1:ncol(Xall))
      Z <- glmnet(Xall[, inds], y, family = "gaussian", intercept = F, lambda = lambda, standardize = F, thresh = 10^(-10), maxit = 10^7)
      first_nonzero <- function(x) match(T, abs(x) > 0)
      indices <- apply(Z$beta, 1, first_nonzero)
      names(indices) <- NULL
      Z = ifelse(is.na(indices), 0, Z$lambda[indices] * n)
      Z[inds] <- Z
      
      Zs[I] <- Z[I]
      for (ik in 1:n_simul_ko){
        tilZs[I, ik] <- t(Z[p+(((ik-1)*lenI+1) : ((ik-1)*lenI+lenI))])
      }
    }
    ranks_SDP[[i_run]] <- scores_to_ranks(cbind(Zs, tilZs))
  }
  
  
  for (alpha in alpha_vec){
    FDP <- rep(0, n_runs)
    power <- rep(0, n_runs)
    prop <- matrix(0, nrow = n_runs, ncol = p)
    for (i_run in 1:n_runs){
      beta = betas_all[[i_run]]
      ranks <- ranks_SDP[[i_run]]
      mdt <- new("mdt", type ='mirandom', BC1=1, poolDecoys = F, randomTieBreak = T, lambda_lt_c = F, threshold = 1/2, lambda = 1/2)
      results <- number_of_discoveries(ranks[,1], ranks[,2:ncol(ranks)], mdt, alpha, F, seeds[i_run])
      if (length(results[[2]]) == 0){
        FDP[i_run] <- 0
        power[i_run] <- 0
      } else {
        FDP[i_run] <- sum(results[[2]] & (beta == 0))/max(results[[1]], 1)
        discovery_list <- results[[2]] & (beta != 0)
        power[i_run] <- sum(discovery_list)/sum(beta != 0)
      }
      prop[i_run, ] <- results[[5]]/pmax(results[[6]], 1)
    }
    est_FDR_SDP <- c(est_FDR_SDP, mean(FDP))
    est_power_SDP <- c(est_power_SDP, mean(power))
    est_FDR_var_SDP <- c(est_FDR_var_SDP, var(FDP)/n)
    est_power_var_SDP <- c(est_power_var_SDP, var(power)/n)
    alpha_overall <- c(alpha_overall, alpha)
    batch_overall <- c(batch_overall, batch)
    if (alpha == 1/9){
      prop_null_target_wins_SDP <- c(prop_null_target_wins_SDP, apply(prop, 2, mean))
    }
  }
  
  for (i_run in 1:n_runs){
    X <- Xs[[i_run]]
    for (ib in 1:max_n_ko_batches){
      I <- which(clustered_batches == ib)
      lenI <- length(I)
      s0 <- s_entropy[[ib]]
      G <- create_candidate_batch_G(I, Sigma, n_simul_ko, s0)
      Sig_conditional <- G[(p + 1):(p + n_simul_ko*lenI), (p + 1):(p + n_simul_ko*lenI)] - G[(p + 1):(p + n_simul_ko*lenI), 1:p]%*%SigmaInv%*%G[1:p, (p + 1):(p + n_simul_ko*lenI)]
      
      s <- svd(Sig_conditional)
      V <- s$v
      if (length(Sig_conditional) == 1){
        D <- s$d
        D[D < 0] <- 0
        D <- sqrt(D)
        X_ko = as.numeric(Re(V * D)) * t(matrix(rnorm(n*(n_simul_ko*lenI)), (n_simul_ko*lenI)))
      } else{
        D <-diag(c(s$d))
        D[D < 0] <- 0
        D <- sqrt(D)
        X_ko = t(Re(V %*% D) %*% matrix(rnorm(n*(n_simul_ko*lenI)), (n_simul_ko*lenI)))
      }
      
      X_ko = t(Re(V %*% D) %*% matrix(rnorm(n*(n_simul_ko*lenI)), (n_simul_ko*lenI)))
      
      X_ko_entropy[[ib, i_run]] <- X_ko
      for (i in 1:n){
        mu_conditional <- matrix(rep(mu[I], n_simul_ko), ncol = 1) + G[(p + 1):(p + n_simul_ko*lenI), 1:p]%*%SigmaInv%*%(X[i,] - matrix(mu, ncol = 1))
        X_ko_entropy[[ib, i_run]][i, ] <- X_ko_entropy[[ib, i_run]][i, ] + mu_conditional
      }
    }
  }
  
  for (i_run in 1:n_runs){
    y <- ys[[i_run]]
    X <- Xs[[i_run]]
    Zs = zeros(p, 1)
    tilZs <- zeros(p, n_simul_ko)
    for (ib in 1:max_n_ko_batches){
      I <- which(clustered_batches == ib)
      lenI <- length(I)
      Xall <- scale(cbind(Xs[[i_run]], X_ko_entropy[[ib, i_run]]))
      
      nlambda = 5*(1 + n_simul_ko)*p
      lambda_max = max(abs(t(Xall) %*% y))/n
      lambda_min = lambda_max/2000
      k = (0:(nlambda - 1))/nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
      
      inds <- sample(1:ncol(Xall))
      Z <- glmnet(Xall[, inds], y, family = "gaussian", intercept = F, lambda = lambda, standardize = F, thresh = 10^(-10), maxit = 10^7)
      first_nonzero <- function(x) match(T, abs(x) > 0)
      indices <- apply(Z$beta, 1, first_nonzero)
      names(indices) <- NULL
      Z = ifelse(is.na(indices), 0, Z$lambda[indices] * n)
      Z[inds] <- Z
      
      Zs[I] <- Z[I]
      for (ik in 1:n_simul_ko){
        tilZs[I, ik] <- t(Z[p+(((ik-1)*lenI+1) : ((ik-1)*lenI+lenI))])
      }
    }
    ranks_entropy[[i_run]] <- scores_to_ranks(cbind(Zs, tilZs))
  }
  
  for (alpha in alpha_vec){
    FDP <- rep(0, n_runs)
    power <- rep(0, n_runs)
    prop <- matrix(0, nrow = n_runs, ncol = p)
    for (i_run in 1:n_runs){
      beta = betas_all[[i_run]]
      ranks <- ranks_entropy[[i_run]]
      mdt <- new("mdt", type ='mirandom', BC1=1, poolDecoys = F, randomTieBreak = T, lambda_lt_c = F, threshold = 1/2, lambda = 1/2)
      results <- number_of_discoveries(ranks[,1], ranks[,2:ncol(ranks)], mdt, alpha, F, seeds[i_run])
      if (length(results[[2]]) == 0){
        FDP[i_run] <- 0
        power[i_run] <- 0
      } else {
        FDP[i_run] <- sum(results[[2]] & (beta == 0))/max(results[[1]], 1)
        discovery_list <- results[[2]] & (beta != 0)
        power[i_run] <- sum(discovery_list)/sum(beta != 0)
      }
      prop[i_run, ] <- results[[5]]/pmax(results[[6]], 1)
    }
    est_FDR_entropy <- c(est_FDR_entropy, mean(FDP))
    est_power_entropy <- c(est_power_entropy, mean(power))
    est_FDR_var_entropy <- c(est_FDR_var_entropy, var(FDP)/n)
    est_power_var_entropy <- c(est_power_var_entropy, var(power)/n)
    if (alpha == 1/9){
      prop_null_target_wins_entropy <- c(prop_null_target_wins_entropy, prop(asdf, 2, mean))
    }
  }
}

dataset = data.frame(batch_overall, alpha_overall, est_FDR_SDP, est_power_SDP, est_FDR_var_SDP, est_power_var_SDP, est_FDR_entropy, est_power_entropy, est_FDR_var_entropy, est_power_var_entropy)

save(dataset,file="Example.Rda")
