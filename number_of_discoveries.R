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
  
  target_dist <- cumsum(as.matrix(obs_wins[sort_idxs0] & ((beta[tie_perm])[sort_idxs0] != 0)))
  decoy_dist <- cumsum(as.matrix(obs_wins[sort_idxs0] & ((beta[tie_perm])[sort_idxs0] == 0)))
  asdf <- cumsum(as.matrix((beta[tie_perm])[sort_idxs0] == 0))
  return(list(total_no, discovery_inds, fdps, target_dist, decoy_dist, asdf, obs_wins[sort_idxs0]))
}