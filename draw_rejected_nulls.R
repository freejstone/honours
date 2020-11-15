draw_rejected_nulls = function(obs_wins, cntd_decoy_wins, p_obs, sort_idxs, decoy_factor, options){
  p <- numel(sort_idxs)
  if (options@n_hyps_per_seg > 0){
    n_hyps_per_seg <- options@n_hyps_per_seg
  } else{
    n_hyps_per_seg <- 1
  }
  is_rejected_null <- matrix(FALSE, p, 1)
  n_obs_wins_so_far <- 0
  n_decoy_wins_so_far <- 0
  n_rejected_so_far <- 0
  i0 <- 1
  last_rjct_ind <- 0
  if (options@bound_n_rejects){
    est_n_rejects <- max(-1, sum(obs_wins) - sum(cntd_decoy_wins))
  }
  
  while (i0 <= p){
    if (options@bound_n_rejects && n_rejected_so_far >= est_n_rejects){
      break
    }
    
    i1 <- min(p, i0 + n_hyps_per_seg - 1)
    indxs <- sort_idxs[i0:i1]
    n_obs_wins_in_segment <- sum(obs_wins[indxs])
    n_obs_wins_so_far <- n_obs_wins_so_far + n_obs_wins_in_segment
    n_decoy_wins_in_segment <- sum(cntd_decoy_wins[indxs])
    n_decoy_wins_so_far <- n_decoy_wins_so_far + n_decoy_wins_in_segment
    est_n_rejects_in_segment <- max( 0, n_obs_wins_so_far - round2((options@BC1 + n_decoy_wins_so_far) * decoy_factor - 1e-10, n = 0) -  n_rejected_so_far)	# -1e-10 for compatibility with previous version: how is x.5 rounded
    if (est_n_rejects_in_segment <= 0){	# do nothing
    n_newly_rejected <- 0
    } else if (numel(indxs) <= est_n_rejects_in_segment && (!isempty(strfind(options@sample, '_in_seg')) || last_rjct_ind + 1 == i0)){
      is_rejected_null[indxs] <- T
      last_rjct_ind <- i1
      n_newly_rejected <- i1 - i0 + 1
    } else if (options@sample == 'uni_in_seg'){
      sampled_reject_inds <- indxs(sample(numel(indxs), est_n_rejects_in_segment))	# uniformly sample indices in the segment
      is_rejected_null[sampled_reject_inds] <- T
      n_newly_rejected <- est_n_rejects_in_segment
    } else if (options@sample == 'uni_in_seg_tgt_wins'){
      samp_indxs <- indxs[obs_wins[indxs]]	# indices of target wins in the segment
      sampled_reject_inds <- samp_indxs[sample(numel(samp_indxs), est_n_rejects_in_segment)]
      is_rejected_null[sampled_reject_inds] <- T
      n_newly_rejected <- est_n_rejects_in_segment
    } else if (options@sample == 'uni_last_rjct'){
      indxs2 <- sort_idxs[last_rjct_ind+1:i1]		# indices since last reject
      newly_rejected_inds <- sample(numel(indxs2), est_n_rejects_in_segment)
      is_rejected_null[indxs2[newly_rejected_inds]] <- T
      last_rjct_ind <- last_rjct_ind + max(newly_rejected_inds)		
      n_newly_rejected <- est_n_rejects_in_segment
    } else if (options@sample == 'uni_last_rjct_tgt_wins'){
      indxs2 <- sort_idxs[last_rjct_ind+1:i1]		# indices since last reject
      tgt_wins_in_indxs2 <- which(obs_wins[indxs2])	# indices of target wins since last reject
      newly_rejected_inds <- tgt_wins_in_indxs2[sample(numel(tgt_wins_in_indxs2), est_n_rejects_in_segment)]
      is_rejected_null[indxs2[newly_rejected_inds]] = T
      last_rjct_ind <- last_rjct_ind + max(newly_rejected_inds)
      n_newly_rejected <- est_n_rejects_in_segment
    } else if (options@sample == 'wtd_p_obs_in_seg'){
      newly_rejected_inds <- zeros(est_n_rejects_in_segment, 1)
      omp <- 1 - p_obs[indxs]
      for (ir in 1 : est_n_rejects_in_segment){
        mps <- omp / sum(omp)
        newly_rejected_inds[ir] <- which(rmultinom(n = 1, size = 1, prob = mps) == 1)
        omp[newly_rejected_inds[ir]] <- 0
      }
      last_rjct_ind <- last_rjct_ind + max(newly_rejected_inds)
      is_rejected_null[indxs[newly_rejected_inds]] <- T
      n_newly_rejected <- est_n_rejects_in_segment
    } else if (options@sample == 'wtd_p_obs_last_rjct_ind'){
      indxs2 <- sort_idxs[(last_rjct_ind+1):i1]		# indices since last reject
      newly_rejected_inds <- zeros(est_n_rejects_in_segment, 1)
      omp <- 1 - p_obs[indxs2]
      for (ir in 1 : est_n_rejects_in_segment){
        mps <- omp / sum(omp)
        newly_rejected_inds[ir] <- which(rmultinom(n = 1, size = 1, prob = mps) == 1)
        omp[newly_rejected_inds[ir]] <- 0
      }
      last_rjct_ind <- last_rjct_ind + max(newly_rejected_inds)
      is_rejected_null[indxs2[newly_rejected_inds]] <-  T
      n_newly_rejected <- est_n_rejects_in_segment
    } else {
      print('how do you want to draw rejected nulls')
    }
    n_rejected_so_far <- n_rejected_so_far + n_newly_rejected
    if (options@n_hyps_per_seg == 0){
      n_hyps_per_seg <- max(n_hyps_per_seg, round2(i1 / n_rejected_so_far, n = 0))
    } else if (options@n_hyps_per_seg == -1){
      n_hyps_per_seg <- max(n_hyps_per_seg, round2((i1-i0+2) / (n_newly_rejected+1), n = 0))
    }
    i0 <- i1+1
  }
    return(is_rejected_null)
}
