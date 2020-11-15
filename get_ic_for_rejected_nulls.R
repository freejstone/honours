get_ic_for_rejected_hyp = function(cs, c_for_candidate, pv_thr_for_HB_c, per_c_obs_wins, per_c_cntd_decoy_wins, p_obs){
  ncs <- numel(cs)
  cs0 <- t(cbind(0, t(cs[2:length(cs)])))
  if (c_for_candidate >= 1){
    ic0 <- which(cs0 <= 1/2)[length(which(cs0 <= 1/2))]
    c_for_candidate <- min(1/c_for_candidate, max(0, (sum(per_c_obs_wins[,ic0]) - sum(per_c_cntd_decoy_wins[,ic0]))/(sum(per_c_obs_wins[,ic0]) + sum(per_c_cntd_decoy_wins[,ic0])))) 
  } else if (c_for_candidate > 0){
    c_for_candidate <- c_for_candidate
  } else if (c_for_candidate == 0 || c_for_candidate == -1 || c_for_candidate == -2){
    c_for_candidate <- cs[length(cs)]
    c_denom <- round2(1/cs[1], n=1)
    p_obs_cnts <- accumarray(rbind(matrix(round2(p_obs * c_denom, n=1)), matrix(1:c_denom)), val = matrix(1, ncol = 1, nrow = length(p_obs) + c_denom))
    p_obs_cnts <- p_obs_cnts - 1
    for (i in 1 : min(ncs-1, c_denom-2)){
      midpoint = ceil((i + c_denom) / 2)
      if (mod(i + c_denom, 2) == 0 || c_for_candidate > -2){
        pv <- pbinom(sum(p_obs_cnts[(i+1):midpoint])-1, sum(p_obs_cnts[(i+1):c_denom]), (midpoint-i) / (c_denom-i), lower.tail = F)
        cond_break <- pv > pv_thr_for_HB_c
      }
      if (c_for_candidate == -1 && mod(i + c_denom, 2) != 0){
        pv2 <- pbinom(sum(p_obs_cnts[i+1:midpoint-1])-1, sum(p_obs_cnts[i+1:c_denom]), (midpoint-1-i) / (c_denom-i), lower.tail = F)
        cond_break <- cond_break || pv2 > pv_thr_for_HB_c
      } else if (c_for_candidate == -2 && mod(i + c_denom, 2) != 0){
        pv <- pbinom(sum(p_obs_cnts[i+1:midpoint-1])-1, sum(p_obs_cnts[i+1:midpoint-1]) + sum(p_obs_cnts[midpoint+1:c_denom]), 1/2, lower.tail = F)
        cond_break <- pv > pv_thr_for_HB_c
      }
      if (cond_break){
        c_for_candidate = cs[i]
        break
      }
    }
  } else {
    print('oooooh no')
  }
  ic <- which(cs0 <= c_for_candidate)[length(which(cs0 <= c_for_candidate))]
  return(ic)
}