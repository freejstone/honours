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