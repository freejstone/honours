optimal_batch_G_analytic <- function(I, Sig, nsko){
  p <- ncol(Sig)
  nI <- setdiff(1:p, I)
  if (length(nI) == 0){
    s_equi = ((nsko + 1)/nsko)*min(eigen(Sig)$values)
  } else {
    s_equi <- ((nsko + 1)/nsko)*min(eigen(Sig[I, I] - Sig[I, nI]%*%as.matrix(solve(Sig[nI, nI]))%*%Sig[nI, I])$values)
  }
  G <- create_candidate_batch_G(I, Sig, nsko, min(1, s_equi))
  return(list(G, s_equi))
}