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
 if (minEig(s0*(nsko + 1)/nsko) == 0){
    s_equi <- s0*(nsko + 1)/nsko
 } else {
    s_equi <- as.numeric(fzero(minEig, c(s0*(nsko + 1)/nsko, 2))$x)
 }
 G <- create_candidate_batch_G(I, Sig, nsko, min(1, s_equi))
 return(list(G, s_equi))
}

