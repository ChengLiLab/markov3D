# Rcpp::sourceCpp("computeSSD.cpp")
ComputeSSD <- function(transition.mat, iter.epsilon = 1e-8, iter.max = 1e5) {
  #' compute steady-state distribution from transition matrix
  #'
  #' @param transition.mat Matrix, transition matrix.
  #' @param iter.epsilon Float, the tolerance of iterative computation of SSD
  #' @param iter.max Int, the max number of iteration
  #' @return Named vector, steady-state distribution.
  ssd <- Re(rARPACK::eigs(t(transition.mat), 1)$vector)
  if (length(ssd) == 0 | any(ssd < 0)) {
    ssd <- as.vector(computeSSD(transition.mat, iter.epsilon, iter.max))
  }  # compute ssd through iteration
  ssd <- ssd / sum(ssd)
  names(ssd) <- colnames(transition.mat)
  return(ssd)
}
