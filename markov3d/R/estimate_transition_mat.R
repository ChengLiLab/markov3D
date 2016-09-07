EstimateTransitionMat <- function(ice.matrix, distance.matrix,
                                  method = 'binomial',
                                  t = 1e-3, distance_scale = 1) {
  #' estimate transition matrix from distance matrix
  #'
  #' @param ice.matrix Matrix ICE normalized matrix.
  #' @param distance.matrix Matrix, shortest-path corrected matrix.
  #' @param method, String c('binomial', 'brownian'), binomial or
  #' Wiener Process method to estimate transition matrix,
  #' in the previous manuscript we use binomial method.
  #' @return Matrix, transition matrix.
  if (method != 'binomial' & method != 'brownian') {
    stop('Please choose the right method')
  }
  if (method == 'binomial') {
    diag.p <- diag(ice.matrix) / rowSums(ice.matrix)
    diag.p <- diag.p[colnames(distance.matrix)]
    p.raw <- (1 / distance.matrix) ^ 3
    diag(p.raw) <- 0
    bin.mapped.reads <- rowSums(p.raw)
    transition.mat <- (p.raw / bin.mapped.reads) * (1 - diag.p)
    diag(transition.mat) <- diag.p
  } else if (method == 'brownian') {
    transition.mat <- exp(-((distance.matrix * distance_scale) ^ 2) / (2 * t))
    transition.mat <- transition.mat / rowSums(transition.mat)
  }
  return(transition.mat)
}
