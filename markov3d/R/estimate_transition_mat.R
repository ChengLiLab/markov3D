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


DistanceCounts <- function(hic_mat, resolution) {
  #' Transform Hi-C full matrix.
  #' @param hic_mat matrix, hic full matrix
  #' @param resolution int, the matrix's resolution
  #' @return a data.table with two columns. The first column is the linear 
  #' genomic distance between two bins and the second is the number of counts
  #' mapped.
  counts_dist <- NULL
  mat_dim <- ncol(hic_mat)
  print('Transforming...')
  for (i in 1:(mat_dim - 1)) {
    if (i %% 100 == 0) {
      sprintf('----------%f%%', i / mat_dim * 100)
    }
    genomic_distance <- resolution * (i - 1)
    sub_mat <- hic_mat[1:(mat_dim + 1 - i), i:mat_dim]
    counts <- diag(sub_mat)
    tmp_dt <- data.table(
      distance = genomic_distance,
      norm_counts = counts
    )
    counts_dist <- rbind(counts_dist, tmp_dt)
  }
  return(counts_dist)
}


PowerLawFitting <- function(dt) {
  #' Fit a power law function between the linear genomic distance and counts
  #' @param dt data.table, the first column is genomic distance and the second
  #' is counts
  #' @return the coefficient of the power law function
  
  # remove 0
  dt <- dt + 1
  
  # fit the model and return the coefficient
  power_model <- lm(log(dt[[1]]) ~ log(dt[[2]]))
  return(power_model$coefficients[2])
}


