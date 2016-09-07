ICENormalization <- function(htc_obj, max_iter = 200) {
  #' ICE normalization of Hi-C matrix
  #'
  #' @param htc_obj HTC object.
  #' @param max_iter Integer, maximum number of iteration
  #' @return HTC object of the ICE normalized matrix
  return(normICE(htc_obj, max_iter = max_iter))
}


ShortestPathCorrection <- function(htc_obj, resolution, alpha = 1) {
  #' shortest-path correction of contact matrix
  #'
  #' @param htc_obj HTC object.
  #' @param Hi-C contact matrix resolution
  #' @param alpha scale parameter for the transformation from contact matrix
  #' to distance matrix 
  #' @return matrix, the shortest-path corrected matrix.
  hic_matrix <- as.matrix(intdata(htc_obj))
  hic_dist <- 1 / (hic_matrix ^ alpha)
  # hic_graph_weight <- melt(hic_dist, id.vars = 'region')
  hic_graph_weight <- melt(hic_dist)
  hic_edges <- rbind(as.numeric(hic_graph_weight[[1]]),
                     as.numeric(as.character(hic_graph_weight[[2]])))
  if (resolution) {
    hic_graph <- add_edges(make_empty_graph(nrow(hic_dist)),
                           hic_edges / resolution + 1,
                           weight = hic_graph_weight[[3]])
  } else {
    hic_graph <- add_edges(make_empty_graph(nrow(hic_dist)), hic_edges,
                           weight = hic_graph_weight[[3]])
  }
  hic_shortest_path <- distances(hic_graph, algorithm = 'dijkstra')
  colnames(hic_shortest_path) <- colnames(hic_matrix)
  rownames(hic_shortest_path) <- colnames(hic_matrix)
  return(hic_shortest_path)
}


RemoveRegions <- function(hic_matrix, centromere_region, unmapped_egion,
                          resolution, n_mask = 10) {
  #' Remove regions (regions without reads mapped, centromere,
  #' regions near centromere)
  #'
  #' @param hic_matrix Matrix, contact matrix.
  #' @param centromere_region Numeric vector with length 2, the first element
  #' should be the start of centromere and the second should be the end.
  #' @param unmapped_egion Numeric vector, regions without mapped reads.
  #' @param resolution Integer, resolution of Hi-C.
  #' @param n_mask Integer, the number of bins to be removed around centromere.
  #' @return Matrix, Hi-C matrix without centromere and unmapped bins
  near_centromere_region <- c(centromere_region[1] - n_mask * resolution,
                              centromere_region[2] + n_mask * resolution)
  rm_region <- seq(floor(near_centromere_region[1] / resolution) * resolution,
                   floor(near_centromere_region[2] / resolution) * resolution,
                   by = resolution)
  rm_region <- union(rm_region, unmapped_egion)
  remain_region_row <- setdiff(rownames(hic_matrix), rm_region)
  remain_region_col <- setdiff(colnames(hic_matrix), rm_region)
  return(hic_matrix[remain_region_row, remain_region_col])
}
