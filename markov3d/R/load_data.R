LoadRawMatrix <- function(raw_matrix_path, chr_size, resolution) {
  #' Load Hi-C raw matrix stored in sparse format and transform it to full matrix
  #'
  #' @param raw_matrix_path String, file path for sparse matrix, contain Hi-C
  #'  data for one chromosome.
  #' @param chr_size Integer, the length of chromosome.
  #' @param resolution Integer, Hi-C resolution.
  #' @return Matrix, Hi-C raw matrix
  raw_data <- fread(raw_matrix_path)
  if (ncol(raw_data) != 3) {
    stop('The input matrix is not sparse')
  }
  colnames(raw_data) <- c('start_bin', 'end_bin', 'counts')
  matrix_dim <- ceiling(chr_size / resolution)
  dim_name <- seq(0, matrix_dim - 1) * resolution
  dim_complement <- setdiff(dim_name,
                            union(raw_data$start_bin, raw_data$end_bin))
  raw_data.complement <- data.table(
    start_bin = dim_complement,
    end_bin = dim_complement,
    counts = 0
  )
  raw_data <- rbind(raw_data, raw_data.complement)

  raw_data_rev <- data.table(
    start_bin = raw_data$end_bin,
    end_bin = raw_data$start_bin,
    counts = raw_data$counts
  )
  raw_data <- rbind(raw_data, raw_data_rev)
  raw_data <- raw_data[!duplicated(raw_data), ]
  raw_matrix <- dcast(raw_data, start_bin ~ end_bin, fill = 0)
  rownames(raw_matrix) <- raw_matrix[[1]]
  raw_matrix[, 'start_bin'] <- NULL
  return(raw_matrix)
}


MatrixToHTCobj <- function(contact_matrix, chr, chr_size, resolution) {
  #' transform Hi-C matrix to HTC object
  #'
  #' @param contact_matrix Matrix.
  #' @param chr Character, the chromosome of the matrix.
  #' @param chr_size Integer, the length of chromosome.
  #' @param resolution Integer, Hi-C resolution.
  #' @return HTC object of the contact matrix
  if (ncol(contact_matrix) != nrow(contact_matrix)) {
    stop('The matrix is not symmetric')
  }
  gr_obj <- GRanges(
    seqnames = chr,
    ranges = IRanges(seq(0, chr_size, resolution),
                     seq(0, chr_size, resolution) + resolution),
    name = seq(0, chr_size, resolution)
  )
  htc_mat <- Matrix(as.matrix(contact_matrix))
  rownames(htc_mat) <- seq(0, chr_size, resolution)
  colnames(htc_mat) <- seq(0, chr_size, resolution)
  htc_obj <- new("HTCexp", htc_mat, gr_obj, gr_obj)
  return(htc_obj)
}
