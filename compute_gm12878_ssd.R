# Compute SSD from GM12878
library(data.table)
library(dplyr)
library(ggplot2)

# Set libpath
current_libpath <- .libPaths()
new_libpath <- '~/lustre/markov_3d/source/'
.libPaths(unique(c(current_libpath, new_libpath)))

# Load package markov3d
library(devtools)
load_all('~/lustre/markov_3d/source/markov3d/')
data("chrom.size")
data("centromere.pos")

# Load other packages
library(data.table)
library(dplyr)
library(RPostgreSQL)

# Read data: GM12878 primary chr1 50kb resolution
resolution <- 50e3
n_mask <- 5
chrom <- 'chr1'
gm12878_chr1_path <- '~/lustre/hic_2014_cell/GM12878/GM12878_primary/50kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_50kb.RAWobserved'
chrom_size <- filter(chrom.size, chr == chrom)$size
gm12878_chr1 <- LoadRawMatrix(gm12878_chr1_path, chrom_size, resolution)

# ICE normalization
gm12878_chr1 <- MatrixToHTCobj(gm12878_chr1, chrom, chrom_size, resolution)
gm12878_chr1_ice <- ICENormalization(gm12878_chr1)
gm12878_chr1_ice_mat <- intdata(gm12878_chr1_ice)
gm12878_chr1_ice_mat <- as.matrix(gm12878_chr1_ice_mat)
plot(colSums(gm12878_chr1_ice_mat))

# Shortest-path correction
gm12878_chr1_sp <- ShortestPathCorrection(gm12878_chr1_ice, resolution)

# Remove unmapped and centromere regions
ummaped_region <- colnames(gm12878_chr1_ice_mat[
  colSums(gm12878_chr1_ice_mat) == 0
  ])
centromere_region <- c(filter(centromere.pos, chr == chrom)$start, 
                       filter(centromere.pos, chr == chrom)$end)
gm12878_chr1_sp_rmna <- RemoveRegions(gm12878_chr1_sp, centromere_region,
                                      ummaped_region ,resolution)


# Load H3K4me1 chr1
pg <- dbDriver("PostgreSQL")
con <- dbConnect(pg, dbname = "markov3d_db",
                 host = "localhost", port = 5432, user = "wangyn")
table_name <- 'gm12878_h3k4me1_50kb_chr1'
if (dbExistsTable(con, table_name)) {
  query_cmd <- paste0("SELECT * FROM ", table_name)
  h3k4me1_chr1 <- as.data.table(dbGetQuery(con, query_cmd))
} else {
  sprintf('The table %s does not exist', table_name)
}
dbDisconnect(con)


# Estimate transition matrix
gm12878_chr1_tmat <- EstimateTransitionMat(
  gm12878_chr1_ice_mat, gm12878_chr1_sp_rmna)


# Compute SSD
gm12878_chr1_ssd <- ComputeSSD(gm12878_chr1_tmat)
gm12878_chr1_ssd_dt <- data.table(
  bin_start = as.numeric(names(gm12878_chr1_ssd)),
  ssd = gm12878_chr1_ssd
)


# Relationship between SSD and ChIP-seq signals
gm12878_merge <- merge(gm12878_chr1_ssd_dt, h3k4me1_chr1, by = 'bin_start')
plot(gm12878_merge$ssd, gm12878_merge$rpm)
cor(gm12878_merge$ssd, gm12878_merge$rpm, method = 'spearman')
