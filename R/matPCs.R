#' matPCs function
#'
#' Performs PCA on a given matrix
#'
#' @param data Prediction model generated from "multiAdaSampling".
#' @param percentVar A dimension reduced matrix.
#' @param seed seed for randomisation.
#' @details This function performs PCA to reduce the dimension of the gene expression matrix limited from 10 to 20 PCs.
#' @return Dimensionally reduced matrix.
#' @usage matPCs(dat)
#' @author Pengyi Yang, Taiyun Kim
#' @export matPCs
#' @examples
#' data("GSE87795_liver.development.data")
#'
#' mat.expr = GSE87795_liver.development.data$data
#'
#' mat.pc = matPCs(mat.expr)
#'
#' # to capture at least 70% of overall variability in the dataset,
#' mat.dim.reduct.70 = matPCs(mat.expr, 0.7)
matPCs <- function(data, percentVar=0.8, seed=1) {
  # genes are rows, cells are cols
  pcs <- c()
  set.seed(seed)
  pca <- prcomp(t(data), center = TRUE, scale. = TRUE)
  eigs <- pca$sdev^2
  top <- which(cumsum(eigs)/sum(eigs) > percentVar)[1]
  m <- ifelse(top > 20, 20, top)
  m <- ifelse(m < 10, 10, m)
  pcs <- pca$x[,1:m]

  return(t(pcs))
}
