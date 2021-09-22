#' matPCs function
#'
#' @description Performs PCA on a given matrix and returns a dimension reduced 
#' matrix which captures at least 80% (default) of overall variability.
#'
#' @param data An expression matrix or a SingleCellExperiment object.
#' @param assay An assay to select if \code{data} is a SingleCellExperiment 
#' object
#' @param percentVar The percentage of variance threshold. This is used to 
#' select number of Principal Components.
#' @details This function performs PCA to reduce the dimension of the gene
#' expression matrix limited from 10 to 20 PCs.
#' @return Dimensionally reduced matrix.
#' @author Pengyi Yang, Taiyun Kim
#' @export matPCs
#' @examples
#' data("gse87795_subset_sce")
#'
#' mat.expr <- gse87795_subset_sce
#'
#' mat.pc <- matPCs(mat.expr, assay = "logNorm")
#'
#' # to capture at least 70% of overall variability in the dataset,
#' mat.dim.reduct.70 <- matPCs(mat.expr, assay = "logNorm", 0.7)
#'
#' @importFrom stats prcomp sd
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @import SingleCellExperiment
matPCs <- function(data, assay = NULL, percentVar=0.8) {
    
    if (percentVar < 0 | percentVar > 1) {
        stop("percentVar must be a number between 0 and 1")
    }
    # If SCE
    if (is(data, "SingleCellExperiment")) {
        if (is.null(assay)) {
            stop("assay parameter cannot be NULL")
        }
        mat <- SummarizedExperiment::assay(data, assay)
    } else {
        mat <- data
    }
    
    if (is(data, "matrix") | is(data, "data.frame")) {
        pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
    } else if ((is(data, "SingleCellExperiment")) & 
        ("PCA" %in% reducedDimNames(data))) {
        pca <- SingleCellExperiment::reducedDim(data, "PCA")
    } else if (is(data, "SingleCellExperiment")) {
        pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
        if (length(SingleCellExperiment::reducedDimNames(data))) {
            rd <- SingleCellExperiment::reducedDims(data)
            rd$PCA = pca$x
        } else {
            SingleCellExperiment::reducedDims(data) <- list(PCA = pca$x)
        }
    }

    # genes are rows, cells are cols
    pcs <- c()
    sdev <- apply(pca$x, 2, sd)
    eigs <- sdev^2
    top <- which(cumsum(eigs)/sum(eigs) > percentVar)[1]
    m <- ifelse(top > 20, 20, top)
    m <- ifelse(m < 10, 10, m)
    pcs <- pca$x[,seq_len(m)]

    return(t(pcs))
}
