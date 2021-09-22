#' bAccuracy
#'
#' This function calculates the accuracy of the prediction to the true label.
#'
#' @param cls.truth A character vector of true class label.
#' @param final A vector of final classified label prediction from 
#' \code{multiAdaSampling}.
#' @return An accuracy value.
#' @author Pengyi Yang, Taiyun Kim
#'
#' @examples
#' data("gse87795_subset_sce")
#'
#' mat.expr <- gse87795_subset_sce
#' cellTypes <- gse87795_subset_sce$cellTypes
#'
#' # Get dimension reduced matrix. We are using `logNorm` assay from `mat.expr`.
#' mat.pc <- matPCs(mat.expr, assay = "logNorm")
#'
#' # Here we are using Support Vector Machine as a base classifier.
#' result <- multiAdaSampling(mat.pc, cellTypes, classifier = "svm",
#' percent = 1, L = 10)
#'
#' final <- result$final
#'
#' # Balanced accuracy
#' bacc <- bAccuracy(cellTypes, final)
#'
#' @export bAccuracy
bAccuracy <- function(cls.truth, final) {
    
    if (length(cls.truth) != length(final)) {
        stop("`cls.truth` and `final` objects have different length.")
    }
    
    gs <- names(table(cls.truth))
    acc <- c()
    acc <- unlist(lapply(seq_len(length(gs)), function(i) {
        sum((cls.truth == final)[cls.truth==gs[i]]) / sum(cls.truth==gs[i])
    }))
    mean(acc)
}
