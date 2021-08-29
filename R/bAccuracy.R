#' bAccuracy
#'
#' This function calculates the accuracy of the prediction to the true label.
#'
#' @param cls.truth A vector true class label.
#' @param final A vector of final classified label prediction from XXX.
#' @return An accuracy value.
#' @usage bAccuracy(cls.truth, final)
#' @author Pengyi Yang, Taiyun Kim
#'
#' @examples
#' data("gse87795_subset")
#'
#' mat.expr = dat
#'
#' mat.pc = matPCs(mat.expr)
#'
#' result = multiAdaSampling(mat.pc, cellTypes, classifier = "svm",
#' percent = 1, L = 10)
#'
#' final = result$final
#'
#' # Balanced accuracy
#' bacc = bAccuracy(cellTypes, final)
#'
#' @export bAccuracy
bAccuracy <- function(cls.truth, final) {
    gs <- names(table(cls.truth))
    acc <- c()
    for (i in seq_len(length(gs))) {
        acc <- c(acc, sum((cls.truth == final)[cls.truth==gs[i]]) /
                sum(cls.truth==gs[i]))
    }
    mean(acc)
}
