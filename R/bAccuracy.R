#' bAccuracy
#'
#' This function calculates the accuracy of the prediction to the true label.
#'
#' @param cls.truth A vector true class label.
#' @param final A vector of final classified label prediction from XXX.
#' @return An accuracy value.
#' @usage bAccuracy(truth, final)
#' @author Pengyi Yang
#' @export bAccuracy
bAccuracy <- function(cls.truth, final) {
  gs <- names(table(cls.truth))
  acc <- c()
  for (i in 1:length(gs)) {
    acc <- c(acc, sum((cls.truth == final)[cls.truth==gs[i]]) / sum(cls.truth==gs[i]))
  }
  mean(acc)
}
