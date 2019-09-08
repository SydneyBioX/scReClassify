#' multi Adaptive Sampling function
#'
#' Performs multiple Adaptive sampling to train a classifier model
#'
#' @param data A dimension reduced matrix.
#' @param label A vector of label information for each sample.
#' @param seed Seed before base classifier model.
#' @param classifier Base classifier model, either  "SVM" or "RF"
#' @param percent Percentage of samples to select at each iteration during.
#' @param L Number of ensembles. Default to 10.
#' @param prob logical flag to return sample's probabilities to each class.
#' @param balance logical flag to down sample large cell types classes to the median of all class sizes.
#' @param iter A number of iterations to perform adaSampling.
#' @return A final prediction, probabilities for each cell type and the model are returned as a list.
#' @usage multiAdaSampling(dat, lab)
#' @author Pengyi Yang
#' @importFrom randomForest randomForest
#' @importFrom e1071 svm
#' @examples
#' data("GSE87795_liver.development.data")
#'
#' mat.expr = GSE87795_liver.development.data$data
#' cellTypes = GSE87795_liver.development.data$cellTypes
#'
#' mat.pc = matPCs(mat.expr)
#'
#' result = multiAdaSampling(mat.pc, cellTypes, seed = 1, classifier = "svm", percent = 1, L = 10)
#' @export multiAdaSampling
multiAdaSampling <- function(data, label, seed=1, classifier="svm", percent=1, L=10, prob=FALSE, balance=FALSE, iter=3) {
  models <- list()
  for(l in 1:L) {
    set.seed(seed+l)
    X <- data
    Y <- label

    model <- c()
    prob.mat <- c()
    for (i in 1:iter) {

      if (classifier == "rf") {
        model <- randomForest::randomForest(t(X), factor(Y), ntree = 100)
        prob.mat <- predict(model, newdata=t(data), type="prob")
      }
      if (classifier == "svm") {
        tmp <- t(X)
        rownames(tmp) <- NULL
        model <- e1071::svm(tmp, factor(Y), probability = TRUE)
        prob.mat <- attr(predict(model, t(data), decision.values = FALSE, probability = TRUE), "probabilities")
      }

      X <- c()
      Y <- c()
      for(j in 1:ncol(prob.mat)) {
        voteClass <- prob.mat[label==colnames(prob.mat)[j],]
        idx <- c()
        if (balance == FALSE) {
          idx <- sample(1:nrow(voteClass), size=nrow(voteClass)*percent, replace = TRUE, prob=voteClass[,j])
        } else {
          sampleSize <- round(median(table(label)))
          if (nrow(voteClass) > sampleSize) {
            idx <- sample(1:nrow(voteClass), size=sampleSize*percent, replace = TRUE, prob=voteClass[,j])
          } else {
            idx <- sample(1:nrow(voteClass), size=nrow(voteClass)*percent, replace = TRUE, prob=voteClass[,j])
          }
        }

        X <- cbind(X, data[, rownames(voteClass)[idx]])
        Y <- c(Y, label[rownames(voteClass)[idx]])
      }
    }
    models[[l]] <- model
  }

  predictMat <- matrix(0, nrow=ncol(data), ncol=length(table(label)))
  final <- c()
  for (l in 1:L) {
    if (classifier == "svm") {
      tmp <- attr(predict(models[[l]], newdata=t(data), probability = TRUE), "prob")[,names(table(label))]
      predictMat <- predictMat + tmp
    } else {
      tmp <- predict(models[[l]], newdata=t(data), type="prob")[,names(table(label))]
      predictMat <- predictMat + tmp
    }
  }

  if(prob==TRUE) {
    final <- apply(predictMat, 1, max)
    names(final) <- names(table(label))[apply(predictMat, 1, which.max)]
  } else {
    final <- names(table(label))[apply(predictMat, 1, which.max)]
  }

  # return (final)
  return(
    list(
      final = final,
      models = models,
      prob = predictMat
    )
  )
}
