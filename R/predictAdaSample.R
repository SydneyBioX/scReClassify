#' @import stats
predictAdaSample <- function(models, data, labelNames, classifier = "svm",
                            prob = FALSE) {
    L = length(models)

    predictMat <- matrix(0, nrow=ncol(data), ncol=length(labelNames))
    final <- c()
    for (l in seq_len(L)) {
        if (classifier == "svm") {
            tmp <- attr(predict(models[[l]],newdata=t(data),probability = TRUE),
                        "prob")[,names(labelNames)]
            predictMat <- predictMat + tmp
        } else {
            tmp <- predict(models[[l]], newdata=t(data),
                            type="prob")[,names(labelNames)]
            predictMat <- predictMat + tmp
        }
    }
    if(prob==TRUE) {
        final <- apply(predictMat, 1, max)
        names(final) <- names(labelNames)[apply(predictMat, 1, which.max)]
    } else {
        final <- names(labelNames)[apply(predictMat, 1, which.max)]
    }

    return(final)
}
