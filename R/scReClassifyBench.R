# list.of.packages <- c("randomForest", "e1071", "gmodels", "mclust")
# new.package.list <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if (length(new.package.list)) install.packages(new.package.list)
#
# # potential packages
# # "amap", "aricode", "edgeR", "biomaRt"
#
# library(randomForest)
# library(e1071)
# library(gmodels)
#
#
# scReClassifyBench <- function(mat, label, seed, classifier="svm", L=10, prob=FALSE, percentVar=0.8) {
#   dat.selected <- matPCs(dat.processed, percentVar, seed = seed)
#
#   cellTypes <- label
#
#   colnames(dat.selected$pcs) <- paste("V", 1:ncol(dat.selected$pcs), sep = "")
#   names(cellTypes) <- paste("V", 1:ncol(dat.selected$pcs), sep = "")
#
#   pred <- multiAdaSampling(dat.selected$pcs, cellTypes, seed=seed, classifier=classifier, percent=1, L=L)$final
#   # rf.pred <- multiAdaSampling(dat.selected$pcs, cellTypes, seed=1, classifier="rf", percent=1, L=10)$final
#
#   adjustedRandIndex(pred, cellTypes)
#   # adjustedRandIndex(rf.pred, cellTypes)
#
#
#
# }
