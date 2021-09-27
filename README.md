
# scReClassify <img src='img/scReClassify_sticker.png' align="right" height="138.5" />

`scReClassify` is a post hoc cell type classification of single-cell
RNA-sequencing data. Using semi-supervised learning algorithm,
adaSampling, to correct cell type annotation from noise.

# Getting started

## Vignette

For more detailed instuction, you can find the vignette at our website:
<https://sydneybiox.github.io/scdney/>.

### Installation

#### Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scReClassify")
```

#### GitHub

``` r
devtools::install_github("SydneyBioX/scReClassify", build_opts = c("--no-resave-data", "--no-manual"))
library(scReClassify)
```

For `devtool version
< 2.0.0`,

``` r
devtools::install_github("SydneyBioX/scReClassify", build_vignettes = TRUE)
library(scReClassify)
```


``` r
suppressPackageStartupMessages({
    library(scReClassify)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
})
```


# Usage

![alt text](https://github.com/SydneyBioX/scReClassify/raw/master/img/scReClassify.jpg)


Current version of this package is implemented to run with `svm` and `randomForest` classifiers.

## Load data

``` r
data(gse87795_subset_sce)
dat <- gse87795_subset_sce
cellTypes <- gse87795_subset_sce$cellTypes

# number of clusters
nCs <- length(table(cellTypes))

# This demo dataset is already pre-processed
dim(dat)
```

## Part A. scReClassify (Demonstration)

### Dimension reduction

``` r
reducedDim(dat, "matPCs") = matPCs(dat, assay = "logNorm", 0.7)
```

### Synthetic noise (Demonstration purpose)

Here in this example, we will synthetically generate varying degree of noise in sample labels.

``` r
lab <- cellTypes

set.seed(1)
noisyCls <- function(dat, rho, cls.truth){
    cls.noisy <- cls.truth
    names(cls.noisy) <- colnames(dat)
    for(i in seq_len(length(table(cls.noisy)))) {
        # class label starts from 0
        if (i != length(table(cls.noisy))) {
            cls.noisy[sample(which(cls.truth == names(table(cls.noisy))[i]), floor(sum(cls.truth == names(table(cls.noisy))[i]) * rho))] <- names(table(cls.noisy))[i+1]
        } else {
            cls.noisy[sample(which(cls.truth == names(table(cls.noisy))[i]), floor(sum(cls.truth == names(table(cls.noisy))[i]) * rho))] <- names(table(cls.noisy))[1]
        }
    }
  
    print(sum(cls.truth != cls.noisy))
    return(cls.noisy)
}

cls.noisy01 <- noisyCls(t(reducedDim(dat, "matPCs")), rho=0.1, lab)
cls.noisy02 <- noisyCls(t(reducedDim(dat, "matPCs")), rho=0.2, lab)
cls.noisy03 <- noisyCls(t(reducedDim(dat, "matPCs")), rho=0.3, lab)
cls.noisy04 <- noisyCls(t(reducedDim(dat, "matPCs")), rho=0.4, lab)
cls.noisy05 <- noisyCls(t(reducedDim(dat, "matPCs")), rho=0.5, lab)
```

### Use scReClassify to correct mislabeled cell types.

Here in this example, we will only use `Support Vector machine (svm)` as base classifier.

#### Benchmark evaluation

``` r
###################################
# SVM
###################################
base <- "svm"
set.seed(1)
result = lapply(seq_len(10), function(j) {
    final <- multiAdaSampling(dat, cls.noisy01, reducedDimName = "matPCs", 
                                classifier=base, percent=1, L=10)$final
    ari01 <- mclust::adjustedRandIndex(lab, final)
    acc01 <- bAccuracy(lab, final)
    
    final <- multiAdaSampling(dat, cls.noisy02, reducedDimName = "matPCs", 
                                classifier=base, percent=1, L=10)$final
    ari02 <- mclust::adjustedRandIndex(lab, final)
    acc02 <- bAccuracy(lab, final)
    
    final <- multiAdaSampling(dat, cls.noisy03, reducedDimName = "matPCs", 
                                classifier=base, percent=1, L=10)$final
    ari03 <- mclust::adjustedRandIndex(lab, final)
    acc03 <- bAccuracy(lab, final)
    
    final <- multiAdaSampling(dat, cls.noisy04, reducedDimName = "matPCs", 
                                classifier=base, percent=1, L=10)$final
    ari04 <- mclust::adjustedRandIndex(lab, final)
    acc04 <- bAccuracy(lab, final)
    
    final <- multiAdaSampling(dat, cls.noisy05, reducedDimName = "matPCs", 
                                classifier=base, percent=1, L=10)$final
    ari05 <- mclust::adjustedRandIndex(lab, final)
    acc05 <- bAccuracy(lab, final)
    
    c(
      acc01 = acc01,
      acc02 = acc02,
      acc03 = acc03,
      acc04 = acc04,
      acc05 = acc05,
      ari01 = ari01,
      ari02 = ari02,
      ari03 = ari03,
      ari04 = ari04,
      ari05 = ari05
    )
})

result = do.call(rbind, result)
acc = result[,seq_len(5)]
colnames(acc) = seq(from=0.1,to=0.5,by=0.1)

ari = result[,seq(from= 6, to = 10)]
colnames(ari) = seq(from=0.1,to=0.5,by=0.1)


plot.new()
par(mfrow = c(1,2))
boxplot(acc, col="lightblue", main="SVM Accuracy", 
        ylim=c(0.45, 1), xlab = "rho", ylab = "Accuracy")
points(x=seq_len(5), y=c(
    bAccuracy(lab, cls.noisy01), 
    bAccuracy(lab, cls.noisy02),
    bAccuracy(lab, cls.noisy03), 
    bAccuracy(lab, cls.noisy04),
    bAccuracy(lab, cls.noisy05)), 
    col="red3", pch=c(2,3,4,5,6), cex=1)
boxplot(ari, col="lightblue", main="SVM ARI", 
        ylim=c(0.25, 1), xlab = "rho", ylab = "ARI")
points(x=seq_len(5), y=c(
    mclust::adjustedRandIndex(lab, cls.noisy01), 
    mclust::adjustedRandIndex(lab, cls.noisy02),
    mclust::adjustedRandIndex(lab, cls.noisy03), 
    mclust::adjustedRandIndex(lab, cls.noisy04),
    mclust::adjustedRandIndex(lab, cls.noisy05)), 
    col="red3", pch=c(2,3,4,5,6), cex=1)
```


## Part B. scReClassify (mislabeled cell type correction)

``` r
# PCA procedure
reducedDim(dat, "matPCs") = matPCs(dat, assay = "logNorm", 0.7)


# run scReClassify
set.seed(1)
cellTypes.reclassify <- multiAdaSampling(dat, cellTypes, 
                                        reducedDimName = "matPCs",
                                        classifier = "svm", percent = 1, L = 10)

# Verification by marker genes
End <- c("ITGA2B", "ITGB3")

# check examples
idx <- which(cellTypes.reclassify$final != cellTypes)

cbind(original=cellTypes[idx], reclassify=cellTypes.reclassify$final[idx]) %>%
    DT::datatable()

mat <- assay(dat, "logNorm")

c1 <- mat[, which(cellTypes=="Endothelial Cell")]
c2 <- mat[, which(cellTypes=="Erythrocyte")]
c3 <- mat[, which(cellTypes=="Hepatoblast")]
c4 <- mat[, which(cellTypes=="Macrophage")]
c5 <- mat[, which(cellTypes=="Megakaryocyte")]
c6 <- mat[, which(cellTypes=="Mesenchymal Cell")]
cs <- rainbow(length(table(cellTypes)))


# (example 1 E13.5_C14)
#####
par(mfrow=c(1,2))
marker <- End[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], 
        c4[marker,], c5[marker,], c6[marker,], 
        col=cs, main=marker, 
        names=c("Others", "Others", "Others", "Orignal", 
                "Reclassified", "Others"), las=2, xlab = "Labels",
        ylab = "log2FPKM")
points(5, mat[marker, which(colnames(mat) %in% "E13.5_C14")], 
        pch=16, col="red", cex=2)

marker <- End[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], 
        c4[marker,], c5[marker,], c6[marker,], 
        col=cs, main=marker, 
        names=c("Others", "Others", "Others", "Orignal", 
                "Reclassified", "Others"), las=2, xlab = "Labels", 
        ylab = "log2FPKM")
points(5, mat[marker, which(colnames(mat) %in% "E13.5_C14")], 
        pch=16, col="red", cex=2)
```


# Reference
scReClassify is published in BMC Genomics. Please refer to the following article for more details on method implementation and evaluation.

Taiyun Kim, Kitty Lo, Thomas A. Geddes, Hani Jieun Kim, Jean Yee Hwa Yang & Pengyi Yang (2019) scReClassify: post hoc cell type classification of single-cell RNA-seq data. BMC Genomics, 20:913. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6305-x
