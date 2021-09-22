test_that(
    "Case multiAdaSampling.1: Missing parameters",
    {
        dat = matrix(rnorm(1200), nrow= 12)
        label = 1:ncol(dat)
        expect_error(multiAdaSampling())
        expect_error(multiAdaSampling(data = dat))
        expect_error(multiAdaSampling(label = label))
    }
)

test_that(
    "Case multiAdaSampling.2: Expected output",
    {
        set.seed(123)
        # very basic data with strong signal with 2 classes
        dat = matrix(rnorm(1200), ncol = 12)
        colnames(dat) = paste0("Cell", seq_len(12))
        rownames(dat) = paste0("Gene", seq_len(100))
        label = rep(c("CT_A", "CT_B"), each = 6)
        names(label) = colnames(dat)
        noise_label = label
        
        pc = matPCs(dat, percentVar = 0.8)
        result = multiAdaSampling(
            pc, 
            label = noise_label, 
            classifier = "svm"
        )
        
        object_names = c("final", "models", "prob")
        prob_dim = c(12L, 2L) # 100 cells, 2 cell types
        
        expect_identical(object_names, names(result))
        expect_identical(prob_dim, dim(result$prob))
        expect_s3_class(result$models[[1]], "svm")
        expect_identical(12L, length(result$final))
    }
)



