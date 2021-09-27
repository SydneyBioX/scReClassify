test_that(
    "Case matPCs.1: Missing parameters",
    {
        expect_error(matPCs())
        expect_error(matPCs(percentVar = 0.8))
    }
)

test_that(
    "Case matPCs.2: Expected outputs",
    {
        set.seed(123)
        # 50 cells 400 genes
        mat = matrix(rnorm(20000), nrow = 50)
        pc = prcomp(mat)
        eigs = pc$sdev^2
        p = cumsum(eigs)/sum(eigs)
        
        result1 = min(which(p > 0.5))
        result2 = min(which(p > 0.8))
        # Since result2 > 20
        result2 = 20L
        
        # Number of PCs returned
        expect_identical(result1, ncol(matPCs(t(mat), percentVar = 0.5)))
        expect_identical(result2, ncol(matPCs(t(mat), percentVar = 0.8)))
    }
)
