test_that(
    "Case bAccuracy.1: Missing parameters",
    {
        cls = 1:12
        label = 1:12
        
        expect_error(bAccuracy())
        expect_error(bAccuracy(cls.truth = cls))
        expect_error(bAccuracy(final = label))
    }
)


test_that(
    "Case bAccuracy.2: Different parameter lengths",
    {
        cls = 1:12
        label = 1:10
        
        expect_error(bAccuracy(cls, label))
    }
)

test_that(
    "Case bAccuracy.3: Expected outputs",
    {
        cls = 1:12
        label = 12:1
        exp_result1 = 0
        exp_result2 = 1
        
        expect_identical(exp_result1, bAccuracy(cls, label))
        expect_identical(exp_result2, bAccuracy(cls, cls))
        
    }
)
