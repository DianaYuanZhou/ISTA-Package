library(MASS)
library(Metrics)

test_that("test if main function of ISTA is working properly", {
  expect_output(ISTA.main(as.matrix(bodyfat[,-1]), as.matrix(bodyfat[,1]))
)
}
)
