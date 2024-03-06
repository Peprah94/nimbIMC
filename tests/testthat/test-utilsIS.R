devtools::load_all()

test_that("Gamma estimation function", {
  res <- gammaEstimation(iNode=1, ecoParamsProposal = "binomial", obsParamsProposal = "normal", dfTdist = 1, timeIndex = 10)
  sigmaBeta <- array(NA, dim = c(2,2,3))
  sigmaBeta[,,1] <- diag(2)
  expect_equal(class(res$run(meanBeta = diag(2), sigmaBeta = sigmaBeta, meanDisc = c(0.5, 0.5), ecoParamsVals = c(1,1,1,0,1,0,0), obsParamsVals = c(-1, 0.5))), "numeric")
})
