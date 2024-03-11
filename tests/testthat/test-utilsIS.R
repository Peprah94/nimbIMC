devtools::load_all()



test_that("Gamma estimation function", {
  code <- nimble::nimbleCode({

    for(i in 1:2){
      beta[i] ~ dnorm(0, tau = 0.001)
    }

    a ~ dnorm(0, tau = 0.001)
    #Fitting the inla with the simulated parameters
    # inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)


    for(i in 1:N){
      logit(linpred[i]) <-  a + beta[1]*x[i, 1] + beta[2]*x[i, 2]
    }

    #Bivariate linear model specification
    for(i in 1:N){
      y[i] ~ dbinom(prob = linpred[i], size = 1)
    }

  })

  N <- 10; m = 5
  idm_data <- list(y= c(rbinom(N-5, 1, prob = 0.5), rep(NA, N-5)),
                   x = matrix(rnorm(N*2, 0, 1), ncol = 2))

  constants = list(N = N)

  inits <-  function(){list(beta =c(1,1),
                            a = 1)
  }

  initsList <- inits()


  mwtc <- nimble::nimbleModel(code,
                              data = idm_data,
                              constants =constants,
                              inits = initsList)


cmwtc <- compileNimble(mwtc)

ecoParams <- mwtc$expandNodeNames("y")


mvEWSamples <- modelValues(modelValuesConf(vars = "y",
                                           types = "double",
                                           sizes = 10))
mvEWSamples[["y"]][[1]] <- rbinom(N, 1, prob = 0.5)

rt <- adaptingFunction(model = mwtc,
mvEWSamples,
m = 5,
iNode = 1,
dfTdist = 1,
ecoParams = ecoParams,
obsParams = ecoParams,
ecoParamsProposal = "binomial",
obsParamsProposal = "normal",
adaptive = FALSE)

meanBeta = matrix(NA, ncol = N, nrow = 2)
meanBeta[1, ] <- rep(0,N)
sigmaBeta = array(NA, dim = c(N,N,2))
sigmaBeta[,,1] <- diag(N)
meanDisc = rep(0, 2)
obsParamsEstUpd <- matrix(0, nrow = m, ncol = N + 1)
ecoParamsEstUpd <- matrix(0, nrow = m, ncol = N + 1)

rt$run(meanBeta = meanBeta, sigmaBeta = sigmaBeta, meanDisc = meanDisc, ecoParamsEst = ecoParamsEstUpd, obsParamsEst = obsParamsEstUpd)
rt$updateMeanEcoParams()
rt$updateMeanObsParams()
rt$updateSigmaObsParams()

#IF we don't adapt, I expect the mean and sd to be the same as before
expect_equal(rt$updateMeanEcoParams(), meanDisc[1])
expect_equal(rt$updateMeanObsParams(), meanBeta[1, ])
expect_equal(rt$updateSigmaObsParams(), sigmaBeta[,,1])

#if we adapt

rt <- adaptingFunction(model = mwtc,
                       mvEWSamples,
                       m = 5,
                       iNode = 1,
                       dfTdist = 1,
                       ecoParams = ecoParams,
                       obsParams = ecoParams,
                       ecoParamsProposal = "binomial",
                       obsParamsProposal = "normal",
                       adaptive = TRUE)

meanBeta = matrix(NA, ncol = N, nrow = 2)
meanBeta[1, ] <- rep(0,N)
sigmaBeta = array(NA, dim = c(N,N,2))
sigmaBeta[,,1] <- diag(N)
meanDisc = rep(0, 2)

rt$run(meanBeta = meanBeta, sigmaBeta = sigmaBeta, meanDisc = meanDisc, ecoParamsEst = ecoParamsEstUpd, obsParamsEst = obsParamsEstUpd)
rt$updateMeanEcoParams()
rt$updateMeanObsParams()
rt$updateSigmaObsParams()

expect_equal(length(rt$updateMeanEcoParams()), length(meanDisc[1]))
expect_equal(length(rt$updateMeanObsParams()), length(meanBeta[1, ]))
expect_equal(dim(rt$updateSigmaObsParams()), dim(sigmaBeta[,,1]))

})
