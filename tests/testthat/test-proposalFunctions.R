devtools::load_all()

test_that("Testing Binomial proposal distribution",{
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

N <- 10
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

ecoParams <- mwtc$expandNodeNames("y")
binNodeToSimulate <- ecoParams[!mwtc$isData(ecoParams)]
binNodeToFix <- ecoParams[mwtc$isData(ecoParams)]
binNodesToFixVals <- rep(1, 5)

binProposaltest <- binProposal(mwtc, binNodeToSimulate, binNodeToFix,binNodesToFixVals, ecoParams, size = 1)

cmwtc <- compileNimble(binProposaltest, mwtc)
expect_equal(class(cmwtc$binProposaltest$run(meanDisc = 0.5, n = 5, m = 10))[1], "matrix")
})


test_that("Testing Poisson proposal distribution",{
  code <- nimble::nimbleCode({

    for(i in 1:2){
      beta[i] ~ dnorm(0, tau = 0.001)
    }

    a ~ dnorm(0, tau = 0.001)
    #Fitting the inla with the simulated parameters
    # inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)


    for(i in 1:N){
      log(linpred[i]) <-  a + beta[1]*x[i, 1] + beta[2]*x[i, 2]
    }

    #Bivariate linear model specification
    for(i in 1:N){
      y[i] ~ dpois(lambda = linpred[i])
    }

  })

  N <- 10
  idm_data <- list(y= rpois(n = N, lambda = 2),
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

  ecoParams <- mwtc$expandNodeNames("y")

  poisProposaltest <- poisProposal(mwtc, ecoParams, lowerBound = rep(0, 10), includeLowerBound = 0)

  cmwtc <- compileNimble(poisProposaltest, mwtc)
  expect_equal(class(cmwtc$poisProposaltest$run(n = 10, meanDisc = 2, m = 10))[1], "matrix")
  expect_error(cmwtc$poisProposaltest$run(n = 5, meanDisc = 2, lowerBound = rep(0, 10), includeLowerBound = 0))
})

test_that("Testing prior proposal distribution",{
  code <- nimble::nimbleCode({

    for(i in 1:2){
      beta[i] ~ dnorm(0, tau = 0.001)
    }

    a ~ dnorm(0, tau = 0.001)
    #Fitting the inla with the simulated parameters
    # inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)


    for(i in 1:N){
      log(linpred[i]) <-  a + beta[1]*x[i, 1] + beta[2]*x[i, 2]
    }

    #Bivariate linear model specification
    for(i in 1:N){
      y[i] ~ dpois(lambda = linpred[i])
    }

  })

  N <- 10
  idm_data <- list(y= rpois(n = N, lambda = 2),
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

  ecoParams <- mwtc$expandNodeNames("y")

  priorProposaltest <- priorProposal(mwtc, ecoParams)

  cmwtc <- compileNimble(priorProposaltest, mwtc)
  expect_equal(class(cmwtc$priorProposaltest$run(n = 10, meanDisc = 2, m = 10))[1], "matrix")
})

test_that("Testing normal proposal distribution",{
  code <- nimble::nimbleCode({

    for(i in 1:2){
      beta[i] ~ dnorm(0, tau = 0.001)
    }

    a ~ dnorm(0, tau = 0.001)
    #Fitting the inla with the simulated parameters
    # inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)


    for(i in 1:N){
      linpred[i] <-  a + beta[1]*x[i, 1] + beta[2]*x[i, 2]
    }

    #Bivariate linear model specification
    for(i in 1:N){
      y[i] ~ dnorm(linpred[i], sd = 1)
    }

  })

  N <- 10
  idm_data <- list(y= rnorm(n = N),
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

  ecoParams <- mwtc$expandNodeNames("y")

  normalProposaltest <- normalProposal(mwtc, ecoParams)

  cmwtc <- compileNimble(normalProposaltest, mwtc)

  expect_equal(class(cmwtc$normalProposaltest$run(mean = rep(0, 10), sigma=diag(10), m = 10))[1], "matrix")
})

test_that("Testing student-T proposal distribution",{
  code <- nimble::nimbleCode({

    for(i in 1:2){
      beta[i] ~ dnorm(0, tau = 0.001)
    }

    a ~ dnorm(0, tau = 0.001)
    #Fitting the inla with the simulated parameters
    # inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)


    for(i in 1:N){
      linpred[i] <-  a + beta[1]*x[i, 1] + beta[2]*x[i, 2]
    }

    #Bivariate linear model specification
    for(i in 1:N){
      y[i] ~ dnorm(linpred[i], sd = 1)
    }

  })

  N <- 10
  idm_data <- list(y= rnorm(n = N),
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

  ecoParams <- mwtc$expandNodeNames("y")

  studentProposaltest <- studentProposal(mwtc, ecoParams, df = 2)

  cmwtc <- compileNimble(studentProposaltest, mwtc)

  expect_equal(class(cmwtc$studentProposaltest$run(mean = rep(0, 10), sigma=diag(10), m = 10))[1], "matrix")
})
