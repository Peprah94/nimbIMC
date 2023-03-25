#' Create an auxiliary particle filter algorithm to estimate log-likelihood.
#'
#' @description Create an auxiliary particle filter algorithm for a given NIMBLE state space model.
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model.
#' @param nodes  A character vector specifying the latent model nodes
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author  Nicholas Michaud
#' @details
#' Each of the \code{control()} list options are described in detail here:
#' \describe{
#' \item{lookahead}{The lookahead function used to calculate auxiliary weights.  Can choose between \code{'mean'} and \code{'simulate'}.
#'  Defaults to \code{'simulate'}.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter.  Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function), \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}. Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (\code{TRUE}), or only for the most recent time point (\code{FALSE})}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}. \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time. This need only be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#'  }
#'
#' The auxiliary particle filter modifies the bootstrap filter (\code{\link{buildBootstrapFilter}})
#' by adding a lookahead step to the algorithm: before propagating particles from one time
#' point to the next via the transition equation, the auxiliary filter calculates a weight
#' for each pre-propogated particle by predicting how well the particle will agree with the
#' next data point.  These pre-weights are used to conduct an initial resampling step before
#' propagation.
#'
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWsamp} modelValues object.
#'
#'   The auxiliary particle filter uses a lookahead function to select promising particles before propagation.  This function can eithre be the expected
#'   value of the latent state at the next time point (\code{lookahead = 'mean'}) or a simulation from the distribution of the latent state at the next time point (\code{lookahead = 'simulate'}), conditioned on the current particle.
#'
#'  @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of an auxiliary particle filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.

#' @export
#'
#' @family particle filtering methods
#' @references Pitt, M.K., and Shephard, N. (1999). Filtering via simulation: Auxiliary particle filters. \emph{Journal of the American Statistical Association} 94(446): 590-599.
#'
#' @examples
#' ## For illustration only.
#' exampleCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(.8 * x0, var = 1)
#'   y[1] ~ dnorm(x[1], var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(.8 * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t], var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#' my_AuxF <- buildAuxiliaryFilter(model, 'x',
#'                 control = list(saveAll = TRUE, lookahead = 'mean'))
#' ## Now compile and run, e.g.,
#' ## Cmodel <- compileNimble(model)
#' ## Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' ## logLik <- Cmy_AuxF$run(m = 1000)
#' ## ESS <- Cmy_AuxF$returnESS()
#' ## aux_X <- as.matrix(Cmy_AuxF$mvEWSamples, 'x')

compareModelsPars <- function(models = list(),
                              n.chains = NULL,
                              modelNames = NULL,
                              method = "bm",
                              metrics = list(timeRun = TRUE,
                                             ESS = TRUE,
                                             efficiency = TRUE,
                                             MCse = TRUE)){
  #retrieving parameters needed for the function
  timeRun <- metrics[["timeRun"]]
  ESS <- metrics[["ESS"]]
  efficiency <- metrics[["efficiency"]]
  MCse <- metrics[["MCse"]]

  #function for g
  # interesting in estimating second moments
  #g needs to be defined
  g <- function(x){
  return(sum(x^2))
  }

  #assign names for models
  modelsLength <- length(models)
  if(is.null(modelNames)) modelNames = paste("Model", 1:modelsLength)

  if(is.null(n.chains)) n.chains = 1

  ######################
  # Times run
  ##########################

  if(timeRun) timesRet <- lapply(models, function(x){
    nDim <- length(x$timeRun)
    if(nDim == 1){
      x$timeRun
    }else{
      x$timeRun$all.chains
    }
  })%>%
    do.call('c', .)

  if(ESS) ESSret <- lapply(models, function(x){
    nDim <- length(x$samples[[1]])
    if(nDim > 2){
      ret <- lapply(as.list(1:n.chains), function(y){mcmcse::multiESS(as.matrix(x$samples),
                                                                      method = method,
                                                                      r = 1,
                       size = NULL, g = NULL, adjust = TRUE,
                       blather = TRUE)})%>%
        do.call('c', .)%>%
        mean(.)
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){mcmcse::multiESS(as.matrix(x$samples[[y]][[2]]),
                                                                      method = method,
                                                                      r = 1,
                                                               size = NULL,
                                                               g = NULL,
                                                               adjust = TRUE,
                                                               blather = TRUE)})%>%
        do.call('c', .)%>%
        mean(.)
    }
  })%>%
    do.call('c', .)

  #setting names for efficiency
  colnames(ESSret) <- modelNames

  #####################
  # Monte Carlo Sample Error
  ######################
  if(MCse) MCseret <- lapply(models, function(x){
    nDim <- length(x$samples[[1]])
    if(nDim > 2){
      ret <- lapply(as.list(1:n.chains), function(y){
        N <- nrow(as.matrix(x$samples))
        mcse <- mcmcse::mcse.multi(as.matrix(x$samples),
                                   method = "bm",
                                   g = g,
                         blather = TRUE)
        mcseRet <- c(mcse$cov/N)
        return(mcseRet)
        }%>%
          do.call("c", .)%>%
          mean(.))
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){
        N <- nrow(as.matrix(x$samples[[y]]))
        mcse <- mcmcse::mcse.multi(as.matrix(x$samples[[y]]),
                                   method = "bm",
                                   g = g,
                                   blather = TRUE)
        mcseRet <- c(mcse$cov/N)
        return(mcseRet)
      }%>%
        do.call("c", .)%>%
        mean(.))
    }
  })%>%
    do.call('rbind', .)

  #setting names for model
  rownames(MCseret) <- modelNames

#########################
# estimating efficiency
#######################

  if(efficiency) efficiencyRet <- ESSret/timesRet
  #setting names for efficiency
  colnames(efficiencyRet) <- modelNames

  retDataFrame <- data.frame(timesRun = timesRet,
                             ess = ESSret,
                             efficiency = efficiencyRet,
                             mcse = MCseret)
  return(retDataFrame)
}


#' Create an auxiliary particle filter algorithm to estimate log-likelihood.
#'
#' @description Returns Monte Carlo Standard error, effective sample size and efficiency of individual parameters
#'
#' @param models A list of MCMC models we are interested in returning results for
#' @param nodes A vector name of nodes we are interested in returning results for
#' @param method  The method to use for estimating the Monte Carlo standard error. For details, check the vignette of mcmcse package.
#' @author  Kwaku Peprah Adjei
#' @export
#'
#' @family updating particle filters
#' @references Pitt, M.K., and Shephard, N. (1999). Filtering via simulation: Auxiliary particle filters. \emph{Journal of the American Statistical Association} 94(446): 590-599.
#'
#' @examples
#' ## For illustration only.
#' exampleCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(.8 * x0, var = 1)
#'   y[1] ~ dnorm(x[1], var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(.8 * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t], var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#' my_AuxF <- buildAuxiliaryFilter(model, 'x',
#'                 control = list(saveAll = TRUE, lookahead = 'mean'))
#' ## Now compile and run, e.g.,
#' ## Cmodel <- compileNimble(model)
#' ## Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' ## logLik <- Cmy_AuxF$run(m = 1000)
#' ## ESS <- Cmy_AuxF$returnESS()
#' ## aux_X <- as.matrix(Cmy_AuxF$mvEWSamples, 'x')

compareModelsIndividualPars <- function(models = list(),
                                        modelNames = NULL,
                                        n.chains = NULL,
                                        nodes = c(),
                                        method = "bm", #parameterisations for mcse.mat
                              metrics = list(timeRun = TRUE,
                                             ESS = TRUE,
                                             efficiency = TRUE,
                                             MCse = TRUE)){
  timeRun <- metrics[["timeRun"]]
  ESS <- metrics[["ESS"]]
  efficiency <- metrics[["efficiency"]]
  MCse = metrics[["MCse"]]

  #assign names for models
  modelsLength <- length(models)
if(is.null(modelNames)) modelNames = paste("Model", 1:modelsLength)

  #assign default number of chains
  if(is.null(n.chains)) n.chains = 1

  #assign chain names
  chainNames <- paste0("chain", 1:n.chains)

  ##################
  # Times run
  ############
  if(timeRun) timesRet <- lapply(models, function(x){
    nDim <- length(x$timeRun)
    if(nDim == 1){
      x$timeRun
    }else{
      x$timeRun$all.chains
    }
  })%>%
    do.call('c', .)
  ####################
  # estimate Monte Carlo Standard error
  ##################
  if(MCse) MCseRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    nDim <- length(x$samples[[1]])
    if(nDim >2 ){
      ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][,nodes]),
                                                               method = "bm",
                                                               g = NULL)%>%
                                                              as.data.frame()%>%
                                                              dplyr::select(se)
                    colnames(seEst) <- modelNames[i]
                    return(seEst)
                    })

      names(ret) <- chainNames

      ret$all.chains <- do.call("cbind", ret)%>%
        rowMeans(.)
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][[2]][,nodes]),
                                                                               method = "bm",
                                                                               g = NULL)%>%
        as.data.frame()%>%
        dplyr::select(se)
      colnames(seEst) <- modelNames[i]
      return(seEst)
      })

      names(ret) <- chainNames

      ret$all.chains <- do.call("cbind", ret)%>%
        rowMeans(.)
    }

    return(ret)
  })

  if(ESS) ESSret <- lapply(models, function(x){
    nDim <- length(x$samples[[1]])
    if(nDim > 2){
        ret <- coda::effectiveSize(x$samples[,nodes])
    }else{
      ret <- coda::effectiveSize(x$samples[[2]][,nodes])
    }
  })%>%
    do.call('c', .)

  ggmcmc::ggs_effective(ggs(x$samples))

  if(efficiency) efficiencyRet <- ESSret/timesRet

  retDataFrame <- data.frame(timesRun = timesRet,
                             ess = ESSret,
                             efficiency = efficiencyRet)
  return(retDataFrame)
}
