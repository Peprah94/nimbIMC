##  Contains code to run bootstrap particle filters.
##  We have a build function (buildBootstrapFilter),
##  and step function.

inlaStepVirtual <- nimbleFunctionVirtual(
  run = function(beta = double(1)
                 ) {
    returnType(double(0))
  }
)

# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.
inlaFStep <- nimbleFunction(
  name = 'inlaFStep',
  contains = inlaStepVirtual,
  setup = function(model,
                   mvEWSamples,
                   fixedVals,
                  # inlaModel,
                   x,
                   y,
                  # interInModel,
                   fam
                   ) {
    #res <- inlaModel(x, y, beta, fixedVals, interInModel, family = fam)
    #copy(model, mvEWSamples, nodes = fixedVals, row = 1)
    N <- length(fixedVals)

#     if(N >1){
#       mult <- TRUE
#       #vals <- rep(0, N)
#     }else{
#       mult <- FALSE
#      # vals <- 0
#     }
#     print(mult)
# print(N)

    # mult <- TRUE
    # mult1 <- FALSE
    # if(length(fixedVals) == 1){
    #   mult <- FALSE
    #   mult1 <- TRUE}
    # print(mult)
    #
    # if(N ==1){
    #   mult <- FALSE
    # }else{
    #   mult <- TRUE
    # }
    #
    # if(N == 1){
    #   vals <- 0
    # }else{
    #   vals <- rep(0, N)
    # }
  },
  run = function(beta = double(1)
                 ) {
    returnType(double(0))

    #vals <- numeric(N, init=FALSE)
   # vals <-
    #ll <-
    #copy(mvEWSamples, model, nodes = fixedVals, row = 1)
    res <- try(nimbleINLA(x, y, beta= beta, fixedVals,  family = fam), silent = TRUE)

#interInModel= interInModel,
    # save results
   # if(!mult){
      #vals <- res[1, 2:N]
   # }
#     mult <<- FALSE
#     #if(mult == TRUE)
# print(mult)
#     if(mult == FALSE){
#       vals <- res[1, 2]
#     }
#     else if(mult == TRUE){
#       vals <- numeric(N, init = FALSE)
#       for(i in seq_along(fixedVals)) vals[i] <- res[1, i+1]
#     }


#if(mult){
#for(i in 1:N){
#vals <- res[1, 2 :N]
#     print(mult)
#     if(!mult) vals <<- res[1, 2]
# print(vals)
#   if(mult) for(i in seq_along(fixedVals)) vals[i] <<- res[1, i+1]



#print(vals)
#}
#}

    lll <- res[1,1]
   #if(mult) values(model, fixedVals) <<- res[1, 2:length(fixedVals)]
    # if(mult == TRUE){
    # for(i in seq_along(fixedVals)){
    #   mvEWSamples[fixedVals[i], 1] <<- res[1, 2]
    # }
    # }

   # if(mult == FALSE) mvEWSamples[fixedVals, 1] <<- res[1, 2]
     #values(model, fixedVals) <<- res[1, 2]

        saveResults(fixedVals, res)
      copy( model, mvEWSamples, nodes = fixedVals, row = 1)
      out <- lll
    #values(model, fixedVals) <<- vals
    #model[[fixedVals]] <<- vals
    #}
    #ßvalues(model, fixedVals) #<<- vals

    return(out)
  },
methods = list(
saveResults = function(fixedVals = character(0),
                       res = double(2)){
  #n <- length(fixedVals)
  vals <- res[1, 2]
  values(model, fixedVals) <<- c(vals)
  #if(n > 1)
    #for(i in seq_along(fixedVals)) vals[i] <- res[1, i + 1]
  #values(model, fixedVals) <<- vals
  #else{
   #vals <<- res[1, 2]
  #}
  return(vals)
  returnType(double(0))
}
)
)


# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.
inlaFStepMultiple <- nimbleFunction(
  name = 'inlaFStepMultiple',
  contains = inlaStepVirtual,
  setup = function(model,
                   mvEWSamples,
                   fixedVals,
                   # inlaModel,
                   x,
                   y,
                   # interInModel,
                   fam
  ) {
    #res <- inlaModel(x, y, beta, fixedVals, interInModel, family = fam)
    #copy(model, mvEWSamples, nodes = fixedVals, row = 1)
    N <- length(fixedVals)

    #     if(N >1){
    #       mult <- TRUE
    #       #vals <- rep(0, N)
    #     }else{
    #       mult <- FALSE
    #      # vals <- 0
    #     }
    #     print(mult)
    # print(N)

    # mult <- TRUE
    # mult1 <- FALSE
    # if(length(fixedVals) == 1){
    #   mult <- FALSE
    #   mult1 <- TRUE}
    # print(mult)
    #
    # if(N ==1){
    #   mult <- FALSE
    # }else{
    #   mult <- TRUE
    # }
    #
    # if(N == 1){
    #   vals <- 0
    # }else{
    #   vals <- rep(0, N)
    # }
  },
  run = function(beta = double(1) # beta is a vector

  ) {
    returnType(double(0))

    #vals <- numeric(N, init=FALSE)
    # vals <-
    #ll <-
    #copy(mvEWSamples, model, nodes = fixedVals, row = 1)
    res <- nimbleINLA(x, y, beta= beta, fixedVals,  family = fam)

    #interInModel= interInModel,
    # save results
    # if(!mult){
    #vals <- res[1, 2:N]
    # }
    #     mult <<- FALSE
    #     #if(mult == TRUE)
    # print(mult)
    #     if(mult == FALSE){
    #       vals <- res[1, 2]
    #     }
    #     else if(mult == TRUE){
    #       vals <- numeric(N, init = FALSE)
    #       for(i in seq_along(fixedVals)) vals[i] <- res[1, i+1]
    #     }


    #if(mult){
    #for(i in 1:N){
    #vals <- res[1, 2 :N]
    #     print(mult)
    #     if(!mult) vals <<- res[1, 2]
    # print(vals)
    #   if(mult) for(i in seq_along(fixedVals)) vals[i] <<- res[1, i+1]



    #print(vals)
    #}
    #}

    lll <- res[1,1]
    #if(mult) values(model, fixedVals) <<- res[1, 2:length(fixedVals)]
    # if(mult == TRUE){
    # for(i in seq_along(fixedVals)){
    #   mvEWSamples[fixedVals[i], 1] <<- res[1, 2]
    # }
    # }

    # if(mult == FALSE) mvEWSamples[fixedVals, 1] <<- res[1, 2]
    #values(model, fixedVals) <<- res[1, 2]

      saveResults(fixedVals, res)
      copy( model, mvEWSamples, nodes = fixedVals, row = 1)
      out <- lll

    #values(model, fixedVals) <<- vals
    #model[[fixedVals]] <<- vals
    #}
    #ßvalues(model, fixedVals) #<<- vals

    return(out)
  },
  methods = list(
    saveResults = function(fixedVals = character(1),
                           res = double(2)){
      n <- length(fixedVals)
      vals <- numeric(n, init = FALSE)
      #r <- character(0)
      #if(n > 1){
      for(i in seq_along(fixedVals)){
        # r <- fixedVals[i]
        vals[i] <- res[1, i + 1]
        #model[[r]] <<- vals[i]
      }
      # }
      #for(i in seq_along(fixedVals)){
      #model[[fixedVals[i]]] <<- vals[i]
      #values(model, fixedVals[i]) <<- vals[i]
      #}
      values(model, fixedVals) <<- c(vals)
      #else{
      # vals <<- res[1, 2]
      #}
      return(vals)
      returnType(double(1))
    }
  )
)





#' Create a bootstrap particle filter algorithm to estimate log-likelihood.
#'
#'@description Create a bootstrap particle filter algorithm for a given NIMBLE state space model.
#'
#' @param model A nimble model object, typically representing a state
#'  space model or a hidden Markov model.
#' @param nodes  A character vector specifying the latent model nodes
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author Daniel Turek and Nicholas Michaud
#' @details
#'
#' Each of the \code{control()} list options are described in detail here:
#' \describe{
#'  \item{thresh}{ A number between 0 and 1 specifying when to resample: the resampling step will occur when the
#'   effective sample size is less than \code{thresh} times the number of particles. Defaults to 0.8. Note that at the last time step, resampling will always occur so that the \code{mvEWsamples} \code{modelValues} contains equally-weighted samples.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter. Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function),  \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}.  Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (TRUE), or only for the most recent time point (FALSE)}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}.  \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#' }
#'
#'  The bootstrap filter starts by generating a sample of estimates from the
#'  prior distribution of the latent states of a state space model.  At each time point, these particles are propagated forward
#'  by the model's transition equation.  Each particle is then given a weight
#'  proportional to the value of the observation equation given that particle.
#'  The weights are used to draw an equally-weighted sample of the particles at this time point.
#'  The algorithm then proceeds
#'  to the next time point.  Neither the transition nor the observation equations are required to
#'  be normal for the bootstrap filter to work.
#'
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWSamples} \code{modelValues} object.
#'
#'  Note that if the \code{thresh} argument is set to a value less than 1, resampling may not take place at every time point.
#'  At time points for which resampling did not take place, \code{mvEWSamples} will not contain equally weighted samples.
#'  To ensure equally weighted samples in the case that \code{thresh < 1}, we recommend resampling from \code{mvWSamples} at each time point
#'  after the filter has been run, rather than using \code{mvEWSamples}.
#'
#' @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of a bootstrap filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.
#'
#' @export
#'
#' @family particle filtering methods
#' @references Gordon, N.J., D.J. Salmond, and A.F.M. Smith. (1993). Novel approach to nonlinear/non-Gaussian Bayesian state estimation. \emph{IEEE Proceedings F (Radar and Signal Processing)}. Vol. 140. No. 2. IET Digital Library, 1993.
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
#' my_BootF <- buildBootstrapFilter(model, 'x',
#'                 control = list(saveAll = TRUE, thresh = 1))
#' ## Now compile and run, e.g.,
#' ## Cmodel <- compileNimble(model)
#' ## Cmy_BootF <- compileNimble(my_BootF, project = model)
#' ## logLik <- Cmy_BootF$run(m = 1000)
#' ## ESS <- Cmy_BootF$returnESS()
#' ## boot_X <- as.matrix(Cmy_BootF$mvEWSamples, 'x')

buildINLAmodel <- nimbleFunction(
  name = 'buildINLAmodel',
  setup = function(model, fam, x, y, control) {
    inlaModel <- extractControlElement(control, 'fit.inla',  NULL)
    fixedVals <- extractControlElement(control, 'fixedVals',  double())

    yExpand <- model$expandNodeNames(y, returnScalarComponents = TRUE)
    y <- model[[y]]
    #y <- c(nimble::values(model, yExpand))

    #save posterior samples
    #modelVals = modelValues(model, m = 1)
    vars <- model$getVarNames(nodes = fixedVals)
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]

      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size <- lapply(size, function(x){
        if(length(x) == 0){
          #length(model$vars)
          return(1)
        }else(
          x
        )
      } )

      #size$beta <- 4

      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))

      fixedVals <- model$expandNodeNames(fixedVals)

      multiple <- TRUE
      if(length(model$expandNodeNames(fixedVals)) == 1) multiple = FALSE
     # vals <- rep(0, length(fixedVals))
      #for(i in seq_along(fixedVals)){
       # values(model, fixedVals[i]) <- c(res[1,fixedVals[i]])
      #}
#vals <- c(res[1,fixedVals])
      inlaStepFunctions <- nimbleFunctionList(inlaStepVirtual)
      #for(iNode in seq_along(nodes)){
      if(multiple == TRUE){
        inlaStepFunctions[[1]] <- inlaFStepMultiple(model,
                                            mvEWSamples,
                                            fixedVals,
                                            # inlaModel,
                                            x,
                                            y,
                                           # interInModel,
                                            fam)
        }else{
          inlaStepFunctions[[1]] <- inlaFStep(model,
                                              mvEWSamples,
                                              fixedVals,
                                              # inlaModel,
                                              x,
                                              y,
                                              # interInModel,
                                              fam)
        }
      #}
      #essVals <- rep(0, length(nodes))

    lastLogLik <- -Inf
  },
  run = function(beta=double(1)#, # beta is a vector
    #interInModel = double()
    ) {
    returnType(double())
    #need this to retuen the saved mvEWSamples
    resize(mvEWSamples, 1)
    # for(i in seq_along(fixedVals)){
    #   mvEWSamples[[fixedVals[i]]][1] <<- vals[i]
    # }
    out <- inlaStepFunctions[[1]]$run(beta)#, interInModel)
    #rr <- inlaStepFunctions[[1]]$mvEWSamples
    logL <- out
    if(logL == -Inf) {lastLogLik <<- logL; return(logL)}
    if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
    if(logL == Inf)  {lastLogLik <<- -Inf; return(-Inf)}
    lastLogLik <<-  logL
    return(logL)
  },
  methods = list(
    getLastLogLik = function() {
      return(lastLogLik)
      returnType(double())
    },
    setLastLogLik = function(lll = double()) {
      lastLogLik <<- lll
    }
  )
)




# ,
# getMVest = function(){
#   return(inlaStepFunctions[[1]]$mvEWSamples)
#   returnType(double())
# }

# buildINLAmodel <- nimbleFunction(
#   name = 'buildINLAmodel',
#   setup = function(x,y, control = list(interInModel, family)) {
#     data <- list(y=y, x=x)
#     data$oset = data$x %*% (beta)
#
#     if(interInModel == 1){
#       formula =  1 + offset(data$oset)
#     }else{
#       formula = - 1 + offset(data$oset)
#     }
#
#     res = INLA::inla(y ~ formula,
#                      data = data,
#                      family = family,
#                      verbose=FALSE,
#                      control.predictor = list(compute = TRUE))
#
#     ret <- res$mlik[1,1]
#
#     lastLogLik <- -Inf
#   },
#   run = function(beta=double(1)) {
#     returnType(double())
#     logL <- ret
#     lastLogLik <<- logL
#     return(logL)
#   },
#   methods = list(
#     getLastLogLik = function() {
#       return(lastLogLik)
#       returnType(double())
#     },
#     setLastLogLik = function(lll = double()) {
#       lastLogLik <<- lll
#     }
#   )
# )
#
