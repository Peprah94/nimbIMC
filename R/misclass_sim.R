#' Set the parameters
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param nspecies,nsites, nvisit, nvalidate, nreported,
#' beta0, beta1, gamma0, gamma1, sigma_cov,
#' sigma_class, sigma_corr
#' @return A list of parameters for simulating data
#' @export
#' @examples
#'  misclassSetPars(nspecies = 2, nsites = 1000,nvisit = 1,nvalidate = 100,
#'  nreported = 3, beta0 = c(1,0), beta1 = c(0,1), gamma0 = matrix(c(2, 1, 0,
#'                 1, 3, 0), nrow = 2,
#'                ncol = 3, byrow = TRUE),
#' gamma1 = matrix(c(-1, 1,0 ,
#'                  -1, 2, 0), nrow = 2,
#'                ncol = 3, byrow = TRUE),
#' sigma_cov = 1,
#'sigma_class = 1,
#' sigma_corr = matrix(c(1,0.8,0.8,1), nrow = 2, byrow= T))
misclassSetPars <- function(nspecies,nsites, nvisit, nvalidate, nreported,
                    beta0, beta1, gamma0, gamma1, sigma_cov,
                    sigma_class, sigma_corr, ncov){
  if(length(nsites) != 1) stop("Dimension of number of sites is wrong")
  if(length(nvisit) != 1) stop("Dimension of number of visits is wrong")
  if(length(nspecies) != 1) stop("Dimension of number of species is wrong")
  if(length(sigma_cov) != 1) stop("Dimension of number of variance of ecological process covariate is wrong")
  if(length(sigma_class) != 1) stop("Dimension of number of variance of ecological process covariate is wrong")
  if(dim(sigma_corr)[1] != dim(sigma_corr)[1] ) stop("Dimension of number of correlation matrix wrong of ecological process covariate is wrong")
  if(length(beta0) != nspecies) stop("Dimension of number of intercept is wrong")
  if(dim(beta1)[1] != nspecies) stop("Dimension of number of covariate effects is wrong")
  #if(length(cov) != nsites) stop("Dimension of number of covariates is wrong")
input <- list(
  constants =  list(
    n.species = nspecies, # The number of species
    n.sites = nsites, # The number of sites
    n.visit= nvisit, # The number of visits
    n.validate = nvalidate,
    n.reported = nreported,
    n.cov = ncov
    ),
fixed_effects = list(
  beta0 = beta0,#,0), # intercept for each species
  beta1 = beta1,#,-2), #Covariate effect for each species,
  gamma0 = gamma0,
  gamma1 = gamma1
),
covariates = list(
  #cov = matrix(rnorm(nsites*ncov, 0, sigma_cov), nrow = ncov, ncol = nsites), #NB: n.site = 400,
  #cov_omega = rnorm(nsites, 0, sigma_class),
  cov = matrix(rnorm(nsites*ncov, 0, sigma_cov), nrow = ncov, ncol = nsites), #NB: n.site = 400,
  cov_omega = rnorm(nsites, 0, sigma_class),
  cov_corr = MASS::mvrnorm(nsites, mu = rep(0, nspecies), Sigma = sigma_corr)
),
true_contrasts = list(
  gamma0 = t(apply(gamma0,1,
                  function(x){x[1:(length(x)-1)] - x[length(x)]})),
  gamma1 = t(apply(gamma1,1,
                   function(x){x[1:(length(x)-1)] - x[length(x)]}))
)
)
return(input)
}



#' Generating data
#'
#' This function takes the returned output of the misclassSetPars function
#' and generates multispecies misclassified data
#'
#' @param input, type, seed
#' @return Simulated data for the misclassification project
#' @export
misclassGenData <- function(input, type, seed){ # made all of the inputs null so can see everything feeds through
 # set seed
set.seed(seed)

   #Retrieving Parameter values
  n.sites = input$constants$n.sites # number of sites
  n.visit=input$constants$n.visit #number of visits
  p.tag = input$constants$p.tag # probability of correct classification
  n.species=input$constants$n.species # Number of species
  beta0 <- input$fixed_effects$beta0 # intercept
  beta1 <- input$fixed_effects$beta1 # Covariate effect
  gamma0 <- input$fixed_effects$gamma0
  gamma1 <- input$fixed_effects$gamma1
  cov <- input$covariates$cov #covariate
  cov_omega <- input$covariates$cov_omega
  cov_corr <- input$covariates$cov_corr
  n.reported <- input$constant$n.reported
  n.validate <- input$constant$n.validate
  n.cov <- input$constant$n.cov
  true1 <- input$true_constants$gamma0
  true2 <- input$true_constants$gamm1

  #Check if parameter values are right
  if(length(n.sites) != 1) stop("Dimension of number of sites is wrong")
  if(length(n.visit) != 1) stop("Dimension of number of visits is wrong")
  if(length(n.species) != 1) stop("Dimension of number of species is wrong")
  if(length(beta0) != n.species) stop("Dimension of number of intercept is wrong")
  if(dim(beta1)[1] != n.cov) stop("Dimension of number of covariate effects is wrong")
  if(dim(cov)[1] != n.cov | dim(cov)[2] != n.sites) stop("Dimension of number of covariates is wrong")

  #Estimating mean intensity
  lambda <- matrix(NA, nrow=n.species, ncol=n.sites)
  if(type == "variable"| type == "constant"){
  #for(spe.tag in 1:n.species){
   # for(site.tag in 1:n.sites){
      lambda <- exp(beta0 + beta1%*%cov)
   # }
  #}
  }else{
    # for(spe.tag in 1:n.species){
    #   for(site.tag in 1:n.sites){
    #     lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag]*cov_corr[site.tag,1])
    #   }
    # }
    cov <- rbind(cov[1,], cov_corr[,1])
    lambda <- exp(beta0 + beta1%*%cov)
  }

#Proportion of each species as defined by Equation (2) in the paper
 prop <-  matrix(NA, ncol=n.sites ,nrow=n.species)
for(site.tag in 1:n.sites){
  prop[,site.tag]<- proportions(lambda[1:n.species, site.tag])
}

#Check if the proportions are correct and sum up to 1
for(site.tag in n.sites){
  for(spe.tag in 1:n.species){
    if(any(prop[spe.tag, site.tag] <0) | any(prop[spe.tag, site.tag] >1) ) stop("Entry in Probability wrong")
    if(round(sum(prop[1:n.species,site.tag ]),0) != 1) stop("Sum of probabilities is not equal to 1")
  }
}

z <- array(NA, dim=c(n.species, n.visit, n.sites)) #Detect/Non-Detection
C <- array(NA, dim=c(n.visit,n.sites)) # Obs detect/non-detect
Y <- numb <- numb_with_NA <- array(NA, dim=c(n.visit, n.sites)) # Verified species with Machine Learning

#Verified species
for(site.tag in 1:n.sites){
  for(visit.tag in 1:n.visit){
    C[visit.tag, site.tag] <- extraDistr::rcat(1, prop[1:n.species,site.tag])
  }
}


alpha <- omega <- array(NA, dim = c(n.species, n.reported, n.sites))
if(type == "variable"){

for(site.tag in 1: n.sites){
  for(spe.tag in 1:n.species){
    for(report.tag in 1:(n.reported)){
      alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + gamma1[spe.tag, report.tag] * cov_omega[site.tag])
    }
  }
}
# Confusion Matrix for the Species
for(site.tag in 1: n.sites){
  for(spe.tag in 1:n.species){
    for(report.tag in 1:n.reported){
      omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
    }
  }
}
}


if(type == "constant"){
  for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      for(report.tag in 1:(n.reported)){
        alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag])
      }
    }
  }
  # Confusion Matrix for the Species
  for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      for(report.tag in 1:n.reported){
        omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
      }
    }
  }
}

if(type == "correlation"){

  for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      for(report.tag in 1:(n.reported)){
        alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + gamma1[spe.tag, report.tag] * cov_corr[site.tag,2])
      }
    }
  }
  # Confusion Matrix for the Species
  for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      for(report.tag in 1:n.reported){
        omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
      }
    }
  }
}


for(visit.tag in 1:n.visit){
  for(site.tag in 1:n.sites){
    Y[visit.tag,site.tag] <- extraDistr::rcat(1,omega[C[visit.tag,site.tag],1:n.reported, site.tag])
  }
}
#}
message("printing the tabulation of the simulated data")
print(table(C[1,],Y[1,]))

#Subsetting training and validation set
index_for_splitting <- sort(sample(1:n.sites, n.validate, replace = FALSE))
Validation_C <- C[, index_for_splitting]
Validation_Y <- Y[, index_for_splitting]
validation_indices_for_mismatch = index_for_splitting[which(Validation_C != Validation_Y)]
validation_indices_for_match = index_for_splitting[which(Validation_C == Validation_Y)]
validation_mismatchC = C[, validation_indices_for_mismatch]
validation_mismatchY = Y[, validation_indices_for_mismatch]
validation_matchC = C[, validation_indices_for_match]
validation_matchY = Y[, validation_indices_for_match]
#table(Validation_C,Validation_Y)
C[1:n.visit, index_for_splitting] <- NA
#Y[1:n.visit, index_for_splitting] <- NA

if(type == "variable" | type == "constant"){
  covariate_group = 1
}else{
  covariate_group = 2
}


  dataout <- list(C=C, #Verified data
                  Y=Y, # Reported data
                  Lam=lambda, #Mean intensity
                  omega=omega, # Confusion matrix
                  cov=cov, # covariate
                  cov_omega = cov_omega,
                  cov_corr = cov_corr,
                  prop=prop, #proportion of intensity,
                  true_values = c(true1, true2),
                  Validation_C = Validation_C,
                  Validation_Y = Validation_Y,
                  validation_indices_for_mismatch = validation_indices_for_mismatch,
                  index_for_splitting = index_for_splitting,
                  validation_mismatchC = validation_mismatchC,
                  validation_mismatchY = validation_mismatchY,
                  validation_indices_for_match = validation_indices_for_match,
                  validation_matchC = validation_matchC,
                  validation_matchY = validation_matchY,
                  covariate_group = covariate_group
                  )
  return(dataout)

}

