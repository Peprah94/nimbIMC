setwd("~/OneDrive - NTNU/GitHub/Thesis/INLA within NIMBLE/Bayesian Lasso")

#get data

#packages
library(ISLR)
library(glmnet)
data("Hitters")

#Check for NAs
Hitters <- na.omit(Hitters)

#Fit lassso regression
#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)
df <- list(y = y, x = x)
n.beta <- ncol(df$x)

# ml estimates
ml = summary(lm(y~-1 + x, data = df))$coefficients[,1:2]

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)

#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit lasso model for several values of lambda
lasso.mod <- glmnet(x[train, ] , y[train], alpha = 1, lambda = grid,intercept = F)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1,intercept=F)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)


#Fitted values
lasso.fitted <- predict(out, s = bestlam, newx = x)

#data for MCMC
save(df, file="hitters_data.RData")
