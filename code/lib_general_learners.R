#!/usr/bin/Rscript --vanilla
library(parallel)
library(doMC)
library(doRNG)
library(foreach)
library(glmnet)
library(ranger)

FitRandomForest <- function(x, y, ntrees=1000){
  set.seed(45322)
  ranger(x=x, y=y, importance="impurity", num.trees = ntrees, 
         verbose=F)
}

FitLasso <- function(x, y, parallel_ = TRUE){
  set.seed(65468)
  lasso.model <- cv.glmnet(x, y,
                           standardize = FALSE,
                           alpha = 1,
                           lambda = 10^seq(10, -3, length=100),
                           parallel = parallel_)

  return(lasso.model)
}

FitRidgeReg <- function(x, y, parallel_ = TRUE){
  set.seed(65462)
  ridge.model <- cv.glmnet(x, y,
                           standardize = FALSE,
                           alpha = 0,
                           lambda = 10^seq(10, -3, length=100),
                           parallel = parallel_)

  return(ridge.model)
}