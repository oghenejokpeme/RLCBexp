#!/usr/bin/Rscript --vanilla

library(parallel)
library(doMC)
library(doRNG)
library(foreach)
library(glmnet)
library(randomForest)
library(kernlab)
library(caret)
library(xgboost)
library(quadprog)

perfCR <- function(x, y){
  Rinv <- solve(chol(t(x) %*% x))
  C <- cbind(rep(1,ncol(Rinv)), diag(ncol(Rinv)))
  b <- c(1, rep(0, ncol(Rinv)))
  d <- t(y) %*% x
  qsolve <- solve.QP(Dmat = Rinv, 
                     factorized = TRUE, 
                     dvec = d, 
                     Amat = C, 
                     bvec = b, 
                     meq = 1)
  weights <- qsolve$solution
    
  return(weights)
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

FitRandomForest <- function(x, y, ntrees=500, parallel_=TRUE){
  if (parallel_){
    ncores = getDoParWorkers()
    set.seed(98364)
    rf.model <- foreach(ntree = rep(floor(ntrees/ncores), ncores),
                        .combine = combine, 
                        .multicombine = TRUE, 
                        .packages = "randomForest") %dorng% {
                randomForest(x, y, importance = TRUE, ntree = ntree)  
                }

    return(rf.model)
  } else {
    set.seed(98364)
    rf.model <- randomForest(x, y, importance = TRUE, ntree = ntrees) 
        
    return(rf.model)
  }
}

FitSvmReg <- function(x, y, parallel_=TRUE){
  ctrl <- trainControl(method = "cv", number = 3, allowParallel = parallel_)
  svmgrid <- expand.grid(sigma = 2^seq(-3, 3, by = 1),
                         C = c(0.01, 0.1, 1))
  set.seed(53475)
  svm.selection <- train(x = x, 
                         y = y,
                         scaled = FALSE,
                         method = "svmRadial",
                         metric = "RMSE", 
                         epsilon = 0.001,
                         trControl = ctrl,
                         tuneGrid = svmgrid) 
  svm.model <- svm.selection$finalModel
    
  return(svm.model)
}
