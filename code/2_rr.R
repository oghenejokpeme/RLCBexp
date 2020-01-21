#!/usr/bin/Rscript --vanilla
source("lib_general_functions.R")
source("lib_general_learners.R")

PerformBaseExperiments <- function(dfs, all.genes, parallel_ = TRUE){
  base.ilp.train <- dfs$base.ilp.train
  base.ilp.test <- dfs$base.ilp.test
  base.fp.train <- dfs$base.fp.train
  base.fp.test <- dfs$base.fp.test
  base.ilp.fp.train <- dfs$base.ilp.fp.train
  base.ilp.fp.test <- dfs$base.ilp.fp.test
  expr.train <- dfs$expr.train
  expr.test <- dfs$expr.test

  for (gene in all.genes) {
    print(gene)
    y.train <- expr.train[, gene]
    y.test <- expr.test[, gene]
    y.train.norm <- (y.train-min(y.train))/(max(y.train)-min(y.train))
    y.test.norm <- (y.test-min(y.train))/(max(y.train)-min(y.train))

    # Train and log model.
    # ILP
    rr.ilp.model <- FitRidgeReg(base.ilp.train, y.train.norm)
    ilp.predictions <- predict(rr.ilp.model, 
                               s = rr.ilp.model$lambda.min, 
                               newx = base.ilp.test)[,1]
    ilp.perf <- GetPerformanceMetrics(y.test.norm, ilp.predictions)
    ilp.metrics <- c(ilp.perf$rsquared, ilp.perf$mse, ilp.perf$rmse)
    ilp.gene.perf <- c(gene, ilp.metrics)
    write(ilp.gene.perf, 
          ncolumns = length(ilp.gene.perf),
          append = TRUE, 
          file = paste0("../output/rr/ilp_results.txt")
          )

    # FP
    rr.fp.model <- FitRidgeReg(base.fp.train, y.train.norm)
    fp.predictions <- predict(rr.fp.model, 
                              s = rr.fp.model$lambda.min,
                              newx = base.fp.test)[,1]
    fp.perf <- GetPerformanceMetrics(y.test.norm, fp.predictions)
    fp.metrics <- c(fp.perf$rsquared, fp.perf$mse, fp.perf$rmse)
    fp.gene.perf <- c(gene, fp.metrics)
    write(fp.gene.perf, 
          ncolumns = length(fp.gene.perf),
          append = TRUE, 
          file = paste0("../output/rr/fp_results.txt")
          )

    # ILP + FP
    rr.ilp.fp.model <- FitRidgeReg(base.ilp.fp.train, y.train.norm)
    ilp.fp.predictions <- predict(rr.ilp.fp.model, 
                                  s = rr.ilp.fp.model$lambda.min,
                                  newx = base.ilp.fp.test)[,1]
    ilp.fp.perf <- GetPerformanceMetrics(y.test.norm, ilp.fp.predictions)
    ilp.fp.metrics <- c(ilp.fp.perf$rsquared, ilp.fp.perf$mse, ilp.fp.perf$rmse)
    ilp.fp.gene.perf <- c(gene, ilp.fp.metrics)
    write(ilp.fp.gene.perf, 
          ncolumns = length(ilp.fp.gene.perf),
          append = TRUE, 
          file = paste0("../output/rr/ilp_fp_results.txt")
          )

    # Perform merging
    # Averaging
    merged.preds <- rowMeans(cbind(ilp.predictions, fp.predictions))
    rf.avg <- GetPerformanceMetrics(y.test.norm, merged.preds)
    perf.metrics <- c(rf.avg$rsquared, rf.avg$mse, rf.avg$rmse)
    avg.gene.perf <- c(gene, perf.metrics) 
    write(avg.gene.perf, 
          ncolumns = length(avg.gene.perf),
          append = TRUE, 
          file = paste0("../output/rr/average.txt"))

    # NNLS
    ilp.y <- predict(rr.ilp.model, 
                     s = rr.ilp.model$lambda.min,
                     newx = base.ilp.train)[,1]
    fp.y <- predict(rr.fp.model, 
                    s = rr.fp.model$lambda.min,
                    newx = base.fp.train)[,1]
    weights <- perfCR(cbind(ilp.y, fp.y), y.train.norm)
    lm.preds <- rowSums(cbind(ilp.predictions, fp.predictions) %*% diag(weights))
    lrf.perf <- GetPerformanceMetrics(y.test.norm, lm.preds)
    lperf.metrics <- c(lrf.perf$rsquared, lrf.perf$mse, lrf.perf$rmse)
    lgene.perf <- c(gene, lperf.metrics)
    write(lgene.perf, 
          ncolumns = length(lgene.perf),
          append = TRUE, 
          file = paste0("../output/rr/stacked_nnls.txt"))
    write(weights, 
          ncolumns = length(weights),
          append = TRUE, 
          file = paste0("../output/rr/nnls_weights.txt"))
  }
}

PerformExperiments <- function(parallel_ = TRUE){
  if (parallel_){ 
    registerDoMC(cores = 50) 
  }
  # Load gene list.
  all.genes <- LoadGeneList()
  dfs <- LoadDatasets()
  PerformBaseExperiments(dfs, all.genes) 
}

PerformExperiments()
