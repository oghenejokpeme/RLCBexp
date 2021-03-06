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

    # Train and log model.
    # ILP
    rf.ilp.model <- FitRandomForest(base.ilp.train, y.train)
    ilp.predictions <- predict(rf.ilp.model, base.ilp.test)$predictions
    ilp.perf <- GetPerformanceMetrics(y.test, ilp.predictions)
    ilp.metrics <- c(ilp.perf$rsquared, ilp.perf$mse, ilp.perf$rmse)
    ilp.gene.perf <- c(gene, ilp.metrics)
    write(ilp.gene.perf, 
          ncolumns = length(ilp.gene.perf),
          append = TRUE, 
          file = paste0("../output/rf/ilp_results.txt")
          )

    # FP
    rf.fp.model <- FitRandomForest(base.fp.train, y.train)
    fp.predictions <- predict(rf.fp.model, base.fp.test)$predictions
    fp.perf <- GetPerformanceMetrics(y.test, fp.predictions)
    fp.metrics <- c(fp.perf$rsquared, fp.perf$mse, fp.perf$rmse)
    fp.gene.perf <- c(gene, fp.metrics)
    write(fp.gene.perf, 
          ncolumns = length(fp.gene.perf),
          append = TRUE, 
          file = paste0("../output/rf/fp_results.txt")
          )

    # ILP + FP
    rf.ilp.fp.model <- FitRandomForest(base.ilp.fp.train, y.train)
    ilp.fp.predictions <- predict(rf.ilp.fp.model, base.ilp.fp.test)$predictions
    ilp.fp.perf <- GetPerformanceMetrics(y.test, ilp.fp.predictions)
    ilp.fp.metrics <- c(ilp.fp.perf$rsquared, ilp.fp.perf$mse, ilp.fp.perf$rmse)
    ilp.fp.gene.perf <- c(gene, ilp.fp.metrics)
    write(ilp.fp.gene.perf, 
          ncolumns = length(ilp.fp.gene.perf),
          append = TRUE, 
          file = paste0("../output/rf/ilp_fp_results.txt")
          )

    # Perform merging
    # Averaging
    merged.preds <- rowMeans(cbind(ilp.predictions, fp.predictions))
    rf.avg <- GetPerformanceMetrics(y.test, merged.preds)
    perf.metrics <- c(rf.avg$rsquared, rf.avg$mse, rf.avg$rmse)
    avg.gene.perf <- c(gene, perf.metrics) 
    write(avg.gene.perf, 
          ncolumns = length(avg.gene.perf),
          append = TRUE, 
          file = paste0("../output/rf/average.txt"))

  }
}

PerformExperiments <- function(parallel_ = TRUE){
  if (parallel_){ 
    registerDoMC(cores = 24) 
  }
  # Load gene list.
  all.genes <- LoadGeneList()
  dfs <- LoadDatasets()
  PerformBaseExperiments(dfs, all.genes) 
}

PerformExperiments()
