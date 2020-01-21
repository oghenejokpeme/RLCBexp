#!/usr/bin/Rscript --vanilla
library(data.table)

# Get general peformance metrics. 
GetPerformanceMetrics <- function(y.actual, y.predicted){
  rsquared <- 1 - (sum((y.actual - y.predicted)^2) / 
                   sum((y.actual - mean(y.actual))^2))
  mse <- sum((y.actual - y.predicted)^2) / length(y.actual)
  rmse <- sqrt(sum((y.actual - y.predicted)^2) / length(y.actual))
    
  performance <- list(rsquared = rsquared,
                      mse = mse,
                      rmse = rmse)

  return(performance)
}

LoadGeneList <- function(){
  all.genes <- read.table(file = "../input/all_genes_ilp.txt",
                          header = FALSE)[,1]
  all.genes <- as.vector(all.genes)

  return(all.genes)
}

LoadDatasets <- function(){
  base.ilp.train <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "base_ilp_train.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  base.ilp.test <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "base_ilp_test.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  base.fp.train <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "base_fp_train.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  base.fp.test <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "base_fp_test.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  base.ilp.fp.train <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "base_ilp_fp_train.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  base.ilp.fp.test <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "base_ilp_fp_test.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  expr.train <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "expression_train.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
  expr.test <- data.matrix(data.frame(fread(paste0("../input/experiment_data/",
                                                  "expression_test.csv"
                                                  ), 
                                           header = TRUE),
                                     row.names = 1))
 
  df <- list(base.ilp.train = base.ilp.train,
             base.ilp.test = base.ilp.test,
             base.fp.train = base.fp.train,
             base.fp.test = base.fp.test,
             base.ilp.fp.train = base.ilp.fp.train,
             base.ilp.fp.test = base.ilp.fp.test,
             expr.train = expr.train,
             expr.test = expr.test)

  return(df)
}
