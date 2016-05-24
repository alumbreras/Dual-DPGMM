# For experiment in the dataset folder:
# - Compute the ARI using the true clusters (in the data folder) and the z samples
# - Compute the negative loglikelihood of the true lengths given some prediction samples.
# store the results in a csv file. 
# author: Alberto Lumbreras

source('R/metrics.r')

pairwise_posterior <- function(traces_z){
  n_samples <- dim(traces_z)[1]
  n_dim <- dim(traces_z)[2]
  pairwise <- matrix(0, n_dim, n_dim)
  for(i in 1:n_samples){
    pairwise <- pairwise + apply(traces_z[i,,drop=FALSE], 2, function(x) as.numeric(x==traces_z[i,,drop=FALSE]))
  }
  return(pairwise/n_samples)
}

least_squares_clustering <- function(traces_z, pairwise){
  best_sample <- traces_z[1,]
  lse_min <- 1000000000
  for(i in 1:nrow(traces_z)){
    lse <- sum((apply(traces_z[i,,drop=FALSE], 2, function(x) x==traces_z[i,,drop=FALSE]) - pairwise)^2)
    if (lse < lse_min){
      LSE_min <-  lse
      best_sample <- traces_z[i,]
    }
  }
  return(best_sample)
}


evaluate <- function(experiment.path, burnin){
  # Get prediction samples for P_test and compute negative loglikelihood by comparing the predictions
  # with the y_tes
  
  P_test <- read.csv(file.path("./data/", dataset, "/test_participations_50.csv"), sep='\t')
  y_test <- read.csv(file.path("./data/", dataset, "/test_lengths_50.csv"), sep='\t')$y
  z_true <- read.csv(file.path("./data/", dataset, "/data_users_50.csv"), sep='\t')$z
  
  traces.path <- file.path(experiment.path, 'traces')
  traces_coefficients <- read.csv(file.path(traces.path, 'traces.coefficients.trc'), sep='')[-c(1:burnin),]
  traces_intercept <- read.csv(file.path(traces.path, 'traces.intercept.trc'), sep='')[-c(1:burnin),]
  traces_noise_inv <- read.csv(file.path(traces.path, 'traces.noise_inv.trc'), sep='')[-c(1:burnin),]
  traces_coefficients <- as.matrix(traces_coefficients)
  traces_noise_inv <- as.matrix(traces_noise_inv)
  
  # Negative loglikelihood
  negloglike <- -loglike(y_test, P_test, traces_coefficients, traces_intercept, traces_noise_inv)
  cat("\n\n", experiments[i], "negloglike:", negloglike)
  
  # Mean squared error
  meanse <- mse(y_test, P_test, traces_coefficients, traces_intercept, traces_noise_inv)
  cat("\n", experiments[i], "mse:", meanse)
  
  # Adjusted Rand Index
  traces_z <- read.csv(file.path(traces.path, 'traces.z.trc'), sep='')[-c(1:burnin),]
  traces_z <- as.matrix(traces_z)
  pairwise <- pairwise_posterior(traces_z[5000:5200,])
  ls_z <- least_squares_clustering(traces_z[5000:5200,], pairwise)
  ari = adjustedRandIndex(z_true, ls_z)
  cat("\n", experiments[i], "ARI:", ari)
  
  # Plot pairwise clustering and predictions
  image(1:dim(pairwise)[1], 1:dim(pairwise)[1], pairwise, col= gray((0:128)/128), 
        xlab="", ylab= "", asp=1, axes=FALSE)
  axis(1,seq(0,50, by=10),seq(0,50, by=10), pos=0.5)
  axis(2,seq(0,50, by=10),seq(0,50, by=10), pos=0.5)
  title(experiments[i])
  
  # Save to file
  model <- strsplit(strsplit(experiment.path, "/")[[1]][3], '_')[[1]][1]
  nthreads <- strsplit(strsplit(strsplit(experiment.path, "/")[[1]][3], '_')[[1]][3], "-")[[1]][1]
  write.table(t(c(experiment.path, model, nthreads, negloglike, ari)), 
              file=file.path("out", dataset, "results.csv"), sep='\t', col.names=FALSE, row.names=FALSE, append=TRUE)
  
}

###################################
library(doParallel)
library(foreach)
library(parallel)

par(mfrow=c(3,3))
burnin <- 10000

# Load test set
dataset <- 'iris'
dataset <- 'overlapped'
dataset <- 'clear'
dataset <- 'confused_features' #mediamining
dataset <- "disagreement"
dataset <- "agreement"

# Experiments that have not been yet evaluated 
experiments <- list.dirs(path = file.path("out", dataset), recursive=F, full.names = T)
experiments.evaluated <- read.csv(file.path('out', dataset, 'results.csv'), sep='\t', head=FALSE)[,1]
experiments <- experiments[! experiments %in% experiments.evaluated]

if(length(experiments)>0){
  ncores <- detectCores() - 2
  cl<-makeCluster(ncores, outfile="", port=11439)
  registerDoParallel(cl)
  pck = c('mclust')
  foreach(i=1:length(experiments), .packages = pck)%dopar%evaluate(experiments[i], burnin)
  stopCluster(cl)
}else{
  cat("\nThere aren't new experiments to evaluate")
}
