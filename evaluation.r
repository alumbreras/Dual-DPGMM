# For experiment in the dataset folder:
# - Compute the ARI using the true clusters (in the data folder) and the z samples
# - Compute the negative loglikelihood of the true lengths given some prediction samples.
# store the results in a csv file. 
# author: Alberto Lumbreras

library(mclust) # contains an ARI function

loglike <- function(y_test, P_test, traces_coefficients, traces_intercept, traces_noise_inv){
  n_sample = dim(traces_coefficients)[1]
  loglike = 0
  locs <-  t(P_test) %*% t(traces_coefficients) # y_pred without intercept
  locs <- t(apply(locs, 1, function(x) x+t(traces_intercept))) # add intercept 
  sd_noise = 1/sqrt(traces_noise_inv)
  for(j in 1:n_sample){
    loglike <- loglike + sum(dnorm(y_test, 
                                   mean = locs[,j],
                                   sd = sd_noise[j],
                                   log=TRUE)
    )
  }
  1.0*loglike/(n_sample*length(y_test))
}

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

burnin <- 7000
# Load test set
dataset <- 'iris'
P_test <- read.csv(file.path("./data/", dataset, "/test_participations_50.csv"), sep='\t')
y_test <- read.csv(file.path("./data/", dataset, "/test_lengths_50.csv"), sep='\t')$y
z_true <- read.csv(file.path("./data/", dataset, "/data_users_50.csv"), sep='\t')$z

# Get prediction samples for P_test and compute negative loglikelihood by comparing the predictions
# with the y_test
experiments <- list.files(path = file.path("out", dataset), recursive=F, full.names = T)
for (i in 1:length(experiments)){
  traces.path <- file.path(experiments[i], 'traces')
  traces_coefficients <- read.csv(file.path(traces.path, 'traces.coefficients.trc'), sep='')[-c(1:burnin),]
  traces_intercept <- read.csv(file.path(traces.path, 'traces.intercept.trc'), sep='')[-c(1:burnin),]
  traces_noise_inv <- read.csv(file.path(traces.path, 'traces.noise_inv.trc'), sep='')[-c(1:burnin),]
  traces_coefficients <- as.matrix(traces_coefficients)
  traces_noise_inv <- as.matrix(traces_noise_inv)

  # Negative loglikelihood
  negloglike <- -loglike(y_test, P_test, traces_coefficients, traces_intercept, traces_noise_inv)
  cat("\n", experiments[i], "negloglike:", negloglike)
  
  # Adjusted Rand Index
  traces_z <- read.csv(file.path(traces.path, 'traces.z.trc'), sep='')[-c(1:burnin),]
  traces_z <- as.matrix(traces_z)
  pairwise <- pairwise_posterior(traces_z[1:100,])
  ls_z <- least_squares_clustering(traces_z[1:100,], pairwise)
  ari = adjustedRandIndex(z_true, ls_z)
  cat("\n", experiments[i], "ARI:", ari)
  
  # Plot pairwise clustering and predictions
  image(1:dim(pairwise)[1], 1:dim(pairwise)[1], pairwise, col= gray((0:128)/128), 
        xlab="", ylab= "", asp=1, axes=FALSE)
  axis(1,seq(0,50, by=10),seq(0,50, by=10), pos=0.5)
  axis(2,seq(0,50, by=10),seq(0,50, by=10), pos=0.5)
  title(experiments[i])

  # Save to file
  model <- strsplit(strsplit(experiments[i], "/")[[1]][3], '_')[[1]][1]
  nthreads <- strsplit(strsplit(strsplit(experiments[i], "/")[[1]][3], '_')[[1]][3], "-")[[1]][1]
  write.table(t(c(model, nthreads, negloglike, ari)), 
              file="results.csv", sep='\t', col.names=FALSE, row.names=FALSE, append=TRUE)
  
  
}