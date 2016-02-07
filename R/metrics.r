library(mclust) # contains an ARI function

mse <- function(y_test, P_test, traces_coefficients, traces_intercept, traces_noise_inv){
  n_sample <- dim(traces_coefficients)[1]
  error <- 0
  locs <-  t(P_test) %*% t(traces_coefficients) # y_pred without intercept
  locs <- t(apply(locs, 1, function(x) x+t(traces_intercept))) # add intercept 
  for(j in 1:n_sample){
    error <- error + sum((y_test-locs[,j])^2)
  }
  error/(n_sample*length(y_test))
}


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
