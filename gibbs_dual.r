# Dual clustering of users
# author: Alberto Lumbreras
###############################
library(abind)
library(coda)
library(mixtools)

source('conditionals_dual.r')
source('ars_alpha.r')
source('ars_beta.r')


# For debugging, this is useful:
if(FALSE){
  plot(mcmc(traces.b[complete.cases(traces.b),]))
  plot(mcmc(1/traces.noise_inv[complete.cases(traces.noise_inv),]))
  plot(mcmc(1/traces.intercept[complete.cases(traces.intercept),]))
  plot(mcmc(traces.coefficients[complete.cases(traces.coefficients),1:10]))
}

PLOTS <- FALSE
PLOTS.FREQUENCY <- 10 # iterations between two plots
data(iris)
data(geyser)

gibbs <- function(A, B, P, y, z_init, iters=100){
  # Arguments:
  #   A: matrix with user features (one column per user)
  #   B: matrix with user coefficients (one column per user and only one row)
  #   P: Participation matrix (one column per thread, one row per user)
  #   y: Thread lengths (as much as columns in P)
  #   z_init: initial assignments (one per user)
  
  # For a Gaussian Mixture Model with a fixed number of clusters
  # set alpha = 0 and don't sample it
  # (alpha controls the probability of opening a new cluster)
  
  ######################################################
  # Choose which view to activate
  # If FEATURES only, the model is a DP-GMM
  # If BEHAVIORS only, the model is a DP - bayesian linear regression 
  # (coefficients assumed to come from clusters)
  # If FEATURES and BEHAVIORS, dual model presented in Lumbreras et al 2016.
  FEATURES <- TRUE
  BEHAVIORS <- TRUE

  threads.idx <- order(y) # for plotting
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  N <- dim(A)[2]
  
  P <- as.matrix(P)
  y <- as.matrix(y)
  
  # traces
  traces.z <- matrix(NA, iters, N)
  traces.coefficients <- matrix(NA, iters, N)
  traces.intercept <- matrix(NA, iters, 3)
  traces.alpha <- matrix(NA, iters, 1)
  traces.a <- matrix(NA, iters, 5)
  traces.b <- matrix(NA, iters, 5)
  traces.noise_inv <- matrix(NA, iters,1)
  
  colnames(traces.a) <- c("n_clusters", "mean_mu_a0", "det(R_a0)", "det(W_a0)", "beta_a")
  colnames(traces.b) <- c("n_clusters", "mean_mu_b0", "det(R_b0)", "det(W_b0)", "beta_b")
  colnames(traces.alpha) <- c("alpha")
  colnames(traces.noise_inv) <- c("noise_inv")
  colnames(traces.intercept) <- c("intercept", "mu_intercept", "s_intercept")
  
  # data mean and covariance
  # used to centers priors on the data
  Sigma_a <- cov(t(A))
  Lambda_a <- solve(Sigma_a)
  mu_a <- as.matrix(rowMeans(A))

  Sigma_b <- cov(t(B))
  Lambda_b <- solve(Sigma_b)
  mu_b <- as.matrix(rowMeans(B))
  
  variance_y <- c(var(y))
  ##########################
  # Init state
  ##########################
  # init common variables
  alpha = 2
  z <- z_init 
  
  # init feature variables
  D <- dim(A)[1]
  mu_a0 <- mu_a
  R_a0 <- solve(Sigma_a)
  W_a0 <- Sigma_a
  beta_a0 <- dim(A)[1]
  S_ar <- array(NA, dim=c(D, D, 100))
  mu_ar <- matrix(NA, D, 100)
  for (k in 1:length(unique(z))){
    mu_ar[,k] <- rowMeans(A[,z==k, drop=FALSE])
  }
  
  # init behavior variables
  D <- dim(B)[1]
  intercept <- 0
  noise_inv <- 1/variance_y
  mu_b0 <- mu_b
  R_b0 <- solve(Sigma_b)
  W_b0 <- Sigma_b
  beta_b0 <- dim(B)[1]
  S_br <- array(NA, dim=c(D, D, 100))
  mu_br <- matrix(NA, D, 100)
  print(k)
  for (k in 1:length(unique(z))){
    mu_br[,k] <- rowMeans(B[,z==k, drop=FALSE])
  }
  mu_intercept <- 0
  
  
  par(mfrow=c(1,3))
  
  for(i in 1:iters){
    
    cat('\n', i)
    
    # Active components
    K <- length(unique(z))
    
    ##############################################################
    # Sample variables in features view
    ##############################################################
    
    if(FEATURES){
      # Sample component parameters
      for (k in 1:length(unique(z))){
        mask <- (z==k)
        S_ar[,,k] <- sample_S_ar(A[,mask, drop=FALSE], mu_ar[,k, drop=FALSE], beta_a0, W_a0)
        mu_ar[,k] <- sample_mu_ar(A[,mask, drop=FALSE], S_ar[,,k], mu_a0, R_a0)
      }
      
      # Sample hyperpriors
      mu_a0 <- sample_mu_a0(Lambda_a, mu_a, mu_ar[,1:K, drop=FALSE], R_a0)
      R_a0 <- sample_R_a0(Sigma_a, mu_ar[,1:K, drop=FALSE], mu_a0) 
      W_a0 <- sample_W_a0(Lambda_a, S_ar[,,1:K, drop=FALSE], beta_a0)
      beta_a0 <- sample_beta_a0(S_ar[,,1:K, drop=FALSE], W_a0, init=beta_a0)
      traces.a[i,] <- c(length(unique(z)), mean(mu_a0), det(R_a0), det(W_a0), beta_a0)
      
      if(PLOTS && (i%%PLOTS.FREQUENCY==0)){
        plot(t(A[c(1,2),]), col=palette()[z+1])
        title(paste("Features", i))
        for(k in 1:K){
          ellipse(mu=mu_a0[c(1,2)], sigma=W_a0[c(1,2),c(1,2)], alpha = .25, lwd=5, npoints = 250, col="black")
          ellipse(mu=mu_ar[c(1,2),k], sigma=solve(S_ar[c(1,2),c(1,2),k]),  alpha = .25, lwd=3, npoints = 250, col=palette()[k+1])
        }
      }
    }
    
    ##############################################################
    # Sample variables in behaviors view
    ##############################################################
    if(BEHAVIORS){
        # Sample component parameters
        for (k in 1:length(unique(z))){
          mask <- (z==k)
          S_br[,,k] <- sample_S_ar(B[,mask, drop=FALSE], mu_br[,k, drop=FALSE], beta_b0, W_b0)
          #S_br[,,k] <- sample_S_ar(B[,mask, drop=FALSE], mu_br[,k, drop=FALSE], beta_b0, as.matrix(1))
          mu_br[,k] <- sample_mu_ar(B[,mask, drop=FALSE], S_br[,,k], mu_b0, R_b0)
          #if(abs(mu_br[,k]) < 0.0001){stop("bad mu_br")}
        }
        # Intercept has its own "component" parameters
        s_intercept <- sample_S_ar(as.matrix(intercept), as.matrix(mu_intercept), 1, 1)
        mu_intercept <- sample_mu_ar(as.matrix(intercept), as.matrix(s_intercept), 0, 1)
        
        # Sample hyperpriors
        mu_b0 <- sample_mu_a0(Lambda_b, mu_b, mu_br[,1:K, drop=FALSE], R_b0)
        R_b0 <- sample_R_a0(Sigma_b, mu_br[,1:K, drop=FALSE], mu_b0) 
        W_b0 <- as.matrix(1) #sample_W_a0(Lambda_b, S_br[,,1:K, drop=FALSE], beta_b0)
        beta_b0 <- sample_beta_a0(S_br[,,1:K, drop=FALSE], W_b0, init=beta_b0)
        
        # This part is the only difference between the behavior view and the feature view
        # once B is sampled, the rest of the sampler sees B as regular features 
        # Sample latent coefficients
        # Sample intercept from a gaussian prior with mean 0 and flat variance
        #noise_inv <- 100
        #B <- sample_b(P, y, z, intercept, mu_br[,1:K, drop=FALSE], S_br[,,1:K, drop=FALSE], s_y=noise_inv)
        B <- tryCatch(sample_b(P, y, z, intercept, mu_br[,1:K, drop=FALSE], S_br[,,1:K, drop=FALSE], s_y=noise_inv),
                      error=function(e){
                        cat(paste("sample_b(): ", e))
                        return(B)
                      }
            )

        intercept <- sample_intercept(P, y, z, B, mu_intercept, s_intercept, s_y=noise_inv)
        noise_inv <- sample_noise_inv(P, y, B, variance_y)
        #if(1/noise_inv>100) stop("bad noise_inv")
        cat("\nintercept:", intercept, "\tnoise:", 1/noise_inv, "W_b0:", W_b0)
  
        
        traces.b[i,] <- c(length(unique(z)), mean(mu_b0), det(R_b0), det(W_b0), beta_b0)
        traces.intercept[i,] <- c(intercept, mu_intercept, s_intercept)
        traces.coefficients[i,] <- c(B)
        traces.noise_inv[i,] <- noise_inv
        
        if(PLOTS && (i%%PLOTS.FREQUENCY==0)){
          plot(t(B), col=palette()[z+1])
          title("Behaviors")
          
          plot(y[threads.idx,], type='l', col="red", lwd=3)
          lines(t(P[,threads.idx])%*%t(B)+intercept)
          title("Predictions")
        }

    }
    
    #############################################################
    # Sample common parameters
    #############################################################
    # Sample assignments and concentration parameter
    if(TRUE){
      for (n in 1:N){
        if(BEHAVIORS && FEATURES){
          res <- sample_z_dual(n, A, B, alpha, z, 
                          mu_ar, S_ar, mu_a0, R_a0,  beta_a0, W_a0,
                          mu_br, S_br, mu_b0, R_b0, beta_b0, W_b0)
          z <- res$z
          mu_ar <- res$mu_ar
          S_ar <- res$S_ar
          mu_br <- res$mu_br
          S_br <- res$S_br
        }
        if(FEATURES && !BEHAVIORS){
          res <- sample_z(n, A, alpha, z, 
                          mu_ar, S_ar, mu_a0, R_a0,  beta_a0, W_a0)
          z <- res$z
          mu_ar <- res$mu_ar
          S_ar <- res$S_ar
        }
        if(BEHAVIORS && !FEATURES){
          res <- sample_z(n, B, alpha, z, 
                          mu_br, S_br, mu_b0, R_b0,  beta_b0, W_b0)
          z <- res$z
          mu_br <- res$mu_ar
          S_br <- res$S_ar
        }
        # Assert
        if(max(z) != length(unique(z))){
          cat('\nuser: ', n)
          cat(tabulate(z))
          cat(z[u])
          stop("Some cluster is not being used")
        }
        #z <- rep(1, length(z)) # debug

      }
      
      K <- length(unique(z))
      alpha <- sample_alpha(K = K, U = N)
      
      traces.z[i,] <- z
      traces.alpha[i,] <- alpha
      
      cat('\n',tabulate(z))
    }
    

  }
  
  res <- list(traces.a=traces.a,
              traces.b=traces.b,
              traces.z=traces.z,
              traces.alpha=traces.alpha,
              traces.coefficients=traces.coefficients,
              traces.intercept=traces.intercept, 
              traces.noise_inv=traces.noise_inv)   
  return(res)
}


if(FALSE){
  
  # Load data
  ######################################
  # Choose one of the datasers
  dataset = 'confused_features'
  dataset = 'overlapped'
  dataset = 'iris'

  # Number of threads used for training
  nthreads.train <- 100
  
  # Load dataset
  data.dir <- paste0('./data/', dataset, '/')
  df <- read.table(paste0(data.dir, 'data_users_50.csv'), sep='\t', header=TRUE)
  z_init <- df$z+1

  if(dataset=='iris'){
    A <- t(df[,2:4])
  }
  else{
    A <- t(df[,2:3])  
  }
  
  B <- t(df$b)
  y <- read.table(paste0(data.dir, 'train_lengths_50.csv'), sep='\t', header=TRUE)
  P <- read.table(paste0(data.dir, 'train_participations_50.csv'), sep='\t', header=TRUE)

  nthreads <- dim(P)[2]
  idx.train <- sample(nthreads, nthreads.train)
  y <- y[idx.train,]
  P <- P[,idx.train]
  
  # Initialization of cluster assignments
  ########################################
  # Non recommended.
  z_init <- rep(1, dim(B)[2])
  
  # Recommended (bit not a good initialization technique either!).
  z_init <- sample(10, dim(B)[2], replace = TRUE)
  
  # guarantees that there are no gaps!
  # http://stackoverflow.com/questions/35141155/create-n-random-integers-with-no-gaps
  tb <- table(z_init); 
  z_init <- rep(seq(tb),tb)
  
  # run !
  #############
  res <- gibbs(A, B, P, y, z_init=z_init, iters=1000)
  
  # Plot traces
  #######################
  traces.a <- res$traces.a
  traces.b <- res$traces.b
  traces.z <- res$traces.z
  chains <- mcmc(traces.a)
  plot(chains)
  print(effectiveSize(chains))
  chains <- mcmc(traces.b)
  plot(chains)
  print(effectiveSize(chains))
}