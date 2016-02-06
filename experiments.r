# The experiments are called form this script
# It loads the data, call the Gibbs sampler,
# and saves the traces in files for further analysis 
# and evaluation of thread length predictions in the test set
# author: Alberto Lumbreras
#######################################################
experiment <- function(nthreads.train){
    source('gibbs_dual.r')
  
    # Load data
    ######################################
    # Choose one of the datasers
    dataset = 'confused_features'
    dataset = 'iris'
    dataset = 'clear'
    dataset = 'overlapped'
    
    data.dir <- paste0('./data/', dataset, '/')
    df <- read.table(paste0(data.dir, 'data_users_50.csv'), sep='\t', header=TRUE)
    z_init <- df$z+1
    if(dataset=='iris'){
      A <- t(df[,2:4])
    }else{
      A <- t(df[,2:3])  
    }
    B <- t(df$b)
  
  
    y <- read.table(paste0(data.dir, 'train_lengths_50.csv'), sep='\t', header=TRUE)
    P <- read.table(paste0(data.dir, 'train_participations_50.csv'), sep='\t', header=TRUE)
    nthreads <- dim(P)[2]
    
    # Select a random subset of threads for training the model
    idx.train <- sample(nthreads, nthreads.train)
    y <- y[idx.train,]
    P <- P[,idx.train]
    
    # Initialization of cluster assignments
    ########################################
    # Non recommended.
    z_init <- rep(1, dim(B)[2])
    
    # Recommended (bit not a good initialization technique either!).
    z_init <- sample(10, dim(B)[2], replace = TRUE)
    
    # A la k-means (more clever)
    z_init <- kmeans(t(A), 10)$cluster
    
    # guarantees that there are no gaps!
    # http://stackoverflow.com/questions/35141155/create-n-random-integers-with-no-gaps
    tb <- table(z_init); 
    z_init <- rep(seq(tb),tb)
    
    # run !
    #############
    res <- gibbs(A, B, P, y, z_init=z_init, iters=30000)
    
    # Save traces to files
    #######################
    # do not overwrite old experiments
    experiment.path <- file.path("out", dataset, paste0("DP_threads_",nthreads.train))
    i <- 1
    while(TRUE){
      if(dir.exists(paste0(experiment.path, '-', i))){
          i <- i+1
      }
      else{
        experiment.path <- paste0(experiment.path, '-', i)
        break
      }
    }
    traces.path <- file.path(experiment.path, "traces")
    dir.create(file.path(traces.path), recursive=TRUE)
    for(i in 1:length(res)){
      tr.name <- paste0(names(res)[i], '.trc')
      write.matrix(as.data.frame(res[[i]]), file=file.path(traces.path,tr.name))
    }
}

# Run experiments for different number of threads in training
#for(nthreads.train in seq(90,100, by=10)){
#  experiment(nthreads.train)
#}

library(doParallel)
library(foreach)
library(parallel)

# See how to use a progress bar with doSNOW
#http://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar

if(TRUE){
  cl<-makeCluster(6, outfile="")
  registerDoParallel(cl)
  pck = c('abind', 'MASS', 'mvtnorm', 'mixtools', 'coda')
  foreach(i=seq(10,100, by=10), .packages = pck)%dopar%experiment(i)
  stopCluster(cl)
}
