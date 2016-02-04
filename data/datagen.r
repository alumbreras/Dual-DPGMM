# Generate synthetic datasets
#
# author: Alberto Lumbreras
##################################
library(MASS)
library(ggplot2)
library(GGally)
library(gridExtra)
library(gridGraphics) 
cov <-  matrix(c(0.01,0.001,
                 0.001,0.01), ncol=2, nrow=2)


clear <- function(U=100, R=5){
  F <- 2
  # Assign roles uniformly
  z <- sort(rep(1:R, length=U))
  
  # Assign features
  A <- matrix(1, F, U)
  for(r in 1:R){
    x = cos(2*pi*r/R)
    y = sin(2*pi*r/R)
    cat(x,y,"\n")
    A[,z==r] <- t(mvrnorm(n=sum(z==r), mu=c(x,y), Sigma=cov))
  }
  
  # Assign behavior coefficients
  mu_br <- seq(-200,200, length.out=R)
  b <- rep(0, U)
  for(r in 1:R){
   b[z==r]<- rnorm(sum(z==r), mean=mu_br[r], sd=8)
  }
  
  list(A=A,
       b=b,
       z=z)
}


confused_features <- function(U=100, R=10){
  F <- 2

  # Assign roles uniformly
  z <- sort(rep(1:R, length=U))
  
  # Assign features
  A <- matrix(1, F, U)
  for(r in 1:(R-2)){
    x = 0.5*cos(2*pi*r/R)
    y = 0.5*sin(2*pi*r/R)
    cat(x,y,"\n")
    A[,z==r] <- t(mvrnorm(n=sum(z==r), mu=c(x,y), Sigma=cov))
  }
  
  for(r in (R-1):R){
    x = 0
    y = 0
    cat(x,y,"\n")
    A[,z==r] <- t(mvrnorm(n=sum(z==r), mu=c(x,y), Sigma=cov))
  }
  
  
  # Assign behavior coefficients
  mu_br <- seq(0,100, length.out=R)
  b <- rep(0, U)
  for(r in 1:R){
    b[z==r]<- rnorm(sum(z==r), mean=mu_br[r], sd=8)
  }
  
  list(A=A,
       b=b,
       z=z)
}


overlapped <- function(U=100, R=10){
  F <- 2
  # Assign roles uniformly
  z <- sort(rep(1:R, length=U))
  
  # Assign features
  A <- matrix(1, F, U)
  for(r in 1:R){
    x = 0.3*cos(2*pi*r/R)
    y = 0.3*sin(2*pi*r/R)
    cat(x,y,"\n")
    A[,z==r] <- t(mvrnorm(n=sum(z==r), mu=c(x,y), Sigma=cov))
  }
  
  # Assign behavior coefficients
  mu_br <- seq(-50,50, length.out=R)
  b <- rep(0, U)
  for(r in 1:R){
    b[z==r]<- rnorm(sum(z==r), mean=mu_br[r], sd=8)
  }
  
  list(A=A,
       b=b,
       z=z)
}

real <- function(){
  data(iris)
  idx <- sample(1:150, 50) # subselect 50 so that feature views does not dominate the inference
  iris <- iris[idx,]
  idx <- as.numeric(row.names(iris))
  iris <- iris[order(idx),]
  z <- as.numeric(iris$Species)
  U <- nrow(iris)
  A <- matrix(1, 3, U)
  
  # real 1
  A[1,] <- iris$Sepal.Width
  A[2,] <- iris$Petal.Length
  A[3,] <- iris$Petal.Width
  rownames(A) <- c("Sepal.Width", "Petal.Length", "Petal.Width")
  
  # too overlapped
  #b <- iris$Sepal.Length
  #b <- 100*(b-min(b))/(max(b)-min(b)) # to a 0,100 range

  # Assign behavior coefficients
  mu_br <- seq(0,100, length.out=3)
  b <- rep(0, U)
  for(r in 1:3){
    b[z==r]<- rnorm(sum(z==r), mean=mu_br[r], sd=8)
  }
  
  list(A=A,
       b=b,
       z=z)
}


#############################################
# Generate participations and thread lengths
#############################################
# from this single dataset we will extract a subset of users and threads according to the scenario.
participations <- function(b=c(1,2,3,4,5), b0=0, s_y=1/5, T=10){
  # Warning: in Iris, b is much smaller and therefore s_y should be greater
  U <- length(b)
  
  # Participations matrix
  P <- matrix(sample(c(0,1), U*T, replace=TRUE, prob=c(0.5,0.5)), nrow=U, ncol=T)
  P <- apply(P,2, function(x){x/sum(x)}) # normalize so that noise s_y has the same effect in every thread
  P <- P[, colSums(is.na(P)) != nrow(P)] # drop empty threads
  
  # Length of threads
  T <- dim(P)[2]
  P_ <- rbind(rep(1,T), P)
  b_ <- c(b0,b)
    
  means <- t(P_)%*%b_
  cat("\ns_y*",s_y, sqrt(1/s_y))
  y <- sapply(1:T, function(i) {rnorm(1, mean=means[i], sd=sqrt(1/s_y))})
  
  list(P=P,
       y=y)
}



##########################################################
##########################################################
save.data <- function(A, b, z){
  # Save dataset to files
  U <- dim(A)[2]
  df <- data.frame(matrix(NA, U, 4))
  names(df) <- c("z", "f1", "f2", "b")
  df$z <- z-1 # start from 0
  df$b <- b
  df$f1 <- A[1,]
  df$f2 <- A[2,]
  write.table(df, file=paste0("data_users_", U, ".csv"), sep= "\t", row.names = FALSE)  
}

load.data <- function(filepath){
  df <- read.table(filepath, sep='\t', header=TRUE)
  df
}

save.real.data <- function(A, b, z){
  # Save dataset to files
  U <- dim(A)[2]
  df <- data.frame(matrix(NA, U, 5))
  names(df) <- c("z", "f1", "f2", "f3", "b")
  df$z <- z-1 # start from 0
  df$b <- b
  df$f1 <- A[1,]
  df$f2 <- A[2,]
  df$f3 <- A[3,]
  write.table(df, file=paste0("data_users_", U, ".csv"), sep= "\t", row.names = FALSE)  
}


save.participations <- function(P, y, suffix="train"){
  # Save participation matrix and thread lengths
  # one big participations file for every scenario and every U
  U <- dim(P)[1]
  write.table(as.data.frame(P), file=paste0(suffix, "_participations_", U, ".csv"), sep= "\t", row.names = FALSE)
  write.table(as.data.frame(y), file=paste0(suffix, "_lengths_", U, ".csv"), sep= "\t", row.names = FALSE) 
}

plot.data <- function(A, b, z){    
  p1 <- qplot(A[1,], A[2,]) +
        geom_point(aes(colour = factor(z), size=b)) +
        xlab('feature 1') +
        ylab('feature 2') +
        theme(text = element_text(size = 15),
              panel.background = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1,
              axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill = "white"), 
              legend.position="none")

  p2 <- qplot(1:length(b), b) +
        geom_point(aes(colour = factor(z))) +
        xlab('users') +
        ylab('b') +
        theme(text = element_text(size = 15),
              panel.background = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1,
              axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill = "white"), 
              legend.position="none")
  
  grid.arrange(p1, p2, ncol=2)
  g <- arrangeGrob(p1, p2, ncol=2)
  g
}

plot.real.data <- function(A, b, z){
  
  A.df <- as.data.frame(cbind(t(A), z))
  A.df$z <- as.factor(A.df$z)
  theme_set(theme_bw())
  p1 <- ggpairs(A.df, columns= 1:3,
          upper = list(continuous = "points"),
          lower = list(continuous = "points"),
          diag = list(continuous = "density"),
          colours='z')+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio=1)
  
  p2 <- qplot(1:length(b), b) +
    geom_point(aes(colour = factor(z))) +
    xlab('users') +
    ylab('b') +
    theme(text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio=1,
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = "white"), 
          legend.position="none")

  print(p1)
  g <- grid.grab() 
  grid.arrange(g, p2, widths=c(0.5,0.5))  
}

##########################################################
# Main
##########################################################
if (FALSE){
  generators <- c(clear, confused_features, overlapped)
  gen <- generators[[2]] # chose one of the three scenarios
  
  # Generate and plot the data
  i = 50 # number of users
  # user features and behaviors
  data <- gen(U=i, R=5)
  g<- plot.data(data$A, data$b, data$z)
  save.data(data$A, data$b, data$z)
  
  # user participations and thread lengths
  # (a big synthetic forum of U users)
  # make T big enough so that training in different experiments get different training instances
  parts <- participations(b=data$b, T=10000)
  save.participations(parts$P,parts$y, suffix="train")
  
  # test set
  parts <- participations(b=data$b, T=1000)
  save.participations(parts$P,parts$y, suffix="test")
}

if (TRUE){
  # Generate a dataset that reproduces the Iris dataset structure
  # s_y adapted to the size of b so that the noise is reasonable
  data <- real()
  plot.real.data(data$A, data$b, data$z)
  save.real.data(data$A, data$b, data$z)
  
  # user participations and thread lengths
  #parts <- participations(b=data$b, T=1000, s_y=1000)
  parts <- participations(b=data$b, T=1000)
  save.participations(parts$P,parts$y, suffix="train")
  
  # test set
  #parts <- participations(b=data$b, T=100, s_y=1000)
  parts <- participations(b=data$b, T=100)
  save.participations(parts$P,parts$y, suffix="test")
  
  # sanity checks
  ###############################################"
  par(mfrow=c(3,1))
  b <- data$b
  P <- parts$P
  y <- parts$y
  U <- dim(P)[1]
  alphaI <- diag(U)*0.01
  
  T <- dim(P)[2]
  P_ <- rbind(rep(1,T), P)
  
  # plot lengths with b and no noise
  b_ <- c(0,b)  
  y_ideal <- t(P_)%*%b_
  idx <- order(y_ideal) # take these lengths idx as reference
  plot(y_ideal[idx])
  title("noise free")
  #plot(b)
  
  # plot the generated lengths
  plot(y[idx])
  title("with noise (real data)")
  
  # what predictiorepeat with mle
  b_mle <- (solve(P%*%t(P)+alphaI)%*%P)%*%y
  b_ <- c(0,b_mle)
  y_mle <- t(P_)%*%b_
  plot(y_mle[idx])
  title("reproduce with estimated coefficients")
}