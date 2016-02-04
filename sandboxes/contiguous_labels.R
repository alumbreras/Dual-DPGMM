z_init <- c(3,2,1, 3,3,7,9)
sort(z_init)

set.seed(1)
K <- 10 # initial number of clusters
N <- 20 # number of data points
z_init <- sample(K,N, replace=TRUE) # initial assignments


idx <- order(z_init)
for (i in 2:length(z_init)){
  if(z_init[idx[i]] > z_init[idx[i-1]]){
    z_init[idx[i]] <- z_init[idx[i-1]]+1
  }
  else{
    z_init[idx[i]] <- z_init[idx[i-1]]  
  }
}

# http://stackoverflow.com/questions/35141155/create-n-random-integers-with-no-gaps
z_init <- c(3,2,1, 3,3,7,9)
y=table(z_init); 
rep(seq(y),y)