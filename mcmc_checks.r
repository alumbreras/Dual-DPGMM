# Checks MCMC traces of a given experiment
#
# author: Alberto Luumbreras
##########################
# http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf
library(coda)
library(ggplot2)
library(reshape2)

traces.dir <- "./out/iris/DP_threads_30-1"
burn <- 1000 # burned samples

# Read trace files
files <- list.files(path = traces.dir, pattern = ".trc", recursive=T, full.names = T)
traces <- NULL
for(i in 1:length(files)){
  filename <- tail(strsplit(files[i], "/", fixed=T)[[1]], 1)
  varname <- strsplit(filename, ".", fixed=T)[[1]][2]
  traces[[varname]] <- read.csv(files[i], sep="")
}

nclusters <- apply(traces[['z']], 1, function(x) {length(unique(x))})

# Plot traces, densities, autocorrelations and Geweke's test for each variabl
# but z and coefficients
traces[['a']]<- traces[['a']][-1]
traces[['b']]<- traces[['b']][-1]
df.traces <- cbind(traces[['alpha']], traces[['intercept']], traces[['noise_inv']], traces[['a']], traces[['b']])
par(mfrow=c(3,4))
for (i in 1:ncol(df.traces)){
  varname <- names(df.traces)[i]
  df <- data.frame(df.traces[,i])
  names(df) <- varname
  chain <- mcmc(df)
  
  # Chain and density plots
  plot(chain, auto.layout = FALSE)
  
  # Autocorrelation plot
  autocorr.plot(chain, auto.layout = FALSE)
  ess <- effectiveSize(chain)
  autocorrelation_time <- (10000-burn)/ess
  title(varname)
  
  # Geweke z-score
  tryCatch(geweke.plot(chain, auto.layout = FALSE), error=function(e) plot(1,1))
  
  cat(varname, ess, autocorrelation_time, "\n")
}


# Report all autocorrelations in a single plot
par(mfrow=c(1,1))
df.temp <- df.traces[,c('alpha', 'mean_mu_a0', 
                        'det.R_a0.', 'det.W_a0.', 'beta_a', 
                        'det.R_b0.', 'det.W_b0.', 'beta_b')]

# Prettier without mean_mu_a0
df.temp <- df.traces[,c('alpha', 
                        'det.R_a0.', 'det.W_a0.', 'beta_a', 
                        'det.R_b0.', 'det.W_b0.', 'beta_b')]
df.ac <- data.frame(matrix(NA, nrow=1001, ncol=ncol(df.temp)))
names(df.ac) <- names(df.temp)
for (i in 1:(ncol(df.temp))){
  df.ac[,i] <- c(acf(df.temp[,i], plot=FALSE, lag.max =1000)$acf)
}

df.ac$lag <- 1:nrow(df.ac) # insert lag variable
df.ac <- melt(df.ac, id='lag')

p1 <- ggplot(df.ac, aes(x=lag, y=value, color=variable, linetype=variable))
p1 <- p1 + geom_line(aes(y=value)) # add lines
p1 <- p1 + geom_point(data=df.ac[seq(1,nrow(df.ac), by=50),], aes(shape=variable), size = 3) # add shape points

p1 <- p1 + theme(panel.background = element_blank(), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.title =  element_blank(),
                 aspect.ratio = 4/10,
                 legend.position=c(0.9, 0.7))
p1 <- p1 + ylab("Autocorrelation")
print(p1)


# Plot z traces in a different way
par(mfrow=c(1,1))
image(as.matrix(traces[['z']]))

# Plot histogram of number of clusters
####################################+
par(mfrow=c(1,2))
nclusters <- apply(traces[['z']][-(1:burn),], 1, function(x) {length(unique(x))})

# plot histogram of number of clusters
hist(nclusters, breaks=100, col="grey", main=NA, xlab="Number of clusters", xlim=c(1,max(nclusters)),  cex.lab=1.2)
box()

# plot size of clusters
sizes <- apply(traces[['z']][-(1:burn),], 1, function(x) {as.vector(sort(table(x)))})
z.matrix <- as.matrix(traces[['z']][-(1:burn),])
sizes <- apply(z.matrix, 1, function(x) {as.vector((sort(table(x), decreasing = TRUE)))})
sizes.all <-  matrix(0, nrow=length(sizes), ncol=max(nclusters))
for (i in 1:nrow(sizes.all)){
  sizes_i <- as.vector((unlist(sizes[i])))
  sizes.all[i,1:length(sizes_i)] <- sizes_i
}
sizes <- data.frame(sizes.all)
means <- colMeans(sizes)
devs <- apply(sizes, 2, sd)^2
plot(means, pch=19, xlab="Clusters", ylab="Number of users", ylim=c(0,21), cex.lab=1.2)
arrows(1:length(means), means-devs/2, 1:length(means), means+devs/2, length=0.075, code=3, angle=90)

dev.copy(png, "ncomponents.png", height=215, unit='mm', res=72)
dev.off()

# plot traces of number of clusters
plot(1:length(nclusters), nclusters, pch=3, type='l', xlab='Iteration', ylab="Number of clusters")