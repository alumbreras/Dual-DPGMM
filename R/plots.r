library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)

sample_predictive_posterior <- function(P_test, traces_coefficients, traces_intercept, traces_noise_inv){
  nsamples <- nrow(traces_coefficients)
  
  # add intercept to the rets of coefficients
  P_test <- rbind(1,P_test)
  traces_coefficients <- cbind(traces_intercept, traces_coefficients)
  T_test = dim(P_test)[2]
  y_preds <- matrix(0, nrow=nsamples, T_test)
  for(i in 1:nsamples){
    # a prediction over all the test threads
    b = traces_coefficients[i,]
    noise = 1/sqrt(traces_noise_inv[i])
    means <- t(P_test)%*%b
    y_preds[i,] <- sapply(1:T_test, function(i) {rnorm(1, mean=means[i], sd=noise)})
  }
  # return samples from the posterior (one sample is a prediction over all the threads in the test set)
  y_preds
}


plot_predictive_posterior <- function(samples){
  # """
  # Plot posterior mean with 50 and 95 confidence bands.
  #
  # Parameters:
  # ==========
  # samples: array with one sample of posterior in every row
  x <- 1:dim(samples)[2]
  
  # Posterior mean and credible intervals
  df <- data.frame(thread=1:ncol(samples), mean=colMeans(samples))
  df$ci_upper50 <- apply(samples, 2, quantile, 0.75)
  df$ci_lower50 <- apply(samples, 2, quantile, 0.25)
  df$ci_upper95 <- apply(samples, 2, quantile, 0.975)
  df$ci_lower95 <- apply(samples, 2, quantile, 0.025)
  
  p <- ggplot(df, aes(x=thread, y=mean))
  p <- p + geom_ribbon(aes(ymin=ci_lower95, ymax=ci_upper95), fill=gray(0.5), alpha=1)
  p <- p + geom_ribbon(aes(ymin=ci_lower50, ymax=ci_upper50), fill=gray(0.8), alpha=1)
  p <- p + geom_line(aes(x=thread, y=mean), color="black")
  p <- p + xlab('Threads')
  p <- p + ylab('Predictive posterior')
  p <- p + labs(color=NULL)
  p <- p + scale_colour_discrete(guide = FALSE)
  p <- p + scale_alpha_identity(guide = FALSE)
  p <- p + theme(text = element_text(size = 15),
                 panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.key = element_rect(fill = "white"),
                 legend.position=c(0.1, 0.8))
  print(p)
  p
}
plot_clustering_posterior <-function(samples){
  # Plot posterior of pairwise co-ocurrence matrix
  # Args:
  # samples: traces of cluster assignments
  coincidences <- sapply(1:ncol(samples), function(i){ colSums(samples[,i]==samples) })
  
  # From X1, X2 to 1,2...
  rownames(coincidences) <- 1:nrow(coincidences)
  p <- ggplot(melt(coincidences), aes(Var1,Var2, fill=value)) + geom_raster(vjust=0.5, hjust=0.5)
  p <- p + scale_fill_continuous(guide = FALSE)
  p <- p + coord_fixed() # preserve aspect ratio
  p <- p + xlab('users')
  p <- p + ylab('users')
  p <- p + scale_y_continuous(trans = "reverse", expand=c(0,0), breaks=c(10,20,3,40,50))
  p <- p + scale_x_continuous(expand=c(0,0), breaks=c(10,20,3,40,50))
  p <- p + theme(text = element_text(size = 15),
                 panel.background = element_blank(),
                 #axis.ticks = element_blank(),
                 #axis.text = element_blank(),
                 axis.line = element_blank())
  p
}
plot_posterior <- function(experiment){
  
  burnin <- 100
  
  # Get traces of the posterior distribution, and test set
  traces.path <- file.path(experiment, 'traces')
  plots.path <- file.path(experiment, 'plots') # store figure here
  P_test <- read.csv(file.path('data', dataset, "test_participations_50.csv"), sep='\t')
  y_test <- read.csv(file.path('data', dataset, "test_lengths_50.csv"), sep='\t')$y
 
  # Posterior over the parameters that allow the prediction
  traces_coefficients <- read.csv(file.path(traces.path, 'traces.coefficients.trc'), sep='')[-c(1:burnin),]
  traces_intercept <- read.csv(file.path(traces.path, 'traces.intercept.trc'), sep='')[-c(1:burnin),]$intercept
  traces_noise_inv <- read.csv(file.path(traces.path, 'traces.noise_inv.trc'), sep='')[-c(1:burnin),]
  traces_coefficients <- as.matrix(traces_coefficients)
  traces_noise_inv <- as.matrix(traces_noise_inv)
  
  # Take samples from posterior predictive distribution
  y_samples <- sample_predictive_posterior(P_test, traces_coefficients, traces_intercept, traces_noise_inv)
  
  # Plot posterior distribution (with threads sorted by length)
  idx_sorted <- order(y_test)
  p1 <- plot_predictive_posterior(y_samples[,idx_sorted])
  p1 <- p1 + geom_line(y=y_test[idx_sorted], color="red") # add true lengths
  
  # Plot confusion matrix
  traces_z <- read.csv(file.path(traces.path, 'traces.z.trc'), sep='')[-c(1:burnin),]
  p2 <- plot_clustering_posterior(traces_z)
  
  # Put both plots together
  g1<-ggplotGrob(p1)
  g2<-ggplotGrob(p2)
  grid.newpage()
  grid.draw(cbind(g1, g2, size = "first"))
  
  #..and save to file
  dir.create(file.path(plots.path), recursive=TRUE)
  dev.copy(png, file.path(plots.path, 'posterior.png'))
  dev.off()
  g <- grid.grab()
  g<- grid.arrange(g1, g2, widths=c(0.5,0.5))
  ggsave(file.path(plots.path, 'posterior_ggsave.png'), plot = g)
}