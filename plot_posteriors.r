if(TRUE){
  burnin <- 15000
  dataset <- 'overlapped'
  experiments <- list.dirs(file.path('out', dataset), recursive=FALSE)
  for(i in 1:length(experiments)){
    experiment <- experiments[i]
    plot_posterior(experiment)
  }
}