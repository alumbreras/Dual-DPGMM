# Simple gmm to see how much ARI do we get in iris dataset
##########################################################
# Check lower and upper bounds for my dataset
library(mclust)
data(iris)
data <- read.csv('./iris/data_users_50.csv', sep='\t')
z_true <- data$z

A <- data[,2:5]
res <- Mclust(A)
z <- res$classification
pairs(A, col=z)
ari <- adjustedRandIndex(z_true, z)
cat("ARI with A+b:", ari)
# ARI: 0.86 ! :)
# ARI 1

A <- data[,2:4]
res <- Mclust(A)
z <- res$classification
pairs(A, col=z)
ari <- adjustedRandIndex(z_true, z)
cat("ARI with A:", ari)
# ARI: 0.5 :)
# ARI: 0.68

###########################################################################
###########################################################################



# Make the fourth dimension more separated so that it helps to find 3 clusters!
# (for the paper)
#data[data$z==0,]$b <- mean(data[data$z==0,]$b)
#data[data$z==1,]$b <- mean(data[data$z==1,]$b)
#data[data$z==2,]$b <- mean(data[data$z==2,]$b)