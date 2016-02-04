#
# Fast way of getting a matrix with the total cooccurrences from a MCMC trace
################################################################################"

cluster.assignments <- t(replicate(100,sample(10, replace=TRUE)))

n_samples <- dim(cluster.assignments)[1]
n_dim <- dim(cluster.assignments)[2]
pairwise <- matrix(0, n_dim, n_dim)
for(i in 1:n_samples){
    pairwise <- pairwise + 
        apply(t(cluster.assignments[i,]), 2, function(x) as.numeric(x==cluster.assignments[i,]))
}


M <- t(replicate(50,sample(50, replace=TRUE)))
M <- pairwise
image(1:dim(M)[1], 1:dim(M)[1], M, col= gray((0:32)/32), xlab="", ylab= "", asp=1, axes=FALSE, xaxt="n", yaxt="n")
axis(1,seq(0,50, by=10),seq(0,50, by=10), pos=0)
axis(2,seq(0,50, by=10),seq(0,50, by=10), pos=0)