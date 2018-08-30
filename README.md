Dual Dirichlet Process - Gaussian Mixture Model
=====

Code for our paper

**Non-parametric clustering over user features and latent behavioral functions with dual-view mixture models**<br>
*Lumbreras A., Guégan M., Velcin J., Jouve B.*<br> 
Computational Statistics (2016) 
  

## Notes

  * DP-GMM: default
  * fixed-GMM: set alpha=0 and dont sample it
  * single: use sample_z() and do not sample feature view


Multi DP-GMM is a double dirichlet process that simultaneously does:

 * Cluster users with respect to their attributes.
 * Cluster users with respect to their behaviors.
 
Inference is done with Monte Carlo methods: Gibbs Sampling and Adaptive Rejection Sampling.

Read the paper for more information.


## Some interesting readings

 - [On the Bogosity of MCMC Diagnostics](http://users.stat.umn.edu/~geyer/mcmc/diag.html)
 - [Sampling from Wisharts with scalar degrees of freedom](http://stats.stackexchange.com/questions/141330/sampling-from-wishart-distributionos-with-scalar-degrees-of-freedom-upsilonp)
 - [The Wishart and Inverse Wishart distributions](http://www.tc.umn.edu/~nydic001/docs/unpubs/Wishart_Distribution.pdf)
   
