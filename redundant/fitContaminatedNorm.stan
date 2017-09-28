data {
  int<lower=1> Ntotal; // number of observations
  real y[Ntotal]; // response variable - normally distributed
}

parameters { // the parameters to track
  real mu; // means to track 
  ordered[2] sigma; // scale parameter for normal distribution  
  simplex[2] iMixWeights; // weights for the number of mixtures (should sum to one)
}
model {
  // see stan manual page 187 for an example of mixture model
  real ps[2]; // temporary variable for log components
  // any priors go here 
  sigma ~ cauchy(0, 2.5); // weak prior
  iMixWeights ~ dirichlet(rep_vector(2.0, 2));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // number of mixture components
    ps[1] = log(iMixWeights[1]) + normal_lpdf(y[n] | mu, sigma[1]);
    ps[2] = log(iMixWeights[2]) + normal_lpdf(y[n] | mu, sigma[2]);
    target += log_sum_exp(ps);
  }
}

