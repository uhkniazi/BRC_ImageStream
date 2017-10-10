data {
  int<lower=1> Ntotal; // number of observations
  real y[Ntotal]; // response variable - normally distributed
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=2> iMixtures; // number of mixture distributions
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each observation to group 2
}

parameters { // the parameters to track
  ordered[iMixtures] mu; // number of means to track Breaking the Labeling Degeneracy by Enforcing an Ordering (population Intercepts)
  real<lower=0> sigma[iMixtures]; // scale parameters for normal distribution (population sigmas) 
  simplex[iMixtures] iMixWeights; // weights for the number of mixtures (should sum to one)
  // regression coefficients and other related parameters
  real<lower=0> sigmaRan1; // random effect standard deviation for group 1
  real<lower=0> sigmaRan2; // random effect standard deviation for group 2
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
}
transformed parameters {
  vector[Ntotal] muFitted; // fitted value from linear predictor
  muFitted =  rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2];
}
model {
  // see stan manual page 187 for an example
  real ps[iMixtures]; // temporary variable for log components
  sigmaRan1 ~ cauchy(0, 2.5);
  sigmaRan2 ~ cauchy(0, 2.5);
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  rGroupsJitter2 ~ normal(0, sigmaRan2);
  // any priors for mixture components go here 
  mu[1] ~ normal(1.2, 1);
  mu[2] ~ normal(1.5, 1);
  sigma ~ normal(0, 1); // standard deviations
  iMixWeights ~ dirichlet(rep_vector(2.0, iMixtures));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // second loop for number of mixture components
    for (k in 1:iMixtures){
      ps[k] = log(iMixWeights[k]) + normal_lpdf(y[n] | mu[k] + muFitted[n], sigma[k]);
    }
    target += log_sum_exp(ps);
  }
}

