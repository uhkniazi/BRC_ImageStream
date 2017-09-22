data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each observation to group 2 
  int<lower=1> Ncol; // total number of columns in model matrix
  real y[Ntotal]; // response variable normally distributed
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real<lower=0> sigmaRan1; // random effect standard deviation for group 1
  real<lower=0> sigmaRan2; // random effect standard deviation for group 2
  real<lower=0> sigmaPop; // population standard deviation
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ntotal] rNewIntercept; // after adding the cluster deflections to population intercept i.e. beta[1]
  // add the random jitters from each group to the population intercept
  // new intercept variable equal to number of observations
  // use the mapping variable to expand the groups intercepts to beta[1] vector
  rNewIntercept = betas[1] + rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2];
  // fitted value
  // add the random intercept
  mu = rNewIntercept;
}
model {
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  sigmaRan2 ~ gamma(gammaShape, gammaRate);
  sigmaPop ~ gamma(gammaShape, gammaRate);
  betas ~ cauchy(0, 10);
  // random effects sample
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  rGroupsJitter2 ~ normal(0, sigmaRan2);
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
