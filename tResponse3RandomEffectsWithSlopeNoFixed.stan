data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each observation to group 2 
  int<lower=1> Ncol; // total number of columns in model matrix
  vector[Ntotal] X; // vector of covariate for slope
  real y[Ntotal]; // response variable normally distributed
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  real intercept; // population intercept
  real slope; // population slope
  real<lower=0.01> sigmaRan1; // random effect standard deviation for group 1
  real<lower=0.01> sigmaRan2; // random effect standard deviation for group 2
  real<lower=0.01> sigmaRanSlope1; // random effect slope
  real<lower=0.01> sigmaPop; // population standard deviation
  vector[Nclusters1] rGroupsJitter1; // number of random intercept jitters for each level of cluster/group 1
  vector[Nclusters1] rGroupsSlope1; // number of random slope jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random intercept jitters for each level of cluster/group 2
  real<lower=1> nu; // normality parameter for t distribution or degree of freedom 
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ntotal] newSlope; // set of slopes after adding random slope jitters
  newSlope = slope + rGroupsSlope1[NgroupMap1];
  // fitted value with slopes
  for (i in 1:Ntotal){
    mu[i] = X[i] * newSlope[i];
  }
  mu = mu + intercept + rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2];
}
model {
  nu ~ exponential(1/29.0);
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  sigmaRan2 ~ gamma(gammaShape, gammaRate);
  sigmaRanSlope1 ~ gamma(gammaShape, gammaRate);
  sigmaPop ~ gamma(gammaShape, gammaRate);
  intercept ~ cauchy(0, 10);
  slope ~ cauchy(0, 2.5);
  // random effects sample
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  rGroupsSlope1 ~ normal(0, sigmaRanSlope1);
  rGroupsJitter2 ~ normal(0, sigmaRan2);
  // likelihood function
  y ~ student_t(nu, mu, sigmaPop);
}
