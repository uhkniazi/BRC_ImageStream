data {
    int<lower=1> Ntotal; // number of observations
    real y[Ntotal]; // observed data
}
parameters {
    real mu; // posterior mean
    real<lower=0> sig; // posterior scale
    real<lower=0.1> nu; // normality parameter for t distribution or degree of freedom 
}
model {
  nu ~ exponential(1/29.0);
  sig ~ cauchy(0, 2.5);
  y ~ student_t(nu, mu, sig);
}

