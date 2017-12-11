# Name: 10_importMisLongitudinalAdalimumab.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 08/12/2017
# Desc: imports the longitudinal data for treatment 1


dfPatient = read.csv('dataExternal/healthyData/merged.files_and_annotations.csv',
                     header=T, sep='\t', na.strings = c('na', 'NA'))

dim(dfPatient)
colnames(dfPatient)
cn = c("Patient.ID", "Treatment", "Visit..Week.", "Transcription.factor", "Median.internalization.score",
       "Stimulation", "Cell.type" , "Rd.score", "Cell.count", "Baseline.PASI", "Week.1.PASI", "Week.4.PASI", "Week.12.PASI",  "Gender", "Age", "Ethnicity", "PsA", "Weight.baseline..Kg.", "ADL.conc", "USK.conc", "relativePASI", "PASI90", "PASI75", "PASI50")
length(cn)
dfData = dfPatient[,cn]

## check structure of data
str(dfData)
dfData = dfData[dfData$Treatment == 'Adalimumab',]
dfData = droplevels.data.frame(dfData)
## check stimulation factor
levels(dfData$Stimulation)

f = is.na(dfData$Median.internalization.score)
table(f)
dfData = dfData[!f,]

## this should clear any empty factor levels
dfData = droplevels.data.frame(dfData)

## make a new time variable
levels(dfData$Visit..Week.)
dfData$Visit..Week. = factor(dfData$Visit..Week., levels=c('Baseline', 
                                                           'Week 1', 'Week 4',
                                                           'Week 12'))
i = rep(NA, length=length(dfData$Visit..Week.))
i[dfData$Visit..Week. == 'Baseline'] = 0;
i[dfData$Visit..Week. == 'Week 1'] = 1;
i[dfData$Visit..Week. == 'Week 4'] = 4;
i[dfData$Visit..Week. == 'Week 12'] = 12;
dfData$time = i

####### data distribution
library(lattice)
library(MASS)
library(car)

xyplot(Median.internalization.score ~ time | Cell.type:Transcription.factor:Stimulation, data=dfData, type=c('g', 'p', 'r'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Median.internalization.score ~ time | Cell.type:Transcription.factor:Stimulation, data=dfData, type=c('smooth'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

## population level effect
xyplot(Median.internalization.score ~ time, data=dfData, type=c('smooth'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

## time as categorical variable
xyplot(Median.internalization.score ~ Visit..Week. | Cell.type:Transcription.factor:Stimulation, data=dfData, type=c('r', 'g', 'p'), pch=19, cex=0.6,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

####### apply the cutoffs
## cell count variable check this 
dim(dfData)
dfData = dfData[dfData$Cell.count >= 10,]
dim(dfData)

xyplot(Median.internalization.score ~ time | Cell.type:Transcription.factor:Stimulation, data=dfData, type=c('g', 'p', 'r'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), groups=dfData$Treatment)

# define the modules as used by previous analyst
# in each cell type, the stimulation and transcription factors are nested 
# apart from where Stimulation = Unstimulated
nlevels(factor(dfData$Stimulation:dfData$Transcription.factor:dfData$Cell.type))
nlevels(factor(dfData$Stimulation:dfData$Cell.type))
nlevels(factor(dfData$Stimulation:dfData$Transcription.factor))
nlevels(dfData$Stimulation)
xtabs(~ dfData$Stimulation+dfData$Transcription.factor+dfData$Cell.type)
# module is cell type + stimulation + transcription factor
## drop modules with average RD score < 0.3
fModule = factor(dfData$Cell.type:dfData$Stimulation:dfData$Transcription.factor)
nlevels(fModule)
## 120 levels
# i = tapply(dfData$Median.internalization.score, fModule, mean)
# summary(i)
# i = which(i < 0.3)
# i = levels(fModule)[i]
# # drop these modules from the dataset
# f = which(fModule %in% i)
# dim(dfData)
# dfData = dfData[-f,]
# dim(dfData)
# dfData = droplevels.data.frame(dfData)

# make modules again
fModule = factor(dfData$Cell.type:dfData$Stimulation:dfData$Transcription.factor)
nlevels(fModule)

dfData$fModule = fModule

xtabs(~ dfData$fModule  + dfData$Transcription.factor)
rm(fModule)
## each module has only one transcription factor
xyplot(Median.internalization.score ~ time | fModule, data=dfData, type=c('g', 'p', 'r'), pch=19, cex=0.6,
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Median.internalization.score ~ time | fModule, data=dfData, type=c('smooth'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

## population level effect
xyplot(Median.internalization.score ~ time, data=dfData, type=c('smooth'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

## time as categorical variable
xyplot(Median.internalization.score ~ Visit..Week. | fModule, data=dfData, type=c('r', 'g', 'p'), pch=19, cex=0.6,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Median.internalization.score ~ Visit..Week. | fModule, data=dfData, type=c('smooth'), pch=19, cex=0.6,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))


########################################################
library(LearnBayes)
set.seed(123) # for replication

ivResp = dfData$Median.internalization.score
# ivResp = ivResp+abs(min(ivResp))+2
# ivResp = log(ivResp)
summary(ivResp)
sd(ivResp)

length(ivResp)

plot(density(ivResp))


############# try a normal distribution
## define a log posterior function
lp = function(theta, data){
  # we define the sigma on a log scale as optimizers work better
  # if scale parameters are well behaved
  s = exp(theta[2])
  m = theta[1]
  d = data$vector # observed data vector
  log.lik = sum(dnorm(d, m, s, log=T))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=mean(ivResp), 'sigma'=log(sd(ivResp)))
lData = list('vector'=ivResp)
lp(start, lData)

op = optim(start, lp, control = list(fnscale = -1), data=lData)
op$par
exp(op$par[2])

## try the laplace function from LearnBayes
fit = laplace(lp, start, lData)
fit
se = sqrt(diag(fit$var))
se
fit$mode+1.96*se
fit$mode-1.96*se

### lets use the results from laplace and SIR sampling with a t proposal density
### to take a sample
tpar = list(m=fit$mode, var=fit$var*2, df=3)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(lp, tpar, 1000, lData)
apply(s, 2, mean)[1]
exp(apply(s, 2, mean)[2])
mSir = s

## lets take a sample from this
sigSample.op = mSir[,'sigma']
muSample.op = mSir[,'mu']
fit$mode['mu']-1.96*se['mu']; fit$mode['mu']+1.96*se['mu']
quantile(muSample.op, c(0.025, 0.975))

########### Model checking
## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
## variance
T1_var = function(Y) return(var(Y))
## is the model adequate except for the extreme tails
T1_symmetry = function(Y, th){
  Yq = quantile(Y, c(0.90, 0.10))
  return(abs(Yq[1]-th) - abs(Yq[2]-th))
} 

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## mChecks
mChecks = matrix(NA, nrow=5, ncol=2)
rownames(mChecks) = c('Variance', 'Symmetry', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('Normal', 'mixture')
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(ivResp), ncol=200)
mThetas = matrix(NA, nrow=200, ncol=2)
colnames(mThetas) = c('mu', 'sd')

for (i in 1:200){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  mDraws[,i] = rnorm(length(ivResp), m, s)
  mThetas[i,] = c(m, s)
}

mDraws.norm = mDraws
### get the test quantity from the test function
t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 1] = getPValue(t1, var(lData$vector))

## test for symmetry
t1 = sapply(seq_along(1:200), function(x) T1_symmetry(mDraws[,x], mThetas[x,'mu']))
t2 = sapply(seq_along(1:200), function(x) T1_symmetry(lData$vector, mThetas[x,'mu']))
plot(t2, t1, pch=20, xlab='Realized Value T(Yobs, Theta)',
     ylab='Test Value T(Yrep, Theta)', main='Symmetry Check (Normal Model)')
abline(0,1)
mChecks['Symmetry', 1] = getPValue(t1, t2) 
## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lData$vector)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lData$vector)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lData$vector)
mChecks['Mean', 1] = getPValue(t1, t2)

mChecks

## normal model seems ok here apart from borderline for outlier and symmetry
yresp = density(ivResp)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2, ylim=c(0, 1.1))
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
hist(ivResp, prob=T)

temp = apply(mDraws, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})

#######################################
############## try a t-distribution
#######################################
lp3 = function(theta, data){
  # function to use to use scale parameter
  ## see here https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
  dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
  ## likelihood function
  lf = function(dat, pred){
    return(log(dt_ls(dat, nu, pred, sigma)))
  }
  nu = exp(theta['nu']) ## normality parameter for t distribution
  sigma = exp(theta['sigma']) # scale parameter for t distribution
  m = theta[1]
  d = data$vector # observed data vector
  if (nu < 1) return(-Inf)
  log.lik = sum(lf(d, m))
  log.prior = dcauchy(sigma, 0, 2.5, log=T) + dexp(nu, 1/29, log=T)
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=mean(ivResp), 'sigma'=log(sd(ivResp)), 'nu'=log(2))
lp3(start, lData)

op = optim(start, lp3, control = list(fnscale = -1), data=lData)
op$par
exp(op$par[2:3])

## try the laplace function from LearnBayes
fit3 = laplace(lp3, start, lData)
fit3
se3 = sqrt(diag(fit3$var))

### lets use the results from laplace and SIR sampling with a t proposal density
### to take a sample
tpar = list(m=fit3$mode, var=fit3$var*2, df=3)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(lp3, tpar, 1000, lData)
apply(s, 2, mean)[1]
exp(apply(s, 2, mean)[c(2,3)])
mSir = s

sigSample.op = mSir[,'sigma']
nuSample = exp(mSir[,'nu'])
# threshold the sample values to above or equal to 1
nuSample[nuSample < 1] = 1

## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
muSample.op = mSir[,'mu']

########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(ivResp), ncol=200)
mThetas = matrix(NA, nrow=200, ncol=3)
colnames(mThetas) = c('mu', 'sd', 'nu')

for (i in 1:200){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  n = nuSample[p]
  mDraws[,i] = rt_ls(length(ivResp), n, m, s)
  mThetas[i,] = c(m, s, n)
}

mDraws.t = mDraws
## get the p-values for the test statistics
t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 2] = getPValue(t1, var(lData$vector))

## test for symmetry
t1 = sapply(seq_along(1:200), function(x) T1_symmetry(mDraws[,x], mThetas[x,'mu']))
t2 = sapply(seq_along(1:200), function(x) T1_symmetry(lData$vector, mThetas[x,'mu']))
plot(t2, t1, pch=20, xlab='Realized Value T(Yobs, Theta)',
     ylab='Test Value T(Yrep, Theta)', main='Symmetry Check (T Distribution)')
abline(0,1)
mChecks['Symmetry', 2] = getPValue(t1, t2) 

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lData$vector)
mChecks['Min', 2] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lData$vector)
mChecks['Max', 2] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lData$vector)
mChecks['Mean', 2] = getPValue(t1, t2)

mChecks

yresp = density(ivResp)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)
temp = apply(mDraws.t, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
lines(yresp, lwd=2)

hist(ivResp, prob=T)
## t samples
temp = apply(mDraws.t, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
## normal samples
temp = apply(mDraws.norm, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})

############## generate an MCMC sample using stan
library(rstan)
stanDso = rstan::stan_model(file='fitTparam.stan')

lStanData = list(Ntotal=length(ivResp), y=ivResp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('mu', 'sig', 'nu'))
print(fit.stan)
m = extract(fit.stan)
muSample.op = m$mu
sigSample.op = m$sig
nuSample = m$nu

mDraws = matrix(NA, nrow = length(ivResp), ncol=200)
mThetas = matrix(NA, nrow=200, ncol=3)
colnames(mThetas) = c('mu', 'sd', 'nu')

for (i in 1:200){
  p = sample(1:5000, size = 1)
  s = sigSample.op[p]
  m = muSample.op[p]
  n = nuSample[p]
  mDraws[,i] = rt_ls(length(ivResp), n, m, s)
  mThetas[i,] = c(m, s, n)
}

mDraws.t = mDraws

yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2, ylim=c(0, 2.4))
temp = apply(mDraws.t, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
lines(yresp, lwd=2)

#####################################################################################
############################### fit a model using stan to estimate mixture parameters
####################################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='fitNormMixture.stan')

## stan data
lStanData = list(Ntotal=length(ivResp), y=ivResp, iMixtures=2)

# ## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(0, 1), iMixWeights=c(0.5, 0.5))
}

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2, cores=2, init=initf)
print(fit.stan, digi=3)
traceplot(fit.stan)

## check if labelling degeneracy has occured
## see here: http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
params1 = as.data.frame(extract(fit.stan, permuted=FALSE)[,1,])
params2 = as.data.frame(extract(fit.stan, permuted=FALSE)[,2,])

## check if the means from different chains overlap
## Labeling Degeneracy by Enforcing an Ordering
par(mfrow=c(1,2))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
plot(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)

par(mfrow=c(1,1))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
points(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)

############# extract the mcmc sample values from stan
mStan = do.call(cbind, extract(fit.stan))
mStan = mStan[,-(ncol(mStan))]
colnames(mStan) = c('mu1', 'mu2', 'sigma1', 'sigma2', 'mix1', 'mix2')
dim(mStan)
## get a sample for this distribution
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(ivResp), ncol=200)

for (i in 1:200){
  p = sample(1:nrow(mStan), size = 1)
  mix = mStan[p,'mix1']
  ## this will take a sample from a normal mixture distribution
  sam = function() {
    ind = rbinom(1, 1, prob = mix)
    return(ind * rnorm(1, mStan[p, 'mu1'], mStan[p, 'sigma1']) + 
             (1-ind) * rnorm(1, mStan[p, 'mu2'], mStan[p, 'sigma2']))
  }
  mDraws[,i] = replicate(length(ivResp), sam())
}

mDraws.normMix = mDraws

plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2)
temp = apply(mDraws, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)

t1 = apply(mDraws, 2, T1_var)
ivTestQuantities = getPValue(t1, var(lStanData$y))

t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lStanData$y)
ivTestQuantities = c(ivTestQuantities, getPValue(t1, t2))

t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lStanData$y)
ivTestQuantities = c(ivTestQuantities, getPValue(t1, t2))

t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lStanData$y)
ivTestQuantities = c(ivTestQuantities, getPValue(t1, t2))

names(ivTestQuantities) = c('variance', 'minimum', 'maximum', 'mean')

ivTestQuantities

hist(ivResp, prob=T)
## t samples
temp = apply(mDraws.t, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
## normal samples
temp = apply(mDraws.norm, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
## mixture samples
## normal samples
temp = apply(mDraws.normMix, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='green', lwd=0.6)
})
lines(density(ivResp))

### save the data for use
write.csv(dfData, file='dataExternal/healthyData/diseasedDataMISAdalimumab.csv', row.names = F)

