# Name: 02_import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 16/09/2017
# Desc: imports and merge the patient and healthy data files

csFiles = list.files('dataExternal/healthyData/', pattern = 'N\\d+.csv$', full.names = T)
l = lapply(csFiles, function(x) read.csv(x, header=T, sep=',', na.strings = 'na'))
dfData = do.call(rbind, l)
dim(dfData)

dfPatient = read.csv('dataExternal/healthyData/merged.files_and_annotations.csv',
                     header=T, sep='\t', na.strings = c('na', 'NA'))

table(colnames(dfPatient) %in% colnames(dfData))
f = colnames(dfPatient)[colnames(dfPatient) %in% colnames(dfData)]
dfPatient = dfPatient[,f]

## check structure of both data frames
str(dfData)
str(dfPatient)

## some levels are coded differently in patient table
levels(dfData$Transcription.factor)
levels(dfPatient$Transcription.factor)
# coded differently correct
i = sapply(dfPatient$Transcription.factor, function(x) gsub("NF-kB", 
                                                            replacement =  "NF-κB", x))
dfPatient$Transcription.factor = factor(i)                                         
identical(levels(dfData$Transcription.factor), levels(dfPatient$Transcription.factor))

## check next factor
levels(dfData$Stimulation)
levels(dfPatient$Stimulation)

i = sapply(dfPatient$Stimulation, function(x) gsub("IFN-gamma", "IFNα", x))
i = sapply(i, function(x) gsub("TNF-alpha", "TNFα", x))
i = factor(i)
identical(levels(dfData$Stimulation), levels(i))
dfPatient$Stimulation = i

## check next factor
identical(levels(dfData$Cell.type), levels(dfPatient$Cell.type))

## merge the 2 tables
dfData = rbind(dfData, dfPatient)
dfData = droplevels.data.frame(dfData)
str(dfData)

dfData$Rd.score = as.numeric(dfData$Rd.score)
f = is.na(dfData$Rd.score)

dfData = dfData[!f,]

dfData = droplevels.data.frame(dfData)

####### data distribution
pairs(dfData[,-c(1, 7, 8, 11)], pch=20)

library(lattice)
library(MASS)
library(car)

# xyplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p', 'r'), 
#        index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
#        par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Rd.score ~ Cell.type | Stimulation, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Rd.score ~ Transcription.factor | Stimulation+Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Treatment,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

stripplot(Rd.score ~ Transcription.factor | Stimulation+Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Treatment,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))



## remove the highest rd score as it may be just an error
f = which(dfData$Rd.score > 4)
dfData[f,]
## cell count also very low here
dfData = dfData[-f,]

xyplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Rd.score ~ Cell.type | Stimulation, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

## cell count variable check this 
x = dfData$Cell.count
hist(x)
quantile(x)
cut.pts = cut(x, breaks = quantile(x, probs = 0:10/10), include.lowest = T, labels = 0:9)
densityplot(~ dfData$Rd.score | cut.pts)
dfData$Cell.count.grps = cut.pts

f = as.character(cut.pts)
f[f %in% c('0')] = 'g1'
f[f != 'g1'] = 'g2'

densityplot(~ dfData$Rd.score | f)

## choose a minimum cutoff of 30ish cells
dfData = dfData[dfData$Cell.count.grps != '0',]

xyplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Transcription.factor,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

xyplot(Rd.score ~ Transcription.factor | Stimulation+Cell.type, data=dfData, type=c('g', 'p'), pch=19, cex=0.5,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))


xyplot(Rd.score ~ Cell.type | Stimulation, data=dfData, type=c('g', 'p'), pch=19, groups=Transcription.factor,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

# define the modules as used by previous analyst
# module is cell type + stimulation
xyplot(Rd.score ~ Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Stimulation,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=5))

stripplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='90 Modules')

barchart(Rd.score ~ Stimulation| Cell.type, data=dfData, type=c('g', 'p'), pch=19,
          #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
          par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='90 Modules')

xyplot(Rd.score ~ Transcription.factor | Stimulation+Cell.type, data=dfData, type=c('g', 'p'), pch=19, cex=0.3,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))


## drop modules with RD score < 0.3
fModule = factor(dfData$Cell.type:dfData$Stimulation)
nlevels(fModule)
i = tapply(dfData$Rd.score, fModule, mean)
summary(i)
i = which(i < 0.3)
i = levels(fModule)[i]
# drop these modules from the dataset
f = which(fModule %in% i)
dfData = dfData[-f,]
dim(dfData)
dfData = droplevels.data.frame(dfData)

# make modules again
fModule = factor(dfData$Cell.type:dfData$Stimulation)
nlevels(fModule)

xyplot(Rd.score ~ Transcription.factor | Stimulation+Cell.type, data=dfData, type=c('g', 'p'), pch=19, cex=0.3,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

xyplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Transcription.factor, cex=0.5,
          #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
          par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules',
       auto.key = list(columns=3))

dfData$fModule = fModule

xyplot(Rd.score ~ Treatment | fModule, data=dfData, type=c('g', 'p'), pch=19, groups=Transcription.factor, cex=0.5,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules',
       auto.key = list(columns=3))

xyplot(Rd.score ~ Treatment | fModule, data=dfData, type=c('g', 'p'), pch=19, groups=Visit..Week., cex=0.5,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules',
       auto.key = list(columns=3))

xyplot(Rd.score ~ Treatment | fModule, data=dfData[dfData$Visit..Week. == 'Baseline',]
       , type=c('g', 'p'), pch=19, groups=Transcription.factor, cex=0.5,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules - baseline',
       auto.key = list(columns=3))


# barchart(Rd.score ~ Stimulation| Cell.type, data=dfData, type=c('g', 'p'), pch=19,
#          par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules')

# check data distribution
x = na.omit(dfData$Rd.score)
hist(x)
fitdistr(x, densfun = 'normal')
qqPlot(x, distribution = 'norm')

########################################################
library(LearnBayes)
set.seed(123) # for replication

ivResp = dfData$Rd.score#[dfData$Visit..Week. == 'Baseline']
ivResp = ivResp+abs(min(ivResp))+1
ivResp = log(ivResp)
summary(ivResp)
sd(ivResp)

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
muSample.op = mSir[,'mu'] #rnorm(1000, mean(lData$vector), exp(sigSample.op)/sqrt(length(lData$vector)))
# # second way of taking the sample
# muSample2.op = rnorm(1000, fit$mode['mu'], se['mu'])

fit$mode['mu']-1.96*se['mu']; fit$mode['mu']+1.96*se['mu']
quantile(muSample.op, c(0.025, 0.975))

### if we look at the histogram of the measurements
hist(ivResp, xlab='RD Score', main='', breaks=10)

########### Model checking
## sample X values, 20 times, each time drawing a fresh draw of sd and mean from the joint posterior
mDraws = matrix(NA, nrow = length(ivResp), ncol=20)

for (i in 1:20){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  mDraws[,i] = rnorm(length(ivResp), m, s)
}

p.old = par(mfrow=c(2, 2))
garbage = apply(mDraws, 2, function(x) hist(x, main='', xlab='', ylab=''))
hist(ivResp, xlab='rd score', main='')

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
colnames(mChecks) = c('Normal', 't')
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
hist(t1, xlab='Test Quantity - Variance (Normal Model)', main='', breaks=10)
abline(v = var(lData$vector), lwd=2)
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

## normal model seems appropriate here
yresp = density(ivResp)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2, ylim=c(0, 1.1))
temp = apply(mDraws, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)

plot(density(ivResp))

########################################3 try second distrubution
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
  log.prior = dcauchy(sigma, 0, 2.5, log=T) + dcauchy(nu, 0, 1, log=T)
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

# sigSample.op = rnorm(1000, fit3$mode['sigma'], se3['sigma'])
# nuSample = exp(rnorm(1000, fit3$mode['nu'], se3['nu']))
sigSample.op = mSir[,'sigma']
nuSample = exp(mSir[,'nu'])
# threshold the sample values to above or equal to 1
nuSample[nuSample < 1] = 1

## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
# muSample.op = rnorm(1000, mean(lData$vector), exp(sigSample.op)/sqrt(length(lData$vector)))
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
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2, ylim=c(0, 2.4))
temp = apply(mDraws.t, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
lines(yresp, lwd=2)














### save the data for use
dfData$fModule = fModule
write.csv(dfData, file='dataExternal/healthyData/importedHealthy.csv', row.names = F)