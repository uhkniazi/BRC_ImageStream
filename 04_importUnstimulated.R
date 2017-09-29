# Name: 04_importUnstimulated.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 28/09/2017
# Desc: second part of the analysis with importing unstimulated data

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

## keep only the unstimulated data
levels(dfData$Stimulation)
dfData = dfData[dfData$Stimulation == 'Unstimulated',]
str(dfData)
## keep only the baseline
levels(dfData$Visit..Week.)
dfData = dfData[dfData$Visit..Week. == 'Baseline',]
## drop the RD score column
dfData = dfData[,-9]

f = is.na(dfData$Median.internalization.score)
table(f)
dfData = dfData[!f,]
dfData = droplevels.data.frame(dfData)
str(dfData)

####### data distribution
as.data.frame(colnames(dfData))
pairs(dfData[,-c(1, 3, 10)], pch=20)

library(lattice)
library(MASS)
library(car)

xyplot(Median.internalization.score ~ Cell.type, data=dfData, type=c('g', 'p', 'r'),
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Median.internalization.score ~ Treatment | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Median.internalization.score ~ Treatment | Cell.type+Transcription.factor, data=dfData, type=c('g', 'p'), pch=19,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

stripplot(Median.internalization.score ~ Treatment | Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Transcription.factor,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

hist(dfData$Median.internalization.score)

## cell count variable check this 
x = dfData$Cell.count
hist(x)
quantile(x)
cut.pts = cut(x, breaks = quantile(x, probs = 0:10/10), include.lowest = T, labels = 0:9)
densityplot(~ dfData$Median.internalization.score | cut.pts)
dfData$Cell.count.grps = cut.pts

f = as.character(cut.pts)
f[f %in% c('0')] = 'g1'
f[f != 'g1'] = 'g2'

densityplot(~ dfData$Median.internalization.score | f)

## choose a minimum cutoff of 30ish cells
dfData = dfData[dfData$Cell.count.grps != '0',]

xyplot(Median.internalization.score ~ Treatment | Cell.type+Transcription.factor, data=dfData, type=c('g', 'p'), pch=19,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

# define the modules as used by previous analyst
# module is cell type + stimulation
# stimulation is not valid for this data set as all cells are unstimulated
# so a module is just a cell type

xyplot(Median.internalization.score ~ Treatment | Cell.type+Transcription.factor, data=dfData, type=c('g', 'p'), pch=19, cex=0.3,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

xyplot(Median.internalization.score ~ Treatment | Cell.type, data=dfData, type=c('g', 'p'), pch=19, cex=0.3, groups=Transcription.factor,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], # layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

## do we drop modules with low average score?
fModule = factor(dfData$Cell.type)
nlevels(fModule)
i = tapply(dfData$Median.internalization.score, fModule, mean)
summary(i)

# i = which(i < 0.3)
# i = levels(fModule)[i]
# # drop these modules from the dataset
# f = which(fModule %in% i)
# dfData = dfData[-f,]
# dim(dfData)
# dfData = droplevels.data.frame(dfData)

# make modules again
fModule = factor(dfData$Cell.type)
nlevels(fModule)

dfData$fModule = fModule

########################################################
library(LearnBayes)
set.seed(123) # for replication

ivResp = dfData$Median.internalization.score
ivResp = ivResp+abs(min(ivResp))+0.1
#ivResp = log(ivResp)
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

## try the laplace function from LearnBayes
fit = laplace(lp, start, lData)
fit

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
# quick sanity check
fit$mode['mu']-1.96*se['mu']; fit$mode['mu']+1.96*se['mu']
quantile(muSample.op, c(0.025, 0.975))

### if we look at the histogram of the measurements
hist(ivResp, xlab='RD Score', main='', breaks=10)

########### Model checking
## calculate bayesian p-value for test statistic
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
colnames(mChecks) = c('Normal', 'gamma')
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

## normal model seems inappropriate here
yresp = density(ivResp)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2, ylim=c(0, 1.1))
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)

plot(density(ivResp))

########################################3 try gamma

getalphabeta.poisson = function(lambda){
  m = mean(lambda)
  v = var(lambda)
  alpha = (m^2)/v
  beta = alpha/m
  return(c(alpha=alpha, beta=beta))
}


lp2 = function(theta, data){
  alpha = exp(theta['alpha']) ## shape
  beta = exp(theta['beta']) # rate
  d = data$vector # observed data vector
  log.lik = sum(dgamma(d, shape = alpha, rate = beta, log = T))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = getalphabeta.poisson(ivResp)
start = log(start)
lData$vector = ivResp
lp2(start, lData)

op = optim(start, lp2, control = list(fnscale = -1), data=lData)
op$par
exp(op$par)

## try the laplace function from LearnBayes
fit2 = laplace(lp2, start, lData)
fit2
se2 = sqrt(diag(fit2$var))

### lets use the results from laplace and SIR sampling with a t proposal density
### to take a sample
tpar = list(m=fit2$mode, var=fit2$var*2, df=3)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(lp2, tpar, 1000, lData)
exp(apply(s, 2, mean))
mSir = s

alphaSample = exp(mSir[,'alpha'])
betaSample = exp(mSir[,'beta'])
plot(density(alphaSample))
plot(density(betaSample))

########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(ivResp), ncol=200)
mThetas = matrix(NA, nrow=200, ncol=3)
colnames(mThetas) = c('mu', 'alpha', 'beta')

for (i in 1:200){
  p = sample(1:1000, size = 1)
  a = alphaSample[p]
  b = betaSample[p]
  mDraws[,i] = rgamma(length(ivResp), shape = a, rate = b)
  mThetas[i,] = c(a/b, a, b)
}

mDraws.g = mDraws
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
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2, ylim=c(0, 1.1))
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)


### save the data for use
write.csv(dfData, file='dataExternal/healthyData/mergedDataUnstimulated.csv', row.names = F)