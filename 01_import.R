# Name: 01_import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 14/09/2017
# Desc: imports and merges the data files

csFiles = list.files('dataExternal/healthyData/', pattern = '*.csv$', full.names = T)
l = lapply(csFiles, function(x) read.csv(x, header=T, sep=',', na.strings = 'na'))
dfData = do.call(rbind, l)
dim(dfData)

dfData$Rd.score = as.numeric(dfData$Rd.score)
f = is.na(dfData$Rd.score)

dfData = dfData[!f,]

dfData = droplevels.data.frame(dfData)

####### data distribution
pairs(dfData[,-c(1, 2, 3, 11)], pch=20)

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
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

xyplot(Rd.score ~ Transcription.factor | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))


xyplot(Rd.score ~ Cell.type | Stimulation, data=dfData, type=c('g', 'p'), pch=19, groups=Transcription.factor,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=3))

# define the modules as used by previous analyst
# module is cell type + stimulation
xyplot(Rd.score ~ Cell.type, data=dfData, type=c('g', 'p'), pch=19, groups=Stimulation,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), auto.key = list(columns=5))

stripplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='90 Modules')

barchart(Rd.score ~ Stimulation| Cell.type, data=dfData, type=c('g', 'p'), pch=19,
          index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
          par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='90 Modules')

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

stripplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
          #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
          par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules')

barchart(Rd.score ~ Stimulation| Cell.type, data=dfData, type=c('g', 'p'), pch=19,
         par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules')

# check data distribution
x = na.omit(dfData$Rd.score)
hist(x)
fitdistr(x, densfun = 'normal')
qqPlot(x, distribution = 'norm')

########################################################
library(LearnBayes)
set.seed(123) # for replication

ivResp = dfData$Rd.score
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

## lets take a sample from this
sigSample.op = rnorm(1000, fit$mode['sigma'], se['sigma'])
muSample.op = rnorm(1000, mean(lData$vector), exp(sigSample.op)/sqrt(length(lData$vector)))
# second way of taking the sample
muSample2.op = rnorm(1000, fit$mode['mu'], se['mu'])

fit$mode['mu']-1.96*se['mu']; fit$mode['mu']+1.96*se['mu']
quantile(muSample.op, c(0.025, 0.975))

### if we look at the histogram of the 66 measurements
hist(ivResp, xlab='RD Score', main='', breaks=50)

########### Model checking
## sample 66 values, 20 times, each time drawing a fresh draw of sd and mean from the joint posterior
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
mChecks = matrix(NA, nrow=5, ncol=3)
rownames(mChecks) = c('Variance', 'Symmetry', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('Normal', 'NormalCont', 'T')
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
hist(t1, xlab='Test Quantity - Variance (Normal Model)', main='', breaks=50)
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

## normal model seems appropriate here
yresp = density(ivResp)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2)
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)

plot(density(ivResp))

