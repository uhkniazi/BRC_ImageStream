# Name: 03_rdScoreHealthyVsPatientBaseline.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 18/01/2019
# Desc: Import and select Baseline data to compare healthy and patients

dfData = read.csv('dataExternal/healthyData/mergedHealthyPatientData.csv',
                     header=T, na.strings = c('na', 'NA'))

dim(dfData)

## use only baseline samples and merge 2 types of treated into one group
dfData = dfData[dfData$Visit..Week. == 'Baseline',]
dfData = droplevels.data.frame(dfData)
t = as.character(dfData$Treatment)
t[t != 'None'] = 'Patient'
t[t == 'None'] = 'Healthy'
dfData$Treatment = factor(t)

## choose only the transcription factor nfkb
levels(dfData$Transcription.factor)
dfData = dfData[dfData$Transcription.factor == 'NF-κB', ]
dfData = droplevels.data.frame(dfData)
dim(dfData)

## cell count variable check this 
dim(dfData)
dfData = dfData[dfData$Cell.count >= 10,]
dim(dfData)
dfData = droplevels.data.frame(dfData)
str(dfData)

## modules of average rd score of less than 0.3 have already been removed previously
## as we are using the data version where this filter has already been applied.

dfData$Modules = factor(dfData$fModule:dfData$Treatment)
nlevels(dfData$Modules); table(dfData$Modules)
dfData = dfData[order(dfData$Patient.ID, dfData$Modules),]
dfData = droplevels.data.frame(dfData)
str(dfData)

## make a plot of the raw data
library(lattice)
## each module has only one transcription factor, and we are only looking at nfkb this time
xyplot(Rd.score ~ Treatment | fModule, data=dfData, type=c('g', 'p'), pch=19, cex=0.6,
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))


library(lme4)
fit.lme1 = lmer(Rd.score ~ 1 + (1 | Modules) + (1 | Patient.ID), data=dfData, REML=F)
summary(fit.lme1)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='tResponse2RandomEffectsNoFixed.stan')

gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

l = gammaShRaFromModeSD(sd(dfData$Rd.score), 2*sd(dfData$Rd.score))
# m = model.matrix(Rd.score ~ Transcription.factor, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Patient.ID), Nclusters2=nlevels(dfData$Modules), 
                 NgroupMap1=as.numeric(dfData$Patient.ID), NgroupMap2=as.numeric(dfData$Modules), 
                 Ncol=1, #X=m,
                 y=dfData$Rd.score, gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=4000, chains=4, pars=c('betas', 'sigmaRan1', 'sigmaRan2',
                                                                         'sigmaPop','nu', 'mu', 'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=4)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)

###### some model checks
### check model fit
m = extract(fit.stan, 'mu')
names(m)
dim(m$mu)
fitted = apply(m$mu, 2, mean)

plot(dfData$Rd.score, fitted, pch=20, cex=0.5)
plot(dfData$Rd.score, dfData$Rd.score - fitted, pch=20, cex=0.5)
iResid = scale(dfData$Rd.score - fitted)
plot(fitted, iResid, pch=20, cex=0.5)
lines(lowess(fitted, iResid), col=2, lwd=2)

### plot the posterior predictive values
m = extract(fit.stan, c('mu', 'nu', 'sigmaPop'))
i = sample(1:8000, 200)
muSample = m$mu[i,]
nuSample = m$nu[i]
sigSample = m$sigmaPop[i]

## t sampling functions
dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu

## use point-wise predictive approach to sample a new value from this data
ivResp = dfData$Rd.score
mDraws = matrix(NA, nrow = length(ivResp), ncol=200)

for (i in 1:ncol(mDraws)){
  mDraws[,i] = rt_ls(length(ivResp), nuSample[i], muSample[i,], sigSample[i])
}

yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)#, ylim=c(0, 1))
temp = apply(mDraws, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})

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
# ## is the model adequate except for the extreme tails
# T1_symmetry = function(Y, th){
#   Yq = quantile(Y, c(0.90, 0.10))
#   return(abs(Yq[1]-th) - abs(Yq[2]-th))
# } 

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
colnames(mChecks) = c('t', 'gamma')

t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 1] = getPValue(t1, var(ivResp))

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(ivResp)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(ivResp)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(ivResp)
mChecks['Mean', 1] = getPValue(t1, t2)

mChecks

## get the coefficient of interest - Modules in our case from the random coefficients section
mModules = extract(fit.stan)$rGroupsJitter2
dim(mModules)
## get the intercept at population level
iIntercept = extract(fit.stan)$betas[,1]
## add the intercept to each random effect variable, to get the full coefficient
mModules = sweep(mModules, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mModules), mods=levels(dfData$Modules), rawAverage=tapply(dfData$Rd.score, 
                                                                                     dfData$Modules, FUN = mean))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:3], c('cells', 'stimulation', 'treatment'))
d$split = factor(d$cells:d$stimulation)

## this data frame is a mapper for each required comparison
ldfMap = split(d, f = d$split)

## get a p-value for each comparison
l = lapply(ldfMap, function(x) {
  c = x$cols
  d = getDifference(ivData = mModules[,c[2]], ivBaseline = mModules[,c[1]])
  r = data.frame(cells= as.character(x$cells[1]), stimulation=as.character(x$stimulation[1]), coef.healthy=mean(mModules[,c[1]]), 
        coef.patient=mean(mModules[,c[2]]), zscore=d$z, pvalue=d$p, rawAverage.healthy=x$rawAverage[1], 
        rawAverage.patient=x$rawAverage[2])
  return(format(r, digi=3))
})

dfResults = do.call(rbind, l)
dfResults$p.adj = format(p.adjust(dfResults$pvalue, method='bonf'), digi=3)
write.csv(dfResults, file='Results/rdScoreHealthyVSPatientsBaseline.csv', row.names = F)

## make the plots for the raw data and fitted data
## format data for plotting
d2 = dfData[,c('Rd.score', 'Modules')]
f = strsplit(as.character(d2$Modules), ':')
d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('cells', 'stimulation', 'treatment'))

dotplot(treatment ~ Rd.score | cells:stimulation, data=d2, groups=treatment, panel=function(x, y, ...) panel.bwplot(x, y, pch='|', ...),
        par.strip.text=list(cex=0.6), main='Raw data 21 Modules at baseline', xlab='RD Score')

## format data for plotting
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Modules)
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('cells', 'stimulation', 'treatment'))

dotplot(treatment ~ m+s1+s2 | cells:stimulation, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=0.6), main='Regression Coeff 21 Modules at baseline', xlab='Model estimated Average RD Score')


