# Name: 19_pasi90vsRDscoreAdalimumab.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 21/01/2019
# Desc: import clean and model for the differences in rd.score for pasi90 trus/false groups

############################## data import and cleaning steps

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
# check the data for na values
dim(dfData)
f = is.na(dfData$Rd.score)
table(f)
dfData = dfData[!f,]

dim(dfData)

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

## choose only the transcription factor nfkb
levels(dfData$Transcription.factor)
dfData = dfData[dfData$Transcription.factor == 'NF-kB', ]
dfData = droplevels.data.frame(dfData)
dim(dfData)

## drop some of the stimulations as decided 
levels(dfData$Stimulation)
dfData = dfData[!(dfData$Stimulation %in% c('IL-17', "TNF-alpha + IL-17")), ]
dfData = droplevels.data.frame(dfData)
dim(dfData)


####### apply the cutoffs
## cell count variable check this 
dim(dfData)
dfData = dfData[dfData$Cell.count >= 10,]
dim(dfData)

# define the blocks as used by previous analyst
# in each cell type, the stimulation and transcription factors are nested 
# apart from where Stimulation = Unstimulated
nlevels(factor(dfData$Stimulation:dfData$Transcription.factor:dfData$Cell.type))
nlevels(factor(dfData$Stimulation:dfData$Cell.type))
nlevels(factor(dfData$Stimulation:dfData$Transcription.factor))
nlevels(dfData$Stimulation)
xtabs(~ dfData$Transcription.factor+dfData$Cell.type)
xtabs(~ dfData$Transcription.factor+dfData$Stimulation)
# module is cell type  + stimulation
fModule = factor(dfData$Cell.type:dfData$Stimulation)
nlevels(fModule)
# block is a combination of modules and time
fBlock = factor(fModule:dfData$Visit..Week.)
nlevels(fBlock)
## 80 blocks
## drop modules with average RD score < 0.3
nlevels(fModule)
i = tapply(dfData$Rd.score, fModule, mean)
summary(i)
i = which(i < 0.3)
i = levels(fModule)[i]
# drop these modules from the dataset
f = which(fModule %in% i)
dim(dfData)
dfData = dfData[-f,]
dim(dfData)
dfData = droplevels.data.frame(dfData)

# make modules and blocks again
fModule = factor(dfData$Cell.type:dfData$Stimulation)
nlevels(fModule)
# block is a combination of modules and time
fBlock = factor(fModule:dfData$Visit..Week.)
nlevels(fBlock)

dfData$fModule = fModule
dfData$fBlock = fBlock

xtabs(~ dfData$fModule  + dfData$Transcription.factor)
xtabs(~ dfData$fBlock + dfData$Transcription.factor)
rm(fModule)
rm(fBlock)

write.csv(dfData, file='dataExternal/healthyData/pasi90VsRdScoreAdalimumab_onlyNFKB_noIL17.csv', row.names = F)

########################### data modelling steps

library(lattice)

gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

#dfData = read.csv('dataExternal/healthyData/pasi90VsRdScoreAdalimumab.csv', header=T)

# dfData$Visit..Week. = factor(dfData$Visit..Week., levels=c('Baseline', 
#                                                            'Week 1', 'Week 4',
#                                                            'Week 12'))

dotplot(PASI90 ~ Rd.score | fBlock, data=dfData, type=c('p'),
       ##index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), pch=20)

dotplot(PASI90 ~ Rd.score | fBlock, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7))

bwplot(PASI90 ~ Rd.score | fBlock, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7), varwidth=T)

bwplot(PASI90 ~ Rd.score | fBlock, data=dfData, panel=panel.violin, type='b',
       par.strip.text=list(cex=0.7), varwidth=F)

## use data on original scale
densityplot(dfData$Rd.score, groups=dfData$fBlock)
densityplot(dfData$Rd.score)

## create a new factor with a combinations of factors of interest
nlevels(dfData$fBlock)
dfData$PASI90 = factor(dfData$PASI90)
levels(dfData$PASI90)
f = factor(dfData$fBlock:dfData$PASI90)
nlevels(f)
table(f)
dfData$Coef = f
densityplot(dfData$Rd.score, groups=dfData$Coef)
densityplot(dfData$Rd.score)

dfData = dfData[order(dfData$Patient.ID, dfData$Coef),]
dfData = droplevels.data.frame(dfData)
str(dfData)

#### fit mixed effect model
library(lme4)
fit.lme1 = lmer(Rd.score ~ 1 + (1 | Coef) + (1 | Patient.ID), data=dfData, REML=F)
summary(fit.lme1)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponse2RandomEffectsNoFixed.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$Rd.score), 2*sd(dfData$Rd.score))

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef), Nclusters2=nlevels(dfData$Patient.ID),
                 NgroupMap1=as.numeric(dfData$Coef), NgroupMap2=as.numeric(dfData$Patient.ID),
                 Ncol=1, 
                 y=dfData$Rd.score, gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=2, pars=c('betas', 'sigmaRan1', 'sigmaRan2', 'nu',
                                                                          'sigmaPop','mu', 'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)


##############################
########## model checks section
###### some model checks
### check model fit
m = extract(fit.stan, 'mu')
names(m)
dim(m$mu)
fitted = apply(m$mu, 2, mean)

plot(dfData$Rd.score, fitted, pch=20, cex=0.5)
plot(dfData$Rd.score, dfData$Rd.score - fitted, pch=20, cex=0.5)
iResid = (dfData$Rd.score - fitted)

par(mfrow=c(1,2))
plot(fitted, iResid, pch=20, cex=0.5, main='t model')
lines(lowess(fitted, iResid), col=2, lwd=2)

plot(predict(fit.lme1), resid(fit.lme1), pch=20, cex=0.5, main='normal')
lines(lowess(predict(fit.lme1), resid(fit.lme1)), col=2, lwd=2)

plot(fitted, predict(fit.lme1), pch=20, cex=0.5)
abline(0, 1, col='red')


### plot the posterior predictive values
m = extract(fit.stan, c('mu', 'nu', 'sigmaPop'))
i = sample(1:10000, 2000)
muSample = m$mu[i,]
nuSample = m$nu[i]
sigSample = m$sigmaPop[i]

## t sampling functions
dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu

## use point-wise predictive approach to sample a new value from this data
ivResp = dfData$Rd.score
mDraws = matrix(NA, nrow = length(ivResp), ncol=2000)

# rppd = function(index){
#   f = muSample[,index]
#   return(rt_ls(length(f), nuSample, f, sigSample))
# }

for (i in 1:ncol(mDraws)){
  mDraws[,i] = rt_ls(length(ivResp), nuSample[i], muSample[i,], sigSample[i])
}

# 
# temp = sapply(1:length(ivResp), function(x) rppd(x))
# mDraws = t(temp)
par(mfrow=c(1,1))
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

################ end model checks
## get the coefficient of interest - block:pasi90 in our case from the random coefficients section
mModules = extract(fit.stan)$rGroupsJitter1
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
d = data.frame(cols=1:ncol(mModules), mods=levels(dfData$Coef), rawAverage=tapply(dfData$Rd.score, dfData$Coef, FUN = mean))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:3], c('cells', 'stimulation', 'time', 'PASI90'))
## create the comparisons/contrasts required
d$split = factor(d$cells:d$stimulation:d$time)
head(d)

########## contrast of interest
### contrasts within each block
## this data frame is a mapper for each required comparison
ldfMap = split(d, f = d$split)

## get a p-value for each comparison
l = lapply(ldfMap, function(x) {
  c = x$cols
  d = getDifference(ivData = mModules[,c[2]], ivBaseline = mModules[,c[1]])
  r = data.frame(block= as.character(x$split[1]), coef.PASI90.FALSE=mean(mModules[,c[1]]), 
                 coef.PASI90.TRUE=mean(mModules[,c[2]]), zscore=d$z, pvalue=d$p, rawAverage.PASI90.FALSE=x$rawAverage[1],
                 rawAverage.PASI90.TRUE=x$rawAverage[2])
  return(format(r, digi=3))
})

dfResults = do.call(rbind, l)
dfResults$p.adj = format(p.adjust(dfResults$pvalue, method='bonf'), digi=3)
write.csv(dfResults, file='Results/pasi90vsRdScoreAdalimumab_onlyNFKB_noIL17.csv', row.names = F)

########################## continue from here to make plots of coefficients
## make the plots for the raw data and coef
bwplot(PASI90 ~ Rd.score | fBlock, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...),
        par.strip.text=list(cex=0.5), varwidth=T)

## format data for plotting
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Coef  )
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('cells', 'stimulation', 'time', 'PASI90'))
d$time = factor(d$time, levels=c('Baseline', 'Week 1', 'Week 4', 'Week 12'))

dotplot(PASI90 ~ m+s1+s2 | cells:stimulation:time, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=0.5), main='48 Blocks with Fitted Coefficients for Rd.score ~ PASI90', xlab='')


## example of xyplot with confidence interval bars
# xyplot(m ~ treatment | cells:stimulation, data=d, 
#        panel=function(x,y,ulim,llim, ...) { 
#          lsegments(x, y-d$s, x, y+d$s, lwd=1)
#          panel.xyplot(x,y, ...) 
#        } 
#        , type=c('p'), pch=19, groups=treatment, cex=0.5,
#        #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
#        par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules - baseline',
#        auto.key = list(columns=3), ylim=c(min(d$s2), max(d$s1)))




