# Name: 03_mergedDataModel.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 19/09/2017
# Desc: imports the merged data and fits a model

library(lattice)

gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

dfData = read.csv('dataExternal/healthyData/mergedData.csv', header=T)

## use only baseline 
dfData = dfData[dfData$Visit..Week. == 'Baseline',]
dfData = droplevels.data.frame(dfData)
dfData$Treatment = relevel(dfData$Treatment, 'None')
dfData$Modules = factor(dfData$fModule:dfData$Treatment)
table(dfData$Modules)
dfData = dfData[order(dfData$Patient.ID, dfData$Modules),]
str(dfData)

## make a plot of the raw data
## format before plotting
d2 = dfData[,c('Rd.score', 'Modules')]
f = strsplit(as.character(d2$Modules), ':')
d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('cells', 'stimulation', 'treatment'))

dotplot(treatment ~ Rd.score | cells:stimulation, data=d2, groups=treatment, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...),
        par.strip.text=list(cex=0.6), main='Raw data 41 Modules at baseline', xlab='Raw RD Score')

## log the data before modelling
ivResp = dfData$Rd.score
iShift = min(ivResp)+1
ivResp = ivResp+abs(min(ivResp))+1
ivResp = log(ivResp)
dfData$Rd.score = ivResp


library(lme4)
fit.lme1 = lmer(Rd.score ~ 1 + (1 | Modules) + (1 | Patient.ID), data=dfData, REML=F)
summary(fit.lme1)

# fit.lme2 = lmer(Rd.score ~ 1 + Transcription.factor + (1 | Patient.ID) + (1 | Modules), data=dfData, REML=F)
# summary(fit.lme2)
# 
# fit.lme3 = lmer(Rd.score ~ 1 + Transcription.factor + Treatment + (1 | Patient.ID) + (1 | Modules), data=dfData, REML=F)
# summary(fit.lme3)
# 
# anova(fit.lme1, fit.lme2, fit.lme3)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='tResponse2RandomEffectsNoFixed.stan')

l = gammaShRaFromModeSD(sd(dfData$Rd.score), 2*sd(dfData$Rd.score))
# m = model.matrix(Rd.score ~ Transcription.factor, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Patient.ID), Nclusters2=nlevels(dfData$Modules), 
                 NgroupMap1=as.numeric(dfData$Patient.ID), NgroupMap2=as.numeric(dfData$Modules), 
                 Ncol=1, #X=m,
                 y=dfData$Rd.score, gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('betas', 'sigmaRan1', 'sigmaRan2',
                                                                         'sigmaPop','nu', 'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=4)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)

## get the coefficient of interest - Modules in our case
mModules = extract(fit.stan)$rGroupsJitter2
dim(mModules)
iIntercept = extract(fit.stan)$betas[,1]
## add the intercept to each random effect variable
mModules = sweep(mModules, 1, iIntercept, '+')

## make the plots for the raw data and fitted data
## format data for plotting
d2 = dfData[,c('Rd.score', 'Modules')]
f = strsplit(as.character(d2$Modules), ':')
d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('cells', 'stimulation', 'treatment'))

dotplot(treatment ~ Rd.score | cells:stimulation, data=d2, groups=treatment, panel=function(x, y, ...) panel.bwplot(x, y, pch='|', ...),
        par.strip.text=list(cex=0.6), main='Raw data 41 Modules at baseline', xlab='Log RD Score')

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
        par.strip.text=list(cex=0.6), main='Regression Coeff 41 Modules at baseline', xlab='Model estimated Average Log RD Score')


## convert coefficients to natural scale / i.e. exponent
## format data for plotting
m2 = exp(mModules)-iShift
m = colMeans(m2)
s = apply(m2, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Modules)
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('cells', 'stimulation', 'treatment'))

dotplot(treatment ~ m+s1+s2 | cells:stimulation, data=d, panel=llines(d$s1, d$s2), cex=0.6,
        par.strip.text=list(cex=0.7), main='Regression Coeff 41 Modules at baseline', xlab='Model estimated Average RD Score')


### get p.values for contrasts of interest





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




