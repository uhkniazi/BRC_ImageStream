# Name: 05_mergedDataModelUnstimulated.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 28/09/2017
# Desc: imports the merged data for unstimulated samples and fits a model

library(lattice)

gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

dfData = read.csv('dataExternal/healthyData/mergedDataUnstimulated.csv', header=T)

## use only baseline samples and merge 2 types of treated into one group
t = as.character(dfData$Treatment)
t[t != 'None'] = 'Patient'
t[t == 'None'] = 'Healthy'
dfData$Treatment = factor(t)
dfData$Modules = factor(dfData$fModule:dfData$Transcription.factor:dfData$Treatment)
nlevels(dfData$Modules); table(dfData$Modules)
dfData = dfData[order(dfData$Patient.ID, dfData$Modules),]
dfData = droplevels.data.frame(dfData)
str(dfData)

## make a plot of the raw data
## format before plotting
d2 = dfData[,c('Median.internalization.score', 'Modules')]
f = strsplit(as.character(d2$Modules), ':')
d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('cells', 'transcription.factor', 'treatment'))

dotplot(treatment ~ Median.internalization.score | cells:transcription.factor, data=d2, groups=treatment, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...),
        par.strip.text=list(cex=0.6), main='Raw data 30 Modules at baseline', xlab='Raw Med int Score')

## shift the data before modelling
ivResp = dfData$Median.internalization.score
iShift = min(ivResp)
ivResp = ivResp+abs(min(ivResp))+0.1
dfData$Median.internalization.score = ivResp


library(lme4)
fit.lme1 = glmer(Median.internalization.score ~ 1 + (1 | Modules) + (1 | Patient.ID), data=dfData, family=Gamma(link='identity'))
summary(fit.lme1)

fit.lme2 = glmer(Median.internalization.score ~ 1 + (1 | Modules) + (1 | Patient.ID), data=dfData, family=Gamma(link='log'))
summary(fit.lme2)

anova(fit.lme1, fit.lme2)


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='gammaResponse2RandomEffectsNoFixed.stan')

l = gammaShRaFromModeSD(sd(log(dfData$Median.internalization.score)), 2*sd(log(dfData$Median.internalization.score)))
# m = model.matrix(Median.internalization.score ~ Transcription.factor, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Patient.ID), Nclusters2=nlevels(dfData$Modules), 
                 NgroupMap1=as.numeric(dfData$Patient.ID), NgroupMap2=as.numeric(dfData$Modules), 
                 Ncol=1, #X=m,
                 y=dfData$Median.internalization.score, gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=4000, chains=4, pars=c('betas', 'sigmaRan1', 'sigmaRan2',
                                                                         'sigmaPop','alpha', 'beta', 'mu', 'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=4)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop'), digits=3)

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
d = data.frame(cols=1:ncol(mModules), mods=levels(dfData$Modules))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:2], c('cells', 'transcription.factor', 'treatment'))
d$split = factor(d$cells:d$transcription.factor)

## this data frame is a mapper for each required comparison
ldfMap = split(d, f = d$split)

## get a p-value for each comparison
l = lapply(ldfMap, function(x) {
  c = x$cols
  d = getDifference(ivData = mModules[,c[2]], ivBaseline = mModules[,c[1]])
  r = data.frame(cells= as.character(x$cells[1]), transcription.factor=as.character(x$transcription.factor[1]), coef.healthy=mean(mModules[,c[1]]), 
        coef.patient=mean(mModules[,c[2]]), zscore=d$z, pvalue=d$p)
  return(format(r, digi=3))
})

dfResults = do.call(rbind, l)
dfResults$p.adj = format(p.adjust(dfResults$pvalue, method='bonf'), digi=3)
write.csv(dfResults, file='Results/mergedDataResultsUnstimulated.csv', row.names = F)

## make the plots for the raw data and fitted data
## format data for plotting
d2 = dfData[,c('Median.internalization.score', 'Modules')]
f = strsplit(as.character(d2$Modules), ':')
d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('cells', 'transcription.factor', 'treatment'))

dotplot(treatment ~ Median.internalization.score | cells:transcription.factor, data=d2, groups=treatment, panel=function(x, y, ...) panel.bwplot(x, y, pch='|', ...),
        par.strip.text=list(cex=0.6), main='Raw data 30 Modules at baseline', xlab='Shifted Med int Score')

## format data for plotting
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Modules)
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('cells', 'transcription.factor', 'treatment'))

dotplot(treatment ~ m+s1+s2 | cells:transcription.factor, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=0.6), main='Regression Coeff 30 Modules at baseline', xlab='Model estimated Average Log Med int Score')


## convert coefficients to natural scale / i.e. exponent
## format data for plotting
m2 = exp(mModules)
m = colMeans(m2)
s = apply(m2, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Modules)
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('cells', 'transcription.factor', 'treatment'))

dotplot(treatment ~ m+s1+s2 | cells:transcription.factor, data=d, panel=llines(d$s1, d$s2), cex=0.6,
        par.strip.text=list(cex=0.7), main='Regression Coeff 30 Modules at baseline', xlab='Model estimated Average Med int Score')







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




