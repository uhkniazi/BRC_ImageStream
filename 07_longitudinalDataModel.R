# Name: 07_longitudinalDataModel.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 09/10/2017
# Desc: imports the longitudinal data and fits a model

library(lattice)

gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

dfData = read.csv('dataExternal/healthyData/diseasedData.csv', header=T)

dfData$Visit..Week. = factor(dfData$Visit..Week., levels=c('Baseline', 
                                                           'Week 1', 'Week 4',
                                                           'Week 12'))

xyplot(Rd.score ~ time | fModule, data=dfData, type=c('smooth'),
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), groups=dfData$Treatment,
       auto.key = list(columns=2))

xyplot(Rd.score ~ time | fModule, data=dfData, type=c('r'),
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), groups=dfData$Treatment,
       auto.key = list(columns=2))

xyplot(Rd.score ~ Visit..Week. | fModule, data=dfData, type=c('p'), pch=20,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), groups=dfData$Treatment,
       auto.key = list(columns=2))

dotplot(Treatment ~ Rd.score | fModule:Visit..Week., data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.5))

dotplot(Treatment ~ Rd.score | fModule:Visit..Week., data=dfData[dfData$Visit..Week. == 'Baseline',], panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7))

dotplot(Treatment ~ Rd.score | fModule:Visit..Week., data=dfData[dfData$Visit..Week. == 'Week 1',], panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7))

dotplot(Treatment ~ Rd.score | fModule:Visit..Week., data=dfData[dfData$Visit..Week. == 'Week 4',], panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7))

dotplot(Treatment ~ Rd.score | fModule:Visit..Week., data=dfData[dfData$Visit..Week. == 'Week 12',], panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7))

# i = sample(1:nrow(dfData), 1000)
# dfData = dfData[i,]
# dim(dfData)

#dotplot(Treatment ~ Rd.score | Visit..Week.:Cell.type, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b')
## make a plot of the raw data
## format before plotting
# d2 = dfData[,c('Rd.score',  'Treatment', 'fModule')]
# f = strsplit(as.character(d2$fModule), ':')
# d2 = cbind(d2, do.call(rbind, f))
# colnames(d2) = c(colnames(d2)[1:3], c('cells', 'stimulation'))
# 
# dotplot(Treatment ~ Rd.score | cells:stimulation, data=d2, groups=Treatment, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...),
#         par.strip.text=list(cex=0.6), main='Raw data 41 Modules at baseline', xlab='Raw RD Score')


## log the data before modelling
ivResp = dfData$Rd.score
iShift = min(ivResp)+2
ivResp = ivResp+abs(min(ivResp))+2
ivResp = log(ivResp)
dfData$Rd.score = ivResp

## create a new factor with a combinations of factors of interest
f = factor(dfData$Treatment:dfData$fModule:dfData$Visit..Week.)
dfData$Coef = f
densityplot(dfData$Rd.score, groups=dfData$Coef)
densityplot(dfData$Rd.score)

library(flexmix)
fit.flex = flexmix(Rd.score ~ Coef, data=dfData, k=2)
#fit.flex = flexmix(.~.|Treatment, data=dfData, k=2, model=FLXMRlmer(Rd.score ~ Treatment, random=~ 1))
summary(fit.flex)
## fitted coefficients
parameters(fit.flex)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='fitNormMixtureRandomEffects.stan')

#l = gammaShRaFromModeSD(sd(dfData$Rd.score), 2*sd(dfData$Rd.score))
# m = model.matrix(Rd.score ~ Transcription.factor, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef), Nclusters2=nlevels(dfData$Patient.ID), 
                 NgroupMap1=as.numeric(dfData$Coef), NgroupMap2=as.numeric(dfData$Patient.ID),  
                 y=dfData$Rd.score, iMixtures=2)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, cores=2, 
                    pars=c('mu', 'sigma', 'iMixWeights',
                           'rGroupsJitter1', 'sigmaRan1',
                           'rGroupsJitter2', 'sigmaRan2'))  ####, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('mu', 'sigma', 'iMixWeights',
                  'sigmaRan1', 'sigmaRan2'), digits=3)

## get the coefficient of interest - Modules in our case from the random coefficients section
mModules = extract(fit.stan)$rGroupsJitter1
dim(mModules)
## get the intercept at population level
#iIntercept = extract(fit.stan)$betas[,1]
## add the intercept to each random effect variable, to get the full coefficient
#mModules = sweep(mModules, 1, iIntercept, '+')

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
d = data.frame(cols=1:ncol(mModules), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:2], c('drug', 'cells', 'stimulation', 'time'))
d$split = factor(d$cells:d$stimulation:d$time)

## this data frame is a mapper for each required comparison
ldfMap = split(d, f = d$split)

## get a p-value for each comparison
l = lapply(ldfMap, function(x) {
  c = x$cols
  d = getDifference(ivData = mModules[,c[2]], ivBaseline = mModules[,c[1]])
  r = data.frame(module= as.character(x$split[1]), coef.adal=mean(mModules[,c[1]]), 
        coef.ustek=mean(mModules[,c[2]]), zscore=d$z, pvalue=d$p)
  return(format(r, digi=3))
})

dfResults = do.call(rbind, l)
dfResults$p.adj = format(p.adjust(dfResults$pvalue, method='bonf'), digi=3)
write.csv(dfResults, file='Results/longitudinalDataResults.csv', row.names = F)

## make the plots for the raw data and fitted data
## format data for plotting
d2 = dfData[,c('Rd.score', 'Coef')]
#f = strsplit(as.character(d2$Modules), ':')
f = strsplit(as.character(d2$Coef), ':')
d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('drug', 'cells', 'stimulation', 'time'))

d2 = cbind(d2, do.call(rbind, f))
colnames(d2) = c(colnames(d2)[1:2], c('cells', 'stimulation', 'treatment'))

dotplot(drug ~ Rd.score | cells:stimulation:time, data=d2, groups=drug, panel=function(x, y, ...) panel.bwplot(x, y, pch='|', ...),
        par.strip.text=list(cex=0.6), main='Raw data 164 Comparisons', xlab='Log RD Score')

## format data for plotting
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Coef)
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('drug', 'cells', 'stimulation', 'time'))

dotplot(drug ~ m+s1+s2 | cells:stimulation:time, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=0.6), main='328 Regression Coeff 164 Comparisons', xlab='Model estimated Coefficients Log RD Score')


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




