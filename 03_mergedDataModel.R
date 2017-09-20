# Name: 03_mergedDataModel.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 19/09/2017
# Desc: imports the merged data and fits a model


gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

dfData = read.csv('dataExternal/healthyData/mergedData.csv', header=T)
dfData = dfData[order(dfData$Patient.ID, dfData$fModule),]
ivResp = dfData$Rd.score
ivResp = ivResp+abs(min(ivResp))+1
ivResp = log(ivResp)
dfData$Rd.score = ivResp
dfData = dfData[dfData$Visit..Week. == 'Baseline',]
dfData = droplevels.data.frame(dfData)
dfData$Treatment = relevel(dfData$Treatment, 'None')
str(dfData)

xtabs( ~ Stimulation + Transcription.factor, data=dfData)
xtabs( ~ fModule + Transcription.factor, data=dfData)


library(lme4)
fit.lme1 = lmer(Rd.score ~ Treatment + (1 | Patient.ID) + (1 | fModule) + (1 | Transcription.factor), data=dfData, REML=F)
summary(fit.lme)

fit.lme2 = lmer(Rd.score ~ Treatment + fModule + (1 | Patient.ID) +  (1 | Transcription.factor), data=dfData, REML=F)
summary(fit.lme2)

fit.lme3 = lmer(Rd.score ~ Treatment:fModule + (1 | Patient.ID), data=dfData, REML=F)
summary(fit.lme3)

fit.lme4 = lmer(Rd.score ~ Treatment+fModule+Treatment:fModule + (1 | Patient.ID), data=dfData, REML=F)
summary(fit.lme4)


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='tResponseRandomEffects.stan')

l = gammaShRaFromModeSD(sd(dfData$Rd.score), 2*sd(dfData$Rd.score))
m = model.matrix(Rd.score ~ Treatment+fModule+Treatment:fModule, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Patient.ID), #Nclusters2=nlevels(dfData$fModule), 
                 NgroupMap1=as.numeric(dfData$Patient.ID), #NgroupMap2=as.numeric(dfData$fModule), 
                 Ncol=ncol(m), X=m,
                 y=dfData$Rd.score, gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('betas', 'sigmaRan1', #'sigmaRan2',
                                                                          'sigmaPop','nu', 'rGroupsJitter1'), #'rGroupsJitter2'),
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, digits=3)






