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

xyplot(Rd.score ~ Stimulation | Cell.type, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

xyplot(Rd.score ~ Cell.type | Stimulation, data=dfData, type=c('g', 'p'), pch=19,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)))

# check data distribution
x = na.omit(dfData$Rd.score)
hist(x)
fitdistr(x, densfun = 't')
qqPlot(x, distribution = 't', df=3)
