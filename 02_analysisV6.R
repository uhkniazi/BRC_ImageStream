## DESIGN OF THE EXPERIMENT:

#- ~50  patients: some from London ( Ustek,  Adalim), some from Queen Elizabeth (Ustek), some from West Middlesex (Ustek)
#- Blood taken from each patient at 4 different time points: before treatment, 1w, 4w and 12w of treatment
#- each sample gives rise to the following experimental conditions:
#  ==> NFkB: unstimulated, TNFkb, TNFkb + IL17, IL17, LPS
#  ==> STAT3: unstimulated, IL6, IL6 + IL23, IL23
#  ==> STAT4: unstimulated, IL12, IFNgamma
#- each of these experimental conditions is tested in 10 different cell types (CD4+ T-cells, CD- T-cells etc.)


#  Q1: WHAT ARE THE TFs/STIMULI/CELL TYPE/TREATMENT WHERE TF TRANSLOCATION (Rd.Score) CORRELATES WITH RESPONSE (relative PASI, PASI75, PASI90, extreme categories)

#  Q2: ARE TREATMENTS HAVING MOLECULAR EFFECTS ?
# 
## Q3: ARE TREATMENTS HAVING DIFFERENT EFFECTS ? IF SO, IN WHICH TF/CELL TYPE/STIMULI 

## Q4: Correlation between our translocation scores and the drug levels (especially for TNFa and adalimumab). 
# 
# The key question regarding the antidrug antibodies is if patients 006012 and 006013 behave differently to the rest of the adalimumab cohort, 
# as those 2 are the only ones with detectable ADAs at week 12 (12 and 30 are the limit of detection to the corresponding technologies: ABT and ARIA).

##EXTREME CASES:
# For the extreme cases, the non-responders would be those who don't meet PASI50 and the super-responders would be those with an absolute PASI below 1.5 at week12.

## ADDITIONAL QUESTIONS

# Cell count: I chose the value of 40 based on my experience and not to compromise too many samples but we need to revisit it: 
# what is the impact of cell count in the median internalization score/MAD/rd score? Is 40 too low or could we lower it further?
# Is the Rd score a suitable measurement for our data? I used it because of the recommendations of the ISX developers but I never had input from a statistician about 
# its suitability or if there is any other better measurement for us.



#install.packages("PerformanceAnalytics")
#install.packages("corrplot")
#install.packages("ggplot2")
#install.packages("plyr")
#install.packages("lattice")
#install.packages("contrast")
#install.packages("corrplot")
#install.packages("devtools")
#install_github("easyGgplot2", "kassambara")
#install.packages("multcomp")
#install.packages("phia")
#install.packages("outliers")



library(phia)
library("multcomp")
library("contrast")
library(lattice)
library(plyr)
library(ggplot2)
library("PerformanceAnalytics")
source("http://www.sthda.com/upload/rquery_cormat.r")
library(nlme)
library(devtools)
library(easyGgplot2)
library(outliers)

rm(list = ls())  
 
setwd("C://Users/derinald/Desktop/ImageStream Data Set - October 2016/data")
main_path = "C://Users/derinald/Desktop/ImageStream Data Set - October 2016/"
input_dir=  paste(main_path, "data", sep="/")
out_dir=paste(main_path, "R/results/", sep="/")

dataset<-read.table(file = "merged.files_and_annotations.txt", header=TRUE)

############################################
# PREPARATION OF WORKING DATASET          ##
############################################
ds=dataset[!is.na(dataset$Rd.score) | dataset$Stimulation=="Unstimulated",]
ds$tp[ds$Visit..Week== "Baseline"] = "one"
ds$tp[ds$Visit..Week== "Week 1"]   = "two"
ds$tp[ds$Visit..Week== "Week 4"]   = "three"
ds$tp[ds$Visit..Week== "Week 12"]  = "four"
ds$tp<-factor(ds$tp)
ds$GenderId <- factor(ds$Gender)      
ds$Ethnicity <- factor(ds$Ethnicity)   
ds$TF_STIMULATION <- factor(paste(ds$Transcription.factor, ds$Stimulation, sep="_"))
ds$ExtremeResponse=NA
ix=ds$Week.12.PASI<1.5
ds$ExtremeResponse[ix] <- 1 
ix <- ds$PASI50==FALSE
ds$ExtremeResponse[ix] <- 0 

ds$Block1 <- factor(paste(ds$tp, ds$Cell.type, ds$TF_STIMULATION, sep="_"))      ##treatment is applied only once in each block
ds$Block2 <- factor(paste(ds$tp, ds$Cell.type, ds$TF_STIMULATION, sep="_"))
ds$Block3 <- factor(paste(ds$tp, ds$Cell.type, ds$TF_STIMULATION, ds$Treatment, sep="_"))  
ds$Block4 <- factor(paste(ds$Cell.type, ds$Stimulation, ds$Transcription.factor, ds$Treatment, ds$Visit..Week ,sep="__"))    
ds$Module <- factor(paste(ds$Cell.type, ds$Stimulation, sep="_"))      
ds$relativePASI.log10<-log10(ds$relativePASI + 1e-16)

      
c<-1
for(i in unique(ds$Block1))
{
  ds$Block1.num[ds$Block1==i]<-c;c<-c+1
}
ds$Block1.num<-as.numeric(ds$Block1.num)
c<-1
for(i in unique(ds$Block2))
{
  ds$Block2.num[ds$Block2==i]<-c;c<-c+1
}
ds$Block2.num<-as.numeric(ds$Block2.num)
c<-1
for(i in unique(ds$Block3))
{
  ds$Block3.num[ds$Block3==i]<-c;c<-c+1
}
ds$Block3.num<-as.numeric(ds$Block3.num)
ds.bkp=ds
 
###
###
CC_THRES=10
#####################################################################################
### "MODULE" ANALYSIS. EACH MODULE CONSISTS OF 8 BLOCKS (CELL.TYPE+STIMULATION).   ##
##   CHECKS THAT MEAN RD.SCORE OF AT LEAST ONE BLOCK OF THE MODULE IS > 0.3       ##
#####################################################################################

ix=which(ds$Cell.count > CC_THRES & ds$Stimulation!="Unstimulated")
ds.tmp=ds[ix,]
module.check<-list()
for(i in  levels(ds.tmp$Module))
{
   ds.m<-ds.tmp[ds.tmp$Module==i,] 
   module.check[i]<-FALSE
   for(j in unique(ds.m$Block4))
   {
    ds.b<-ds.m[ds.m$Block4==j,]
    module.check[i]<-max(module.check[[i]],(mean(ds.b$Rd.score,na.rm = TRUE)>0.3))
    print(i)
    print(j)
    print(module.check[i])
    #print(mean(ds.b$Rd.score,na.rm = TRUE))
    #print(module.check[i])
   } 
}     

print(sum(module.check==FALSE))
print(sum(module.check==TRUE))
print(names(module.check)[module.check==TRUE])
a<-names(module.check)[module.check==TRUE]
b<-as.character(unique(ds[ds$Stimulation=="Unstimulated","Module"]))
selected.modules<-c(a,b)
ds.m.filter<-ds[ds$Module %in% selected.modules, ]
#### 
####
##############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
# 
 
 
   
##############################################################################################################################################################################
## ESTABLISHES POSSIBLE CELL COUNT THRESHOLDS FOR EACH CELL TYPE USING DIFFERENCE CRITERIA E.G. MIN BETWEEN 200 (ARBITRARY VALUE) AND THE 10TH PERCENTILE OF CELL COUNTS   ##
##############################################################################################################################################################################
##
ct.thres1<-list()
ct.thres2<-list()
ix<-list()
for(i in unique(ds$Cell.type))
{
 print(i)
 print(c("mean = ", mean(ds[ds$Cell.type==i,"Cell.count"], na.rm=TRUE)))
 print(c("std = ",sd(ds[ds$Cell.type==i,"Cell.count"], na.rm=TRUE)))
 print(c("threshold = ", min(quantile(ds[ds$Cell.type==i,"Cell.count"], probs=seq(0,1,0.05))[5],200)))
 #ct.thres[i]<-min(quantile(ds[ds$Cell.type==i,"Cell.count"], probs=seq(0,1,0.05))[3],200)
 ct.thres1[i]<-min(quantile(ds[ds$Cell.type==i,"Cell.count"], probs=seq(0,1,0.05))[5],200)
 ct.thres2[i]<-min(quantile(ds[ds$Cell.type==i,"Cell.count"], probs=seq(0,1,0.05))[5],200)
 
 ix[[i]]<-which(ds$Cell.type==i & ds$Cell.count<ct.thres1[i])
 #outl[i]<-outlier(ds[ds$Cell.type==i,"Cell.count"], opposite = FALSE, logical = FALSE)
}

plot(ds$Median.internalization.score,   ds$Rd.score)
cor(ds$Median.internalization.score,   ds$Rd.score, use="complete.obs")

ds.excluded.exps<-ds[unlist(ix),]
ds2<- ds[-p,]
plot(ds2$Median.internalization.score,   ds2$Rd.score)
cor(ds2$Median.internalization.score,   ds2$Rd.score, use="complete.obs")
plot(ds.excluded.exps$Cell.type)

ds3<- ds[ds$Cell.count>10,]
plot(ds3$Median.internalization.score,   ds3$Rd.score)
cor(ds3$Median.internalization.score,   ds3$Rd.score, use="complete.obs")

ds4<- ds[ds$Cell.type=="NKT cells",]
plot(ds4$Median.internalization.score,   ds4$Rd.score)
cor(ds4$Median.internalization.score,   ds4$Rd.score, use="complete.obs")

thres<-c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500)
a<-list()
for(i in thres)
{
  print(i)
  ds.tmp<-ds[ds$Cell.count>i,]
  a[i]<-cor(ds.tmp$Median.internalization.score,   ds.tmp$Rd.score, use="complete.obs")
  print(a[i])
}
b=unlist(a)
plot(b)



###




#### Analysis of how median Internalization Score (MIS) variates for each cell type, as function of the number of cells of each experiment
#### These results are used as a guide to establish a specific cell count threshold for each cell type

par(mfcol=c(1,2))
hist(ds$Cell.count, main = "ALL") ## cell count distribution across the all dataset
plot(ds$Cell.count, main = "ALL") ## cell count distribution across the all dataset

par(mfcol=c(4,3))
for(i in unique(ds$Cell.type))
{
  hist(ds$Cell.count[ds$Cell.type == i], main=i, n=10)
}


for(j in unique(ds$Visit..Week))
{
  dev.new()
  par(mfcol=c(4,3))
  for(i in unique(ds$Cell.type))
  {
    plot(ds$Cell.count[ds$Cell.type == i & ds$Visit..Week == j], ds$Median.internalization.scor[ds$Cell.type == i & ds$Visit..Week == j],main=i, xlab="CELL COUNT", ylab="MIS")
  }
  title(j, outer=TRUE)
}

## PLOTTING MIS VS BINNED CELL COUNTS FOR ALL EXPERIMENTS

v=list()
q=quantile(ds$Cell.count, probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))  
for(k in c(1:length(q)))    ## variability of MIS at each bin, in each cell type, at each time point
{
  if(k==1)
  {
    ix=which(ds$Cell.count < q[k] & ds$Cell.type == i & ds$Visit..Week == j)
  } else{ix=which(ds$Cell.count < q[k] & ds$Cell.count > q[k-1])}
  v[[k]]=ds$Median.internalization.scor[ix]
}
boxplot(v, names=q, col=2, xlab="Cell Count", ylab="MIS", main="ALL CELL TYPES ACROSS ALL TIME POINTS")



## PLOTTING MIS VS BINNED CELL COUNTS (ONE BLOCK AT A TIME)
for(j in unique(ds$Visit..Week))
{
  dev.new()
  par(mfcol=c(4,3))
  for(i in unique(ds$Cell.type))
  {
    v=list()
    q=quantile(ds$Cell.count[ds$Cell.type == i & ds$Visit..Week == j], probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))  
    for(k in c(1:length(q)))    ## variability of MIS at each bin, in each cell type, at each time point
    {
      if(k==1)
      {
        ix=which(ds$Cell.count < q[k] & ds$Cell.type == i & ds$Visit..Week == j)
      } else{ix=which(ds$Cell.count < q[k] & ds$Cell.count > q[k-1] & ds$Cell.type == i & ds$Visit..Week == j)}
      v[[k]]=ds$Median.internalization.scor[ix]
      #print(v[k])
    }
    boxplot(v, names=q, col=2, xlab="Cell Count", ylab="MIS", main=i)
  }
  title(j, outer=TRUE)
}

## PLOTTING MIS VARIANCES VS BINNED CELL COUNTS  (ONE BLOCK AT A TIME)
for(j in unique(ds$Visit..Week))
{
  dev.new()
  par(mfcol=c(4,3))
  for(i in unique(ds$Cell.type))
  {
    v=list()
    q=quantile(ds$Cell.count[ds$Cell.type == i & ds$Visit..Week == j], probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))  
    for(k in c(1:length(q)))    ## variability of MIS at each bin, in each cell type, at each time point
    {
      if(k==1)
      {
        ix=which(ds$Cell.count < q[k] & ds$Cell.type == i & ds$Visit..Week == j)
      } else{ix=which(ds$Cell.count < q[k] & ds$Cell.count > q[k-1] & ds$Cell.type == i & ds$Visit..Week == j)}
      v[[k]]=var(ds$Median.internalization.scor[ix], na.rm=TRUE)
      #print(v[k])
    }
    plot(q,unlist(v), col=2, main=i, xlab="Cell Count", ylab="var(MIS)", ylim=c(0, 0.6))
    lines(q,unlist(v), col=3)
  }
  title(j, outer=TRUE)
}




###################################################
#   CORRELATION PLOTS1    ##
###################################################
p=c("Median.internalization.score", "MAD", "Rd.score", "Cell.count", "X..total.cells", "relativePASI")
ds.t<-dataset[,-c(1,2,3,4,5,6, 16,17,18,20,21,22)]
#ds.t<-dataset
#chart.Correlation(ds, histogram=TRUE, pch=19)
rquery.cormat(ds.t)

###########################################################################
# DISTRUBUTION OF RD AND MIS VALUES USING DIFFERENT CELL COUNT THRESHOLDS #
###########################################################################
plot(ds$Median.internalization.score,   ds$Rd.score)
##
THRES=c(0, 10, 40 , 60, 100, 200)
l=list(MIS=ds$Median.internalization.score, RD_SCORE=ds$Rd.score)
for(j in c(1:length(l)))
{
  dev.new()
  par(mfcol=c(2,3))
  for(i in THRES)
  {
    d=l[[j]]
    ix=which(ds$Cell.count>i)
    hist(d[ix], n=100, ylim=c(0,1500), xlim=c(-5,5), main=paste("cell count > ", i), xlab=names(l)[j])
    lines(density(d[ix], na.rm=TRUE))             # add a density estimate with defaults
    lines(density(d[ix], adjust=2,na.rm=TRUE), lty="dotted")   # add another "smoother" density
  }
}

plot(ds$Median.internalization.score,   ds$Rd.score)
##
THRES=c(10, 40 , 60, 100, 200)

l=list(MIS=ds$Median.internalization.score, RD_SCORE=ds$Rd.score)
for(j in c(1:length(l)))
{
  dev.new()
  par(mfcol=c(2,3))
  for(i in THRES)
  {
    d=l[[j]]
    ix=which(ds$Cell.count<i )
    hist(d[ix], n=100, ylim=c(0,1500), xlim=c(-5,5), main=paste("cell count < ", i), xlab=names(l)[j])
    lines(density(d[ix], na.rm=TRUE))             # add a density estimate with defaults
    lines(density(d[ix], adjust=2,na.rm=TRUE), lty="dotted")   # add another "smoother" density
  }
}



#################################################################################################
# DISTRIBUTION OF RD AND MIS VALUES USING DIFFERENT CELL COUNT THRESHOLDS, EACH CELL TYPE AT A TIME #
#################################################################################################

THRES=c(0, 10, 40 , 60, 100, 200)
l=list(MIS=ds$Median.internalization.score, RD_SCORE=ds$Rd.score)
l2=list()
c=1
for(j in c(1:length(l)))
{
  dev.new()
  par(mfcol=c(5,3))
  for(i in THRES)
  {
    for(kk in levels(ds$Cell.type))
    {
      d=l[[j]]
      ix=which(ds$Cell.count>i & ds$Cell.type==kk)
      hist(d[ix], n=100, ylim=c(0,1500), main=paste("cell count > ", i), xlab=names(l)[j])
      a= shapiro.test(d[ix])
      l2[[c]]= a$p.value
      lines(density(d[ix], na.rm=TRUE))             # add a density estimate with defaults
      lines(density(d[ix], adjust=2,na.rm=TRUE), lty="dotted")   # add another "smoother" density
    }
    c=c+1
  }
}
dev.new()
hist(-log10(unlist(l2)), n=100, col="blue", xlab="-log10 p-values", ylab="frequencies", main="DISTRIBUTION OF SHAPIRO TESTS P-VALUES")

####################################################################
# DISTRIBUTION OF relative PASI VALUES and Rd.score for each block #
###################################################################
hist(ds$relativePASI, n=20)
lines(density(ds$relativePASI, na.rm=TRUE))             # add a density estimate with defaults


ds$Block5 <- factor(paste(ds$Cell.type, ds$Stimulation,sep="__"))    
histogram(~ relativePASI | Block4, data = ds, as.table = TRUE) 
histogram(~ relativePASI | Block5, data = ds, as.table = TRUE) 
histogram(~ relativePASI.log10 | Block5, data = ds, as.table = TRUE) 

histogram(~ Rd.score | Block4, data = ds, as.table = TRUE) 

shapiro.test(ds$relativePASI.log10[1:5000]) 
shapiro.test(ds$Rd.score[1:5000]) 

l=list()
for(i in levels(ds$Block4))
{
  a= shapiro.test(ds$relativePASI[ds$Block4==i])
  l[[i]]= a$p.value
}
hist(unlist(l), n=20, col="blue", xlab="p-values", ylab="frequencies")

l=list()
for(i in levels(ds$Block4))
{
  a= try(shapiro.test(ds$Rd.score[ds$Block4==i]))
  if(class(a)!="try-error")
  {
    l[[i]]= a$p.value
  }else{l[[i]]=NA}
}
hist(unlist(l), n=20, col="blue", xlab="p-values", ylab="frequencies")

##To check which distribution is the closest visually:
##descdist(expression, discrete = FALSE) from package fitdistrplus 

##To fit distributions and check qqplots:

##fit.gamma<-fitdist(expression, "gamma") from package MASS
##plot(fit.gamma)

##fit.norm<-fitdist(expression, "norm") from package MASS
##plot(fit.norm)



##############################################
##  OVERALL ANOVA AND CHECK OF RESIDUALS    ##
##############################################
t1<-lm(relativePASI.log10 ~  1 + Gender +  Ethnicity + Age + Rd.score, data = ds, na.action="na.omit")
t2<-lm(relativePASI ~  1 + Gender +  Ethnicity + Age + Rd.score, data = ds, na.action="na.omit")
t3<-lm(Rd.score ~  1 + Gender +  Ethnicity + Age + relativePASI.log10, data = ds, na.action="na.omit")
t4<-lm(Rd.score ~  1 + Gender +  Ethnicity + Age + relativePASI, data = ds, na.action="na.omit")
t5<-lm(Rd.score ~  1 + Gender +  Ethnicity + Age + Treatment*tp, random = ~ 1|Patient.ID, data = ds, na.action="na.omit")

par(mfcol=c(2,4))
plot(t1,1:2)
plot(t2,1:2)
plot(t3,1:2)
plot(t4,1:2)

Anova(t1)
Anova(t2)
Anova(t3)
Anova(t4)

## TESTS HOMOSCEDASTICITY   ####
ncvTest(t1)
ncvTest(t2)
ncvTest(t3)
ncvTest(t4)
ncvTest(t5)

##################################
##################################


######
#######

######################################################################################################################################################
### Q2: ARE TREATMENTS HAVING MOLECULAR EFFECTS (differences across time points) ? IF SO, IN WHICH TF/CELL TYPE/STIMULI/TIME POINTS/TREATMENTS ?  ####
### Q3: ARE TREATMENTS HAVING DIFFERENT EFFECTS ? IF SO, IN WHICH TF/CELL TYPE/STIMULI                                                             ####
######################################################################################################################################################
CC_THRES=40
ix=which(ds.m.filter$Cell.count > CC_THRES)
ds.tmp=ds.m.filter[ix,]

results.q2<-list()
homos.q2.p.val<-list()
for(cell.type in  levels(ds.tmp$Cell.type))
{
  print(cell.type)
  for(stimulation in levels(ds.tmp$Stimulation))
  {
    for(ts in levels(ds.tmp$Transcription.factor))
    {  
      ##TEST
      #cell.type="NKT cells"
      #stimulation="IFN-gamma"
      #ts="STAT4"
      ##
    
      block<-paste(cell.type, stimulation, ts, sep="__")
      print(block)
      ds.t<-ds.tmp[ds.tmp$Cell.type==cell.type & ds.tmp$Stimulation==stimulation & ds.tmp$Transcription.factor==ts,]
      ds.t$TreatTime <- droplevels(interaction(ds.t$Treatment, ds.t$tp))          #
      ds.t$Ethnicity<-droplevels(ds.t)$Ethnicity
      ds.t$Gender<-droplevels(ds.t)$Gender
      ds.t$Age<-(ds.t$Age-mean(ds.t$Age))/sd(ds.t$Age)           ##z-scores the ages
      if(stimulation == "Unstimulated")
      {
        ds.t$Rd.score<-ds.t$Median.internalization.score   ## for unstimulated cells it uses MIS instead of Rd.Score in the model
        print(length(ds.t$Rd.score))
      }
      if ((max(ds.t$Rd.score, na.rm=TRUE)-min(ds.t$Rd.score, na.rm=TRUE) > 0.3 & (dim(ds.t)[1] > 5)))
      #if(dim(ds.t)[1]>0)
      {
        e1=try(mod1 <- lme(Rd.score ~  1 + Gender +  Ethnicity + Age + Treatment + tp,  random = ~ 1|Patient.ID, data = ds.t,     na.action="na.omit"))   
        e2=try(mod2 <- lme(Rd.score ~  1 + Gender +  Ethnicity + Age + TreatTime,     random = ~ 1|Patient.ID, data=ds.t,       na.action="na.omit"))
        e3=try(mod3 <- lm(Rd.score ~  1 + Gender +  Ethnicity + Age + TreatTime, data=ds.t,       na.action="na.omit"))   
        homos.q2.p.val=try(ncvTest(mod3)$p)
        if(class(e1)!="try-error")
        {
          results.q2[[block]]$time12         <- try(summary(glht(mod1, linfct=mcp(tp=c("two - one = 0")))))
          results.q2[[block]]$time13         <- try(summary(glht(mod1, linfct=mcp(tp=c("three - one = 0")))))
          results.q2[[block]]$time14         <- try(summary(glht(mod1, linfct=mcp(tp=c("four - one = 0"))))) 
          results.q2[[block]]$treatment      <- try(summary(glht(mod1, linfct=mcp(Treatment=c("Ustekinumab - Adalimumab == 0")))))
          rm(mod1)
        }
        else{
          results.q2[[block]]$time12         <- NA
          results.q2[[block]]$time13         <- NA
          results.q2[[block]]$time14         <- NA
          results.q2[[block]]$treatment      <- NA
        }
        if(class(e2)!="try-error")
        {
          results.q2[[block]]$treatA12       <- try(summary(glht(mod2, linfct=mcp(TreatTime = c("Adalimumab.two - Adalimumab.one = 0"))))) 
          results.q2[[block]]$treatA13       <- try(summary(glht(mod2, linfct=mcp(TreatTime = c("Adalimumab.three - Adalimumab.one = 0"))))) 
          results.q2[[block]]$treatA14       <- try(summary(glht(mod2, linfct=mcp(TreatTime = c("Adalimumab.four - Adalimumab.one = 0"))))) 
          results.q2[[block]]$treatU12       <- try(summary(glht(mod2, linfct=mcp(TreatTime = c("Ustekinumab.two - Ustekinumab.one = 0"))))) 
          results.q2[[block]]$treatU13       <- try(summary(glht(mod2, linfct=mcp(TreatTime = c("Ustekinumab.three - Ustekinumab.one = 0"))))) 
          results.q2[[block]]$treatU14       <- try(summary(glht(mod2, linfct=mcp(TreatTime = c("Ustekinumab.four - Ustekinumab.one = 0")))))
          rm(mod2)
        }else{
          results.q2[[block]]$treatA12       <- NA
          results.q2[[block]]$treatA13       <- NA
          results.q2[[block]]$treatA14       <- NA
          results.q2[[block]]$treatU12       <- NA 
          results.q2[[block]]$treatU13       <- NA 
          results.q2[[block]]$treatU14       <- NA
        } 
      } 
    }    
  }
}



f2 <- function(results, field, n_of_tests) 
{  
  p<-data.frame(p_value=numeric(), avg_ratio=numeric())
  c=1
  for(i in names(results.q2))
  {
    #print(i)
    a<-results.q2[[i]][[field]]
    if(max(class(a) != "try-error") & !is.na(a))
    {
      p[c,"p_value"]<-a$test$pvalues
      p[c,"avg_ratio"]<-a$test$coefficients
    }
    else
    {
      p[c,"p_value"]  <- NA
      p[c,"avg_ratio"]<- NA
    }
    c=c+1
  }
  p[,"q_value"]<-p.adjust(p[,"p_value"], method = "fdr", n = length(p[,"p_value"])*n_of_tests)
  return(p)
}

p<-list()
p[[1]]<-f2(results.q2, "time12",1)
p[[2]]<-f2(results.q2, "time13",1)
p[[3]]<-f2(results.q2, "time14",1)
p[[4]]<-f2(results.q2, "treatment",1)
p[[5]]<-f2(results.q2, "treatA12",1)
p[[6]]<-f2(results.q2, "treatA13",1)
p[[7]]<-f2(results.q2, "treatA14",1)
p[[8]]<-f2(results.q2, "treatU12",1)
p[[9]]<-f2(results.q2, "treatU13",1)
p[[10]]<-f2(results.q2, "treatU14",1)
names(p) = c("time12","time13","time14","treatment","treatA12","treatA13","treatA14","treatU12","treatU13","treatU14")

#############################################################
##  PRINTS ALL Q2/Q3 RESULTS IN A TABLE                    ##  
#############################################################
Results.Q2Q3=data.frame(experiment_name=names(results.q2), rep(CC_THRES, length(results.q2)), 
      Time1_2.pval=p$time12[,1],Time1_2.qval=p$time12[,3], Time1_2.avg.ratio=p$time12[,2], 
      Time1_3.pval=p$time13[,1],Time1_3.qval=p$time13[,3],Time1_3.avg.ratio=p$time13[,2], 
      Time1_4.pval=p$time14[,1],Time1_4.qval=p$time14[,3],Time1_4.avg.ratio=p$time14[,2],
      Treatment.pval=p$treatment[,1],Treatment.qval=p$treatment[,3], Treatment.avg.ratio=p$treatment[,2], 
      TreatmentA1_A2.pval=p$treatA12[,1],TreatmentA1_A2.qval=p$treatA12[,3],TreatmentA1_A2.avg.ratio=p$treatA12[,2], 
      TreatmentA1_A3.pval=p$treatA13[,1],TreatmentA1_A3.qval=p$treatA13[,3],TreatmentA1_A3.avg.ratio=p$treatA13[,2], 
      TreatmentA1_A4.pval=p$treatA14[,1],TreatmentA1_A4.qval=p$treatA14[,3],TreatmentA1_A4.avg.ratio=p$treatA14[,2], 
      TreatmentU1_U2.pval=p$treatU12[,1],TreatmentU1_U2.qval=p$treatU12[,3],TreatmentU1_U2.avg.ratio=p$treatU12[,2], 
      TreatmentU1_U3.pval=p$treatU13[,1],TreatmentU1_U3.qval=p$treatU13[,3],TreatmentU1_U3.avg.ratio=p$treatU13[,2],
      TreatmentU1_U4.pval=p$treatU14[,1],TreatmentU1_U4.qval=p$treatU14[,3],TreatmentU1_U4.avg.ratio=p$treatU14[,2],
      homos.q2.p.val=unlist(homos.q2.p.val))
write.table(Results.Q2Q3, file=paste(out_dir,"Q2Q3_",CC_THRES,".csv", sep=""), col.names = TRUE, sep=",", row.names = F)

################################################
##  QQ PLOT AND HISTOGRAMS                    ##
################################################
par(mfrow = c(5,2))
for(a in 1:length(p))
{
  print(a)
  qqplot(p[[a]][,1],runif(100000), main=names(p)[a], xlab="observed p-values", ylab="expected p-values")
  qqline(p[[a]][,1],distribution=qunif)
}

par(mfrow = c(5,2))
for(a in 1:length(p))
{
  print(a)
  hist(p[[a]][,1], n=80, main=names(p)[a], xlab="observed p-values")
}
## checks only results from unstimulated cells
par(mfrow = c(5,2))
for(a in 1:length(p))
{
  print(a)
  rownames(p[[a]])=names((results.q2))
  ix<-grep("Unstimulated",rownames(p[[a]]))
  hist(p[[a]][ix,1], n=20, main=names(p)[a], xlab="observed p-values")
}



###################################################################################################################################################################################
## Q1: Correlation test of all cell types, stimulated with different cytokines, at different time points, between Rd Score for a given transcription factor and relative PASI    ##
###################################################################################################################################################################################
CC_THRES=40
ix=which(ds.m.filter$Cell.count > CC_THRES)
ds.tmp=ds.m.filter[ix,]

results<-list()
homos.q1.relativePASI.p.val<-list()
homos.q1.PASI90.p.val<-list()
homos.q1.PASI75.p.val<-list()
        

for(cell.type in  levels(ds.tmp$Cell.type))
{
  for(stimulation in levels(ds.tmp$Stimulation))
  {
    for(ts in levels(ds.tmp$Transcription.factor))
    {
      for(tr in levels(ds.tmp$Treatment))
      {
        for(time.point in levels(ds.tmp$Visit..Week))
        {
        
        ### TEST    B cells__LPS__NF-kB__Adalimumab__Baseline
        #cell.type <- "B cells"
        #tr <-  "Adalimumab"
        #time.point <- "Baseline"
        #ts <- "NF-kB"
        #stimulation <- "LPS"
        ####
        block<-paste(cell.type, stimulation, ts,tr, time.point, sep="__")
        print(block)
        ds.t<-ds.tmp[ds.tmp$Cell.type==cell.type & ds.tmp$Stimulation==stimulation & ds.tmp$Transcription.factor==ts & ds.tmp$Treatment==tr & ds.tmp$Visit..Week==time.point,]
        if(stimulation == "Unstimulated")
        {
          ds.t$Rd.score<-ds.t$Median.internalization.score   ## for unstimulated cells it uses MIS instead of Rd.Score in the model
          print(length(ds.t$Rd.score))
        }
        if ((max(ds.t$Rd.score, na.rm=TRUE)-min(ds.t$Rd.score, na.rm=TRUE) > 0.3 & (dim(ds.t)[1] > 5))) ## checks that there are at least 5 values in the module and that there is a max-min delta of at least 0.3
        {
          ix<-which(!is.na(ds.t$ExtremeResponse))
          tmp<-ds.t[ix,]
          
          ## UNIVARIATE ANALYSIS
          
          results[[block]]$lm.relativePASI    <- try(summary(lm(relativePASI ~    1 + Rd.score, data = ds.t, na.action="na.omit")))    
          results[[block]]$lm.PASI90          <- try(summary(lm(PASI90 ~          1 + Rd.score, data = ds.t, na.action="na.omit")))  
          results[[block]]$lm.PASI75          <- try(summary(lm(PASI75 ~          1 + Rd.score, data = ds.t, na.action="na.omit")))
          mod1<-lm(relativePASI ~    1 + Rd.score, data = ds.t, na.action="na.omit")
          mod2<-lm(PASI90 ~          1 + Rd.score, data = ds.t, na.action="na.omit")
          mod3<-lm(PASI75 ~          1 + Rd.score, data = ds.t, na.action="na.omit")
          
          homos.q1.relativePASI.p.val[[block]]        <- try(ncvTest(mod1)$p)          ## CHECKS HOMOSCEDASTICITY FOR EACH OF THE THREE MODELS
          homos.q1.PASI90.p.val[[block]]              <- try(ncvTest(mod2)$p)
          homos.q1.PASI75.p.val[[block]]              <- try(ncvTest(mod3)$p)
        
           ## MULTIVARIATE ANALYSIS
          
          if(length(unique(tmp$ExtremeResponse))>1 & length(tmp$ExtremeResponse) > 5 )
          {
             results[[block]]$lm.extreme         <- try(summary(lm(ExtremeResponse ~   1 + Rd.score, data = tmp, na.action="na.omit")))  
          } else results[[block]]$lm.extreme=NA
          e1<-length(unique(ds.t$Ethnicity[!is.na(ds.t$Ethnicity)]))
          e2<-length(unique(ds.t$Gender[!is.na(ds.t$Ethnicity)]))
          
          if(e1 > 1 & e2 >1)
          {
            results[[block]]$lme.relativePASI   <- try(summary(lme(relativePASI ~       1 + Gender +  Ethnicity + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            results[[block]]$lme.PASI90         <- try(summary(lme(PASI90 ~             1 + Gender +  Ethnicity + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            results[[block]]$lme.PASI75         <- try(summary(lme(PASI75 ~             1 + Gender +  Ethnicity + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            if(length(unique(tmp$ExtremeResponse))>1 & length(tmp$ExtremeResponse) > 5 )
            {
              results[[block]]$lme.extreme        <- try(summary(lme(ExtremeResponse ~    1 + Gender +  Ethnicity + Age + Weight.baseline..Kg. + Rd.score, random = ~ 1|Patient.ID, data = tmp, na.action="na.omit")))
            } else results[[block]]$lme.extreme=NA
          }
          else if(e1 > 1 & e2 <=1)
          {
           results[[block]]$lme.relativePASI   <- try(summary(lme(relativePASI ~       1  +  Ethnicity + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
           results[[block]]$lme.PASI90         <- try(summary(lme(PASI90 ~             1  +  Ethnicity + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
           results[[block]]$lme.PASI75         <- try(summary(lme(PASI75 ~             1  +  Ethnicity + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
           if(length(unique(tmp$ExtremeResponse))>1 & length(tmp$ExtremeResponse) > 5 )
           {
            results[[block]]$lme.extreme        <- try(summary(lme(ExtremeResponse ~    1  +  Ethnicity + Age + Weight.baseline..Kg. + Rd.score, random = ~ 1|Patient.ID, data = tmp, na.action="na.omit")))
           } else results[[block]]$lme.extreme=NA
          }
          else if(e1<=1 & e2 > 1)
          {
            results[[block]]$lme.relativePASI   <- try(summary(lme(relativePASI ~       1 + Gender  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            results[[block]]$lme.PASI90         <- try(summary(lme(PASI90 ~             1 + Gender  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            results[[block]]$lme.PASI75         <- try(summary(lme(PASI75 ~             1 + Gender  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            if(length(unique(tmp$ExtremeResponse))>1 & length(tmp$ExtremeResponse) > 5 )
            {
              results[[block]]$lme.extreme        <- try(summary(lme(ExtremeResponse ~    1 + Gender  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = tmp, na.action="na.omit")))
            } else results[[block]]$lme.extreme=NA
          }
          else if(e1<=1 & e2 <=1)
          {
            results[[block]]$lme.relativePASI   <- try(summary(lme(relativePASI ~       1  + Age + Weight.baseline..Kg. +  PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            results[[block]]$lme.PASI90         <- try(summary(lme(PASI90 ~             1  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            results[[block]]$lme.PASI75         <- try(summary(lme(PASI75 ~             1  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")))
            if(length(unique(tmp$ExtremeResponse))>1 & length(tmp$ExtremeResponse) > 5)
            {
              results[[block]]$lme.extreme        <- try(summary(lme(ExtremeResponse ~    1  + Age + Weight.baseline..Kg. + PsA + Rd.score, random = ~ 1|Patient.ID, data = tmp, na.action="na.omit")))
            } else results[[block]]$lme.extreme=NA
          } 
        }
      }
    }
  }
}
}

###
##  
##
f3 <- function(results, field, model,n_of_tests) 
{
  p<-data.frame(p_value=numeric(), coef=numeric())
  c=1
  if(model=="lm")
  {
    for(i in names(results))
    {
      print(i)
      a<-results[[i]][[field]]
      if(!is.na(a) & max(class(a) != "try-error"))
      {
        p[c,"p_value"]<-a$coefficients["Rd.score","Pr(>|t|)"]
        p[c,"coef"]<-a$coefficients["Rd.score", "Estimate"]        
      }
      else
      {
        p[c,"p_value"]<-(NA)
        p[c,"coef"]<-(NA)
      }
      c=c+1
    }
  }
  else if(model=="lme")
  {
   for(i in names(results))
   {
    a<-results[[i]][[field]]
    if(!is.na(a) & max(class(a) != "try-error"))
    {
        p[c,"p_value"]<-a$tTable["Rd.score", "p-value"]
        p[c,"coef"]<-a$tTable["Rd.score", "Value"]    
    }
    else
    {
      p[c,"p_value"]<-NA
      p[c,"coef"]<-(NA)
    }
    c=c+1
   }
  }
  p[,"q_value"]<-p.adjust(p[,"p_value"], method = "fdr", n = length(p[,"p_value"]))
  return(p)
}

p<-list()
p[[1]]<-f3(results, "lm.relativePASI","lm", l)
p[[2]]<-f3(results, "lm.PASI90","lm", l)
p[[3]]<-f3(results, "lm.PASI75","lm", l)
p[[4]]<-f3(results, "lm.extreme","lm", l)
p[[5]]<-f3(results, "lme.relativePASI","lme", l)
p[[6]]<-f3(results, "lme.PASI90","lme", l)
p[[7]]<-f3(results, "lme.PASI75","lme", l)
p[[8]]<-f3(results, "lme.extreme","lme", l)

names(p) = c("lm.relativePASI","lm.PASI90","lm.PASI75","lm.extreme","lme.relativePASI","lme.PASI90","lme.PASI75","lme.extreme")

#############################################################
##  PRINTS ALL Q1 RESULTS IN A TABLE                    ##
#############################################################
Results.Q1=data.frame(experiment_name=names(results), rep(CC_THRES, length(results)), 
      lm.relativePASI.pval=p$lm.relativePASI[,1], lm.relativePASI.qval=p$lm.relativePASI[,3],  lm.relativePASI.coef=p$lm.relativePASI[,2], 
      lm.PASI90.pval=p$lm.PASI90[,1],lm.PASI90.qval=p$lm.PASI90[,3],               lm.PASI90.coef=p$lm.PASI90[,2], 
      lm.PASI75.pval=p$lm.PASI75[,1],lm.PASI75.qval=p$lm.PASI75[,3],               lm.PASI75.coef=p$lm.PASI75[,2],
      lm.extreme.pval=p$lm.extreme[,1],lm.extreme.qval=p$lm.extreme[,3],             lm.extreme.coef=p$lm.extreme[,2], 
      lme.relativePASI.pval=p$lme.relativePASI[,1], lme.relativePASI.qval=p$lme.relativePASI[,3],lme.relativePASI.coef=p$lme.relativePASI[,2], 
      lme.PASI90.pval=p$lme.PASI90[,1], lme.PASI90.qval=p$lme.PASI90[,3],            lme.PASI90.coef=p$lme.PASI90[,2], 
      lme.PASI75.pval=p$lme.PASI75[,1],lme.PASI75.qval=p$lme.PASI75[,3],             lme.PASI75.coef=p$lme.PASI75[,2], 
      lme.extreme.pval=p$lme.extreme[,1],lme.extreme.qval=p$lme.extreme[,3],           lme.extreme.coef=p$lme.extreme[,2],
      homos.q1.relativePASI.p.val=unlist(homos.q1.relativePASI.p.val), homos.q1.PASI90.p.val=unlist(homos.q1.PASI90.p.val),homos.q1.PASI75.p.val=unlist(homos.q1.PASI75.p.val))              
write.table(Results.Q1, file=paste(out_dir,"Q1_",CC_THRES,".csv", sep=""), col.names = TRUE, sep=",", row.names = F)


################################################
##  QQ PLOT AND HISTOGRAMS                    ##
################################################
par(mfrow = c(4,2))
for(a in 1:length(p))
{
  print(a)
  qqplot(p[[a]][,1],runif(100000), main=names(p)[a], xlab="observed p-values", ylab="expected p-values")
  qqline(p[[a]][,1],distribution=qunif)
}

par(mfrow = c(5,2))
for(a in 1:10)
{
  print(a)
  hist(p[[a]][,1], n=80, main=names(p)[a], xlab="observed p-values")
}

###############################################################
##             DRUG CONCENTRATIONS AND RD SCORE              ##
###############################################################
boxplot(dataset$Rd.score~dataset$USK.conc, col="green", xlab="USK drug conc.", ylab="Rd Score", main="TF Internalization and USK Drug Concentration: General Trend")
dev.new()
boxplot(dataset$Rd.score~dataset$ADL.conc, col="green", xlab="USK drug conc.", ylab="Rd Score", main="TF Internalization and ADL Drug Concentration: General Trend")
###
###
###


##############################
###############################
##############################
p.lme.relativePASI<-list()
p.lm.relativePASI<-list()
pv.lme.relativePASI<-array()
pv.lm.relativePASI<-array()

p.lme.PASI90<-list()
p.lm.PASI90<-list()
pv.lme.PASI90<-array()
pv.lm.PASI90<-array()

p.lme.PASI75<-list()
p.lm.PASI75<-list()
pv.lme.PASI75<-array()
pv.lm.PASI75<-array()

p.lme.extreme<-list()
p.lm.extreme<-list()
pv.lme.extreme<-array()
pv.lm.extreme<-array()


c1<-1
c2<-1
c3<-1
c4<-1
c5<-1
c6<-1
c7<-1
c8<-1

p.thres<-0.05

for(i in names(results))
{
 r.lme.relativePASI=results[[i]]$lme.relativePASI 
 r.lm.relativePASI=results[[i]]$lm.relativePASI

 r.lme.PASI90=results[[i]]$lme.PASI90 
 r.lm.PASI90=results[[i]]$lm.PASI90  

 r.lme.PASI75=results[[i]]$lme.PASI75 
 r.lm.PASI75=results[[i]]$lm.PASI75  

 r.lme.extreme=results[[i]]$lme.extreme 
 r.lm.extreme=results[[i]]$lm.extreme  

 
 if(class(r.lme.relativePASI) != "try-error")
 {
  pv.lme.relativePASI[c1]<-r.lme.relativePASI$tTable["Rd.score","p-value"] 
  if(pv.lme.relativePASI[c1]<p.thres){p.lme.relativePASI[[i]]$pval<-pv.lme.relativePASI[c1]}
  c1<-c1+1
 }
 if(class(r.lm.relativePASI) != "try-error")
 {
  pv.lm.relativePASI[c2]<-r.lm.relativePASI$coefficients["Rd.score",4]
  if(pv.lm.relativePASI[c2]<p.thres){p.lm.relativePASI[[i]]$pval<-pv.lm.relativePASI[c2]}
  c2<-c2+1
 }
  if(class(r.lme.PASI90) != "try-error")
 {
  pv.lme.PASI90[c3]<-r.lme.PASI90$tTable["Rd.score","p-value"] 
  if(pv.lme.PASI90[c3]<p.thres){p.lme.PASI90[[i]]$pval<-pv.lme.PASI90[c3]}
  c3<-c3+1
 }
 if(class(r.lm.PASI90) != "try-error")
 {
  pv.lm.PASI90[c4]<-r.lm.PASI90$coefficients["Rd.score",4]
  if(pv.lm.PASI90[c4]<p.thres){p.lm.PASI90[[i]]$pval<-pv.lm.PASI90[c4]}
  c4<-c4+1
 }
  if(class(r.lme.PASI75) != "try-error")
 {
  pv.lme.PASI75[c5]<-r.lme.PASI75$tTable["Rd.score","p-value"] 
  if(pv.lme.PASI75[c5]<p.thres){p.lme.PASI75[[i]]$pval<-pv.lme.PASI75[c5]}
  c5<-c5+1
 }
 if(class(r.lm.PASI75) != "try-error")
 {
  pv.lm.PASI75[c6]<-r.lm.PASI75$coefficients["Rd.score",4]
  if(!is.na(pv.lm.PASI75[c6]) & pv.lm.PASI75[c6]<p.thres){p.lm.PASI75[[i]]$pval<-pv.lm.PASI75[c6]}
  c6<-c6+1
 }
 if(class(r.lme.extreme) != "try-error" & !is.na(r.lme.extreme))
 {
  pv.lme.extreme[c7]<-r.lme.extreme$tTable["Rd.score","p-value"] 
  if(pv.lme.extreme[c7]<p.thres){p.lme.extreme[[i]]$pval<-pv.lme.extreme[c7]}
  c7<-c7+1
 }
 if(class(r.lm.extreme) != "try-error"  & !is.na(r.lm.extreme))
 {
  pv.lm.extreme[c8]<-r.lm.extreme$coefficients["Rd.score",4]
  if(!is.na(pv.lm.extreme[c8]) & pv.lm.extreme[c8]<p.thres){p.lm.extreme[[i]]$pval<-pv.lm.extreme[c8]}
  c8<-c8+1
 }
} 


c=which(unlist(p.lm.relativePASI) == min(unlist(p.lm.relativePASI)))
p.lm.relativePASI[c]

c=which(unlist(p.lme.relativePASI) == min(unlist(p.lme.relativePASI)))
p.lme.relativePASI[c]

#########################################################
## PLOTS RESULTS PASSING THE THRESHOLD      ##
#########################################################
output.dir <- "C:/Users/derinald/Desktop/ImageStream Data Set - October 2016/R/output/"

if(length(p.lme.relativePASI)>0)
{
  #sink(paste(output.dir, "lme.relativePASI.p.values", ".txt", sep=""), append=FALSE)
  print(p.lme.relativePASI)
  #sink()
  #pdf(paste(output.dir, "lme.relativePASI", ".pdf", sep=""),width=7,height=5)
  par(mfrow = c(ceiling(length(p.lme.relativePASI)/3),3))
  for(i in names(p.lme.relativePASI))
  {
    ix<-which(ds$Block4==i)
    plot(ds$Rd.score[ix]~ds$relativePASI[ix], pch = 16, cex = 1.3, col = "blue", main=gsub("\\?", "K", i))
    abline(lm(ds$Rd.score[ix]~ds$relativePASI[ix]))
  }
  #dev.off()
}

if(length(p.lm.relativePASI)>0)
{
  #sink(paste(output.dir, "lm.relativePASI.p.values", ".txt", sep=""), append=FALSE)
  print(p.lm.relativePASI)
  #sink()
  #pdf(paste(output.dir, "lm.relativePASI", ".pdf", sep=""),width=7,height=10)
  par(mfrow = c(ceiling(length(p.lm.relativePASI)/3),3))
  for(i in names(p.lm.relativePASI))
  {
    ix<-which(ds$Block4==i)
    plot(ds$Rd.score[ix]~ds$relativePASI[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i))
    abline(lm(ds$Rd.score[ix]~ds$relativePASI[ix]))
  }
  #dev.off()
}

if(length(p.lme.PASI90)>0)
{
  #sink(paste(output.dir, "p.lme.PASI90.p.values", ".txt", sep=""), append=FALSE)
  print(p.lme.PASI90)
  #sink()
  #pdf(paste(output.dir, "p.lme.PASI90", ".pdf", sep=""),width=7,height=5)
  par(mfrow = c(ceiling(length(p.lme.PASI90)/3),3))
  for(i in names(p.lme.PASI90))
  {
    ix<-which(ds$Block4==i)
    boxplot(ds$Rd.score[ix]~ds$PASI90[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i), xlab="PASI90")
    stripchart(ds$Rd.score[ix]~ds$PASI90[ix], , vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'red',cex=1.6)
  }
  #dev.off()
}

if(length(p.lm.PASI90)>0)
{
  #sink(paste(output.dir, "p.lm.PASI90.p.values", ".txt", sep=""), append=FALSE)
  print(p.lm.PASI90)
  #sink()
  #pdf(paste(output.dir, "p.lm.PASI90", ".pdf", sep=""),width=7,height=5)
  par(mfrow = c(ceiling(length(p.lm.PASI90)/3),3))
  for(i in names(p.lm.PASI90))
  {
    ix<-which(ds$Block4==i)
    boxplot(ds$Rd.score[ix]~ds$PASI90[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i), xlab="PASI90")
    stripchart(ds$Rd.score[ix]~ds$PASI90[ix], , vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'red',cex=1.6)
  }
  #dev.off()
}

if(length(p.lme.PASI75)>0)
{
  #sink(paste(output.dir, "p.lme.PASI75.p.values", ".txt", sep=""), append=FALSE)
  print(p.lme.PASI75)
  #sink()
  #pdf(paste(output.dir, "p.lme.PASI75", ".pdf", sep=""),width=20,height=10)
  par(mfrow = c(ceiling(length(p.lme.PASI75)/3),3))
  for(i in names(p.lme.PASI75))
  {
    ix<-which(ds$Block4==i)
    boxplot(ds$Rd.score[ix]~ds$PASI75[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i), xlab="PASI75")
    stripchart(ds$Rd.score[ix]~ds$PASI75[ix], , vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'red',cex=1.6)
  }
  #dev.off()
}

if(length(p.lm.PASI75)>0)
{
  #sink(paste(output.dir, "p.lm.PASI75.p.values", ".txt", sep=""), append=FALSE)
  print(p.lm.PASI75)
  #sink()
  #pdf(paste(output.dir, "p.lm.PASI75", ".pdf", sep=""),width=7,height=5)
  par(mfrow = c(ceiling(length(p.lm.PASI75)/3),3))
  for(i in names(p.lm.PASI75))
  {
    ix<-which(ds$Block4==i)
    boxplot(ds$Rd.score[ix]~ds$PASI75[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i), xlab="PASI75")
    stripchart(ds$Rd.score[ix]~ds$PASI75[ix], , vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'red',cex=1.6)
  }
  #dev.off()
}

if(length(p.lme.extreme)>0)
{
  #sink(paste(output.dir, "p.lme.extreme.p.values", ".txt", sep=""), append=FALSE)
  print(p.lme.extreme)
  #sink()
  #pdf(paste(output.dir, "p.lme.extreme", ".pdf", sep=""),width=7,height=5)
  par(mfrow = c(ceiling(length(p.lme.extreme)/3),3))
  for(i in names(p.lme.extreme))
  {
    ix<-which(ds$Block4==i)
    boxplot(ds$Rd.score[ix]~ds$ExtremeResponse[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i), xlab="extreme")
    stripchart(ds$Rd.score[ix]~ds$ExtremeResponse[ix], , vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'red',cex=1.6)
  }
  #dev.off()
}

if(length(p.lm.extreme)>0)
{
  #sink(paste(output.dir, "p.lm.extreme.p.values", ".txt", sep=""), append=FALSE)
  print(p.lm.extreme)
  #sink()
  #pdf(paste(output.dir, "p.lm.extreme", ".pdf", sep=""),width=7,height=5)
  par(mfrow = c(ceiling(length(p.lm.extreme)/3),3))
  for(i in names(p.lm.extreme))
  {
    ix<-which(ds$Block4==i)
    boxplot(ds$Rd.score[ix]~ds$ExtremeResponse[ix], pch = 16, cex = 1.3, col = "blue",main=gsub("\\?", "K", i), xlab="extreme")
    stripchart(ds$Rd.score[ix]~ds$ExtremeResponse[ix], , vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'red',cex=1.6)
  }
  #dev.off()
}
#fig1 = bwplot(ds$Rd.score[ix]~ds$relativePASI[ix]|ds$PASI90[ix], pch = 16)
#print(fig1)



print(p.lm.relativePASI)
print(p.lme.PASI90)
print(p.lm.PASI90)
print(p.lme.PASI75)
print(p.lm.PASI75)
print(p.lme.extreme)
print(p.lm.extreme)

##################################

#################################################################
##  DISTRIBUTION OF P-VALUES                                   ##
#################################################################
par(mfrow = c(4,4))
##qqnorm(p.val.array)
qqplot(pv.lme.relativePASI,runif(100000) )
qqline(pv.lme.relativePASI,distribution=qunif)

hist(pv.lme.relativePASI, n=20)            # prob=TRUE for probabilities not counts
lines(density(pv.lme.relativePASI))             # add a density estimate with defaults
lines(density(pv.lme.relativePASI, adjust=2), lty="dotted")   # add another "smoother" density

qqplot(pv.lm.relativePASI,runif(100000) )
qqline(pv.lm.relativePASI,distribution=qunif)

hist(pv.lm.relativePASI, n=20)            # prob=TRUE for probabilities not counts
lines(density(pv.lm.relativePASI))             # add a density estimate with defaults
lines(density(pv.lm.relativePASI, adjust=2), lty="dotted")   # add another "smoother" density


qqplot(pv.lme.PASI90,runif(100000) )
qqline(pv.lme.PASI90,distribution=qunif)

hist(pv.lme.PASI90, n=20)            # prob=TRUE for probabilities not counts
lines(density(pv.lme.PASI90))             # add a density estimate with defaults
lines(density(pv.lme.PASI90, adjust=2), lty="dotted")   # add another "smoother" density

qqplot(pv.lm.PASI90,runif(100000) )
qqline(pv.lm.PASI90,distribution=qunif)

hist(pv.lm.PASI90, n=20)            # prob=TRUE for probabilities not counts
lines(density(pv.lm.PASI90))             # add a density estimate with defaults
lines(density(pv.lm.PASI90, adjust=2), lty="dotted")   # add another "smoother" density


qqplot(pv.lme.PASI75,runif(100000) )
qqline(pv.lme.PASI75,distribution=qunif)

hist(pv.lme.PASI75, n=30)            # prob=TRUE for probabilities not counts
lines(density(pv.lme.PASI75))             # add a density estimate with defaults
lines(density(pv.lme.PASI75, adjust=2), lty="dotted")   # add another "smoother" density

qqplot(pv.lm.PASI75,runif(100000) )
qqline(pv.lm.PASI75,distribution=qunif)

hist(pv.lm.PASI75, n=20)            # prob=TRUE for probabilities not counts
lines(density(pv.lm.PASI75, na.rm=TRUE))             # add a density estimate with defaults
lines(density(pv.lm.PASI75, adjust=2, na.rm=TRUE), lty="dotted")   # add another "smoother" density


qqplot(pv.lme.extreme,runif(100000) )
qqline(pv.lme.extreme,distribution=qunif)

hist(pv.lme.extreme, n=30)            # prob=TRUE for probabilities not counts
lines(density(pv.lme.extreme))             # add a density estimate with defaults
lines(density(pv.lme.extreme, adjust=2), lty="dotted")   # add another "smoother" density



qqplot(pv.lm.extreme,runif(100000) )
qqline(pv.lm.extreme,distribution=qunif)

hist(pv.lm.extreme, n=20)            # prob=TRUE for probabilities not counts
lines(density(pv.lm.extreme, na.rm=TRUE))             # add a density estimate with defaults
lines(density(pv.lm.extreme, adjust=2, na.rm=TRUE), lty="dotted")   # add another "smoother" density



##################################################################
##################################################################

## OLD STUFF
###########################################
#######################################
################################################
#
####################################################


## FOR EACH BLOCK CHECKS THE DIFFERENCE BETWEEN THE 2 TREATMENTS

p=list()
for (i in unique(ds$Block1))
{
  a <- t.test(ds$Rd.score[ds$Block1==i]~ds$Treatment[ds$Block1==i])
  if(a$p.value<0.001)
  {
    print(i)
    print(a$p.value)
    boxplot(ds$Rd.score[ds$Block1==i]~ds$Treatment[ds$Block1==i])
  }
  p[i]=a$p.value
}
p.a <- unlist(p, recursive = TRUE, use.names = TRUE)

ix=which(p.a<0.0001)

dd<-ds[ds$Block1 %in% ix,]

fig1 = bwplot(dd$Rd.score ~ dd$Treatment | dd$Block1, pch = 16)
print(fig1)

hist(p.a)
qqplot((p.a),runif(100000) )
qqline((p.a),distribution=qunif)


#################################################################################
## FOR EACH BLOCK CHECKS THE CORRELATION BETWEEN PASI90(YES/NO) AND RD.SCORE
##################################################################################

p=list()
for (i in unique(ds$Block12))
{
  a <- t.test(ds$Rd.score[ds$Block12==i]~ds$PASI90[ds$Block12==i])
  if(a$p.value<0.05)
  {
    print(i)
    print(a$p.value)                       
    boxplot(ds$Rd.score[ds$Block12==i]~ds$PASI90[ds$Block12==i])
  }
  p[i]=a$p.value
}
p.a <- unlist(p, recursive = TRUE, use.names = TRUE)

ix=which(p.a<0.05)

dd<-ds[ds$Block12 %in% ix,]

fig2 = bwplot(dd$PASI90 ~ dd$Rd.score | dd$Block12, pch = 16)
print(fig2)


 ######
 
 z4 <- lme(PASI90 ~ Rd.score, random = ~ 1|Block12, data = ds)
b=summary(z4)

#######################################################################################
## FOR EACH BLOCK CHECKS THE CORRELATION BETWEEN relative.PASI(YES/NO) AND RD.SCORE
#######################################################################################

z5 <- lme(relativePASI ~ Rd.score, random = ~ 1|Block12, data = ds)
b=summary(z5)


p=list()
for (i in unique(ds$Block12))
{
  a <- cor.test(ds$Rd.score[ds$Block12==i], ds$relativePASI[ds$Block12==i])
  if(a$p.value<0.05)
  {
    print(i)
    print(a$p.value)
    plot(ds$Rd.score[ds$Block12==i], ds$relativePASI[ds$Block12==i])
    
    abline(lm(ds$relativePASI[ds$Block12==i]~ ds$Rd.score[ds$Block12==i]))
  }
  p[i]=a$p.value
}
p.a <- unlist(p, recursive = TRUE, use.names = TRUE)
hist(p.a)
qqplot((p.a),runif(100000) )
qqline((p.a),distribution=qunif)

names(p.a[p.a<0.05])



dd<-ds[ds$Block12 %in% names(p.a[p.a<0.01]),]
fig3 = xyplot(dd$Rd.score ~ dd$relativePASI | dd$Block12, pch = 16)
print(fig3)



###CROSS _CHECK WITH PREVIOUS RESULTS
check.ix <- which(ds$Block12=="1_NK cells_NF-?B_TNF?_Adalimumab")
plot(ds$relativePASI[check.ix], ds$Rd.score[check.ix], )
abline(lm(ds$Rd.score[check.ix] ~ ds$relativePASI[check.ix]))



















a=ddply(ds,~ds$Block1,summarise,mean=mean(ds$Rd.score),sd=sd(ds$Rd.score))

a=ddply(ds,~Block1,summarise,mean=mean(Rd.score),sd=sd(Rd.score))

for (i in unique(ds$Block1))
{
  a$Treatment[a$Block1==i] <- as.character(ds[ds$Block1==i, "Treatment"][1])
}

a=ddply(ds,~Block1,summarise,mean=mean(Rd.score),sd=sd(Rd.score))

dt <- data.frame(age=rchisq(20,10),group=sample(1:2,20,rep=T))
ddply(dt,~group,summarise,mean=mean(age),sd=sd(age))