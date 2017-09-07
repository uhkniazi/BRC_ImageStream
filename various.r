 ### TEST
        cell.type <- "Monocytes"
        tr <-  "Ustekinumab"
        ts <- "NF-?B"
        stimulation <- "TNF?"

        block<-paste(cell.type, stimulation, ts, sep="__")
        print(block)
        ds.t<-ds[ds$Cell.type==cell.type & ds$Stimulation==stimulation & ds$Transcription.factor==ts,]

        ###TEST
        ds.t<-ds[ds$Cell.type==cell.type & ds$Stimulation==stimulation & ds$Transcription.factor==ts,]
        #t4<-lme(Rd.score ~  1 + Gender +  Ethnicity + Age + Treatment*tp, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit")
        t4<-lm(Rd.score ~  1 + Gender +  Ethnicity + Age + factor(tp)*Treatment, data = ds.t, na.action="na.omit")

        par(mfcol=c(1,2))
        plot(t4,1:2)
        Anova(t4)

        t4.means<-interactionMeans(t4)
        plot(t4.means)

        interactionMeans(t4, factors="Treatment")

        plot(t4.means, multiple=FALSE)
        testInteractions(t4, fixed="Treatment", across="tp")


        #t1=summary(lme(Rd.score ~  1 + Gender +  Ethnicity + Age + Treatment + tp + Treatment*tp, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit"))
        #t2=summary(lme(Rd.score ~  1 + Gender +  Ethnicity + Age + Treatment*tp, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit"))
        #t3=summary(lme(Rd.score ~  1 + Gender +  Ethnicity + Age + Treatment + tp + Treatment:tp, random = ~ 1|Patient.ID, data = ds.t, na.action="na.omit"))



        t4<-lm(relativePASI.log10 ~  1 + Gender +  Ethnicity + Age + Rd.score, data = ds, na.action="na.omit")

        par(mfcol=c(1,2))
        plot(t4,1:2)
        Anova(t4)



        tmp <- expand.grid(tp = factor(unique(ds.t$tp)), Treatment = unique(ds.t$Treatment))

X <- model.matrix(~ Treatment * factor(tp), data = tmp)
glht(mod, linfct = X)
Tukey <- contrMat(table(ds.t$tp), "Tukey")
K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
rownames(K1) <- paste(levels(ds.t$Treatment)[1], rownames(K1), sep = ":")
K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
rownames(K2) <- paste(levels(ds.t$Treatment)[2], rownames(K2), sep = ":")
K <- rbind(K1, K2)
colnames(K) <- c(colnames(Tukey), colnames(Tukey))
summary(glht(mod, linfct = K %*% X))
glht(mod, linfct = mcp(x = "Tukey"))


#z3 <- lme(Rd.score ~ GenderId +  Ethnicity + Age + Block2*Treatment, random = ~ 1|Patient.ID, data = ds, na.action="na.omit")
#b=summary(z3)
#interaction.plot(ds$Treatment, ds$Block2, ds$Rd.score, ds$)
#b <- summary(z3)
#b$tTable
#print(contrast(z3, list(Treatment="Adalimumab", Block2="1_Monocytes_NF-?B_IL-17", GenderId ="1", Ethnicity="3", Age=23),       ##AGE, GENDER AND ETHNICITY ARE JUST DUMMY VALUS AND ARE NOT USED IN THE STATISTICS
#                   list(Treatment="Adalimumab",  Block2="1_Monocytes_NF-?B_TNF?", GenderId ="1", Ethnicity="3",Age=23), X=TRUE))
#
#ix1=ds$Treatment=="Adalimumab" & ds$Block2=="1_Monocytes_NF-?B_IL-17"
#ix2=ds$Treatment=="Adalimumab" & ds$Block2=="1_Monocytes_NF-?B_TNF?"
#boxplot(ds$Rd.score[ix1],ds$Rd.score[ix2], notch=TRUE)
#
#fig1 = xyplot(ds$Rd.score ~ ds$Treatment | ds$Block1, pch = 4)
#print(fig1)