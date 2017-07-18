setwd('C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/mixing_models2')
source <- read.csv('cor_sources.csv')

phytoH <- source$mean_d2HA_lake
periH <- source$mean_d2Hperi_lake
terrH <- source$d2HT
phytoC <- source$mean_d13CA_lake
periC <- source$mean_d13Cperi_lake
terrC <- source$d13CT
phytoN <- source$mean_d15NA_lake
periN <- source$mean_d15Nperi_lake
terrN <- source$d15NT

consltypes <- c(1,2,1,3,3,2,1,3,3,1,4,1,2,1,2,1,1,4,4,1,1,2,3,3,2,2,2,3,3,1,4)
ltypes <- factor(c(2,2,1,1,1,1,3,3,3,2,3))

# summary(mod <- lm(phytoN ~ ltypes))
# summary(mod <- lm(terrH ~ ltypes))

#noob style ####
out <- as.data.frame(matrix(NA, nrow=9, ncol=11))
colnames(out) <- c('source','tracer','tMSML','dfMSML','pMSML','tMSAP','dfMSAP','pMSAP','tMLAP','dfMLAP','pMLAP')
counter <- 1
for(i in c('phyto','peri','terr')){
    for(j in c('C','N','H')){
        out[counter,1] <- i
        out[counter,2] <- j

        x <- eval(parse(text=paste0(i,j)))
        test <- t.test(x[ltypes==1], x[ltypes==2])
        out[counter,3] <- test$statistic
        out[counter,4] <- test$parameter
        out[counter,5] <- test$p.value
        test <- t.test(x[ltypes==1], x[ltypes==3])
        out[counter,6] <- test$statistic
        out[counter,7] <- test$parameter
        out[counter,8] <- test$p.value
        test <- t.test(x[ltypes==2], x[ltypes==3])
        out[counter,9] <- test$statistic
        out[counter,10] <- test$parameter
        out[counter,11] <- test$p.value

        counter <- counter + 1
    }
}
out

#d2h2o
d2 <- read.csv('C:/Users/Mike/git/seven_lakes/Analysis/7lmeans.csv')

# t.test(d2$dD_H2O[3:6], d2$dD_H2O[c(1,2,10)]) #sm vs lg, signif
# t.test(d2$dD_H2O[3:6], d2$dD_H2O[c(7:9,11)]) #sm vs hi, signif
# t.test(d2$dD_H2O[c(1,2,10)], d2$dD_H2O[c(7:9,11)]) #lg vs hi, not signif
lkfac = factor(c(1,1,2,2,2,2,3,3,3,1,3), labels=c('lg','sm','hi'))
aovmod = aov(d2$dD_H2O[1:11] ~ lkfac)
summary(aovmod)
TukeyHSD(aovmod)

#calanoids
cons <- read.csv('cor_consumers.csv')
cal_sm_ind <- c(1,3,7,10)
cal_lg_ind <- c(2,6)
cal_hi_ind <- c(4,5,8,9)

t.test(cons$C_true_ind[cal_sm_ind], cons$C_true_ind[cal_lg_ind])
t.test(cons$C_true_ind[cal_hi_ind], cons$C_true_ind[cal_lg_ind])

t.test(cons$N_true_ind[cal_sm_ind], cons$N_true_ind[cal_lg_ind])
t.test(cons$N_true_ind[cal_hi_ind], cons$N_true_ind[cal_lg_ind])

mean(cons$C_true_ind[cal_lg_ind])
sd(cons$C_true_ind[cal_lg_ind])

mean(cons$N_true_ind[cal_lg_ind])
sd(cons$N_true_ind[cal_lg_ind])



#consumers comparisons (1-way anova) ####
allcons <- read.csv('cor_all.csv')
allcons = allcons[-c(11,18,19,31),]
cal = allcons[1:10,]
clad = allcons[11:16,]
caddis = allcons[17:27,]
# conslkfac = factor(c(1,2,1,3,3,2,1,3,3,1,1,2,1,2,1,1,1,1,2,3,3,2,2,2,3,3,1), 
#                    labels=c('sm','lg','hi'))
# constypfac = factor(allcons$sample_type)

#cal
callkfac = factor(c(1,2,1,3,3,2,1,3,3,1), labels=c('sm','lg','hi'))

#N
consN = aov(allcons$cons_N_u[allcons$sample_type == 'Calanoida'] ~ callkfac)
write.csv(summary(consN)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consNcal.csv')
consNTukeyOut = pairwise.t.test(allcons$cons_N_u[allcons$sample_type == 'Calanoida'], 
                                callkfac, p.adjust.method='bonferroni')
write.csv(consNTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consNpwcal.csv')

#C
consC = aov(allcons$cons_C_u[allcons$sample_type == 'Calanoida'] ~ callkfac)
write.csv(summary(consC)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consCcal.csv')
consCTukeyOut = pairwise.t.test(allcons$cons_C_u[allcons$sample_type == 'Calanoida'], 
                                callkfac, p.adjust.method='bonferroni')
write.csv(consCTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consCpwcal.csv')

#H
consH = aov(allcons$cons_H_u[allcons$sample_type == 'Calanoida'] ~ callkfac)
write.csv(summary(consH)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consHcal.csv')
consHTukeyOut = pairwise.t.test(allcons$cons_H_u[allcons$sample_type == 'Calanoida'], 
                                callkfac, p.adjust.method='bonferroni')
write.csv(consHTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consHpwcal.csv')

#clad
cladlkfac = factor(c(1,2,1,2,1,1), labels=c('sm','lg'))

#N
consN = aov(allcons$cons_N_u[allcons$sample_type == 'Cladocera'] ~ cladlkfac)
write.csv(summary(consN)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consNclad.csv')
consNTukeyOut = pairwise.t.test(allcons$cons_N_u[allcons$sample_type == 'Cladocera'], 
                                cladlkfac, p.adjust.method='bonferroni')
write.csv(consNTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consNpwclad.csv')

#C
consC = aov(allcons$cons_C_u[allcons$sample_type == 'Cladocera'] ~ cladlkfac)
write.csv(summary(consC)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consCclad.csv')
consCTukeyOut = pairwise.t.test(allcons$cons_C_u[allcons$sample_type == 'Cladocera'], 
                                cladlkfac, p.adjust.method='bonferroni')
write.csv(consCTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consCpwclad.csv')

#H
consH = aov(allcons$cons_H_u[allcons$sample_type == 'Cladocera'] ~ cladlkfac)
write.csv(summary(consH)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consHclad.csv')
consHTukeyOut = pairwise.t.test(allcons$cons_H_u[allcons$sample_type == 'Cladocera'], 
                                cladlkfac, p.adjust.method='bonferroni')
write.csv(consHTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consHpwclad.csv')

#caddis
caddislkfac = factor(c(1,1,2,3,3,2,2,2,3,3,1), labels=c('sm','lg','hi'))

#N
consN = aov(allcons$cons_N_u[allcons$sample_type == 'Trichoptera'] ~ caddislkfac)
write.csv(summary(consN)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consNcaddis.csv')
consNTukeyOut = pairwise.t.test(allcons$cons_N_u[allcons$sample_type == 'Trichoptera'], 
                                caddislkfac, p.adjust.method='bonferroni')
write.csv(consNTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consNpwcaddis.csv')

#C
consC = aov(allcons$cons_C_u[allcons$sample_type == 'Trichoptera'] ~ caddislkfac)
write.csv(summary(consC)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consCcaddis.csv')
consCTukeyOut = pairwise.t.test(allcons$cons_C_u[allcons$sample_type == 'Trichoptera'], 
                                caddislkfac, p.adjust.method='bonferroni')
write.csv(consCTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consCpwcaddis.csv')

#H
consH = aov(allcons$cons_H_u[allcons$sample_type == 'Trichoptera'] ~ caddislkfac)
write.csv(summary(consH)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consHcaddis.csv')
consHTukeyOut = pairwise.t.test(allcons$cons_H_u[allcons$sample_type == 'Trichoptera'], 
                                caddislkfac, p.adjust.method='bonferroni')
write.csv(consHTukeyOut$p.value,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/consHpwcaddis.csv')


#d2h2o (1 way anova) ####
d2 <- read.csv('C:/Users/Mike/git/seven_lakes/Analysis/7lmeans.csv')

lkfac = factor(c(1,1,2,2,2,2,3,3,3,1,3), labels=c('lg','sm','hi'))
aovmod = aov(d2$dD_H2O[1:11] ~ lkfac)
write.csv(summary(aovmod)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/waterH.csv')
write.csv(TukeyHSD(aovmod)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/waterHpw.csv')

#sources ####
source <- read.csv('cor_sources.csv')

#phyto
#N
phytoN = aov(source$mean_d15NA_lake ~ lkfac)
write.csv(summary(phytoN)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcNphyto.csv')
write.csv(TukeyHSD(phytoN)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcNpwphyto.csv')
#C
phytoC = aov(source$mean_d13CA_lake ~ lkfac)
write.csv(summary(phytoC)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcCphyto.csv')
write.csv(TukeyHSD(phytoC)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcCpwphyto.csv')
#H
phytoH = aov(source$mean_d2HA_lake ~ lkfac)
write.csv(summary(phytoH)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcHphyto.csv')
write.csv(TukeyHSD(phytoH)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcHpwphyto.csv')
#peri
#N
periN = aov(source$mean_d15Nperi_lake ~ lkfac)
write.csv(summary(periN)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcNperi.csv')
write.csv(TukeyHSD(periN)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcNpwperi.csv')
#C
periC = aov(source$mean_d13Cperi_lake ~ lkfac)
write.csv(summary(periC)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcCperi.csv')
write.csv(TukeyHSD(periC)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcCpwperi.csv')
#H
periH = aov(source$mean_d2Hperi_lake ~ lkfac)
write.csv(summary(periH)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcHperi.csv')
write.csv(TukeyHSD(periH)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcHpwperi.csv')
#terr
#N
terrN = aov(source$d15NT ~ lkfac)
write.csv(summary(terrN)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcNterr.csv')
write.csv(TukeyHSD(terrN)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcNpwterr.csv')
#C
terrC = aov(source$d13CT ~ lkfac)
write.csv(summary(terrC)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcCterr.csv')
write.csv(TukeyHSD(terrC)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcCpwterr.csv')
#H
terrH = aov(source$d2HT ~ lkfac)
write.csv(summary(terrH)[[1]],
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcHterr.csv')
write.csv(TukeyHSD(terrH)$lkfac,
          file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/table_2/srcHpwterr.csv')

