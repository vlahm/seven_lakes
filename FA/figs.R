#plots for ASLO conference, 6/2016
rm(list=ls())

#setup ####
# library(viridis)
library(RColorBrewer)
a <- brewer.pal(3, 'Dark2')
a[2:3] <- a[c(3,2)]
setwd("C:/Users/Mike/git/seven_lakes/Analysis")
means<-read.csv("7lmeans.csv")
se<-read.csv("7lse.csv")
error.bars <- function(x, y, upper, lower=upper, cap.length=0.1, horiz=F,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("One or more vectors is not the same length")

    if(horiz==F) {
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
    } else if (horiz==T) {
        arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
    }
}
replacer <- function(x, cond, replacement, rows=1:length(x[,1]), cols=1:length(x[1,]), asinsqrt=FALSE){
    for (i in rows)
    {
        for (j in cols)
        {
            if (cond=='==NA')
            {
                if (is.na(x[i,j]==TRUE))
                {
                    x[i,j] <- replacement
                }
            } else
            {
                if (eval(parse(text=paste('x[i,j]', cond))))
                {
                    x[i,j] <- replacement
                }
            }
            if (asinsqrt == TRUE)
            {
                x[i,j]<-asin(sqrt(x[i,j]/100))*(2/pi)
            }
        }
    }
    return(x)
}
library(vegan)

nmds.scree <-
    function(x,distance='bray',k=6,trymax=50,
             autotransform=FALSE,trace=0,...){

        library(vegan)
        library(MASS)

        old.par<-par(no.readonly=TRUE)

        nmds.stress<-rep(0,k)
        nmds.dim<-c(1:k)
        for(i in 1:k){
            y<-metaMDS(x,distance=distance,k=i,trymax=trymax,
                       autotransform=autotransform,trace=trace,...)
            nmds.stress[i]<-y$stress
        }
        plot(nmds.dim,nmds.stress,type='o',pch=19,col='blue',
             ylab='Stress',xlab='Ordination Axis',
             main='Scree Plot of Stress vs. Dimension',...)

        par(old.par)
    }

nmds.monte <-
    function(x,k,distance='bray',trymax=50,autotransform=FALSE,
             trace=0,zerodist='add',perm=100,col.hist='blue',col.line='red',
             lty=2,las=1,lab=c(5,5,4),...){

        library(vegan)
        library(MASS)

        z<-metaMDS(comm=x,k=k,distance=distance,trymax=trymax,
                   autotransform=autotransform,trace=trace,...) #nmds analysis
        z.stress<-z$stress #get stress
        y.stress<-rep(0,perm)

        for(i in 1:perm){
            y<-apply(x,2,sample) #permute data matrix
            y<-metaMDS(comm=y,k=k,distance=distance,trymax=trymax,
                       autotransform=autotransform,trace=trace,...) #nmds analysis
            y.stress[i]<-y$stress #get stress
        }
        n<-sum(y.stress<=z.stress) #compute number of random runs with stress < observed
        p.value<-(1+n)/(1+perm) #compute p-value

        xmin<-min(z.stress,min(y.stress))
        xmax<-max(z.stress,max(y.stress))
        hist(y.stress,col=col.hist,las=las,lab=lab,
             xaxs='i',yaxs='i',xlim=c(xmin,xmax),xlab='Stress',
             main=paste('Random Permutation Distribution of Stress for',k,'Dimensions',sep=' '),...)
        abline(v=z.stress,col=col.line,lty=lty,lwd=2,...)

        cat('Randomization Test of Stress:\n')
        cat('Permutation stress values:\n')
        print(y.stress)
        z<-rbind('Observed stress'=z.stress,'P-value'=p.value)
        return(z)
    }

pca.eigenval <-
    function(x.pca,dim=length(x.pca$sdev),digits=7){

        #check for dim limit
        if(dim>length(x.pca$sdev)){
            cat("Only",length(x.pca$sdev),"axes available\n")
            dim<-length(x.pca$sdev)
        }

        #calculate some variables
        names<-colnames(x.pca$rotation[,1:dim])
        var<-x.pca$sdev^2
        trace<-sum(var)
        prop.var<-var/trace

        #broken-stick distribution
        p<-length(x.pca$sdev)
        y<-rep(0,p)
        for(i in 1:p) y[i]<-1/i
        for(i in 1:p) y[i]<-sum(y[i:p])
        y<-y[1:dim]

        #print results
        cat('Importance of components:\n')
        z<-rbind('Variance(eigenvalue)'=var[1:dim],
                 'Proportion of Variance'=prop.var[1:dim],
                 'Cumulative Proportion'=cumsum(prop.var[1:dim]),
                 'Broken-stick value'=y)
        colnames(z)<-names
        z<-round(z,digits=digits)
        return(z)
    }

pca.eigenvec <-
    function(x.pca,dim=length(x.pca$sdev),
             digits=7,cutoff=0){

        #check for dim limit
        if(dim>ncol(x.pca$rotation)){
            cat("Only",ncol(x.pca$rotation),"axes available\n")
            dim<-ncol(x.pca$rotation)
        }

        #print results
        cat("\nEigenvectors:\n")
        z<-format(round(x.pca$rotation[,1:dim],digits=digits))
        z[abs(x.pca$rotation[,1:dim])<cutoff]<-substring('',1,nchar(z[1,1]))
        z<-as.data.frame(z)
        return(z)
    }

pca.structure <-
    function(x.pca,x,dim=length(x.pca$sdev),
             digits=3,cutoff=0){

        #check for dim limit
        if(dim>length(x.pca$sdev)){
            cat("Only",length(x.pca$sdev),"axes available\n")
            dim<-length(x.pca$sdev)
        }

        #calculate structure correlations
        z<-cor(x,x.pca$x[,1:dim])

        #print results
        cat("\nStructure Correlations:\n")
        z<-round(z,digits=digits)
        z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
        z<-as.data.frame(z)
        return(z)
    }

#chem plots ####

#elevxdocxlake area
# a <- viridis(3, end=0.8)

par(mfrow=c(1,2), mar=c(4,7,4,0))
plot(means$DOC[c(1,2,7:11)], means$elev[c(1,2,7:11)], xlim=c(0,1.5), ylim=c(1225,1500),
     cex=rescale(means$lake_area[1:11], c(2,15))[c(1,2,7:11)], xaxs='i', yaxs='i',
     bty='n', las=2, xaxt='n', xlab='', yaxt='n', ylab='',
     col=c(a[1],a[1],a[2],a[2],a[2],a[1],a[2]), lwd=4, type='n')
error.bars(x=means$DOC[c(1,2,7:11)], y=means$elev[c(1,2,7:11)],
           upper=se$DOC[c(1,2,7:11)], cap.length=0, horiz=TRUE, lwd=2, col='gray30')
points(means$DOC[c(1,2,7:11)], means$elev[c(1,2,7:11)], xlim=c(0,1.5), ylim=c(1225,1500),
     cex=rescale(means$lake_area[1:11], c(2,15))[c(1,2,7:11)], xaxs='i', yaxs='i',
     bty='n', las=2, xaxt='n', xlab='', yaxt='n', ylab='',
     col=c(a[1],a[1],a[2],a[2],a[2],a[1],a[2]), lwd=4)
axis(side=1, at=c(0,.2,.4,.6,.8,1,1.2,1.4), pos=1250, las=1, cex.axis=1.5)
axis(side=2, las=2, cex.axis=1.5)
mtext('Elevation (m)', font=2, side=2, outer=TRUE, line=-2, cex=2)
par(mar=c(4,0,4,4))
plot(means$DOC[3:6], means$elev[3:6], xlim=c(1.5,11), ylim=c(1225,1500), bty='n',
     cex=rescale(means$lake_area[1:11], c(2,15))[3:6], xaxs='i', yaxs='i', axes=FALSE,
     col=a[3], lwd=4, xlab='', type='n')
error.bars(x=means$DOC[3:6], y=means$elev[3:6],
           upper=se$DOC[3:6], cap.length=0, horiz=TRUE, lwd=2, col='gray30')
points(means$DOC[3:6], means$elev[3:6], xlim=c(1.5,11), ylim=c(1225,1500), bty='n',
     cex=rescale(means$lake_area[1:11], c(2,15))[3:6], xaxs='i', yaxs='i',
     col=a[3], lwd=4, xlab='')
axis(side=1, at=c(2,4,6,8,10), labels=c('2.0', '4.0', '6.0', '8.0', '10.0'), pos=1250,
     las=1, cex.axis=1.5)
mtext('Dissolved Organic Carbon (mg/L)', side=1, outer=TRUE, line=-2, font=2, cex=2)
legend('topright', legend=c('Montane Small', 'Montane Large', 'Alpine'), bty='n',
       border=NA, fill=c(a[3],a[1],a[2]), cex=1.6)

#DOC x color:chla x lake area
par(mfrow=c(2,2), mar=c(0,7,4,0))
plot(means$DOC[5:6], means$color_chlA[5:6], xlim=c(1.5,11.5),  bty='n', ylim=c(.05,.26),
     axes=FALSE, type='n', ylab='', xaxs='i', yaxs='i')
legend('top', legend=c('Montane Small', 'Montane Large', 'Alpine'), bty='n',
       border=NA, fill=c(a[3],a[1],a[2]), cex=1.6)
axis(side=2, las=2, cex.axis=1.5, at=c(.1,.15,.2,.25))
par(mar=c(0,0,4,4))
new <- se$abs_440*sqrt(2)/means$chlA^2 + se$chlA*sqrt(2)*means$abs_440^2/means$chlA^4
new[c(1,5)] <- c(2.471466e-02, 1e-3)
plot(means$DOC[5:6], means$color_chlA[5:6], xlim=c(1.5,11.5),  bty='n', ylim=c(.05,.26),
     cex=rescale(means$vegshed.lakearea[1:11], c(2,15))[5:6], xaxs='i', yaxs='i', axes=FALSE,
     col=a[3], lwd=4, xlab='', type='n')
error.bars(x=means$DOC[5:6], y=means$color_chlA[5:6],
           upper=new[5:6], cap.length=0, lwd=2, col='gray30')
error.bars(x=means$DOC[5:6], y=means$color_chlA[5:6],
           upper=se$DOC[5:6], cap.length=0, horiz=TRUE, lwd=2, col='gray30')
points(means$DOC[5:6], means$color_chlA[5:6], xlim=c(1.5,11.5),  bty='n', ylim=c(.05,.26),
     cex=rescale(means$vegshed.lakearea[1:11], c(2,15))[5:6], xaxs='i', yaxs='i',
     col=a[3], lwd=4, xlab='')
par(mar=c(7,7,0,0))
plot(means$DOC[c(1,2,7:11)], means$color_chlA[c(1,2,7:11)], xlim=c(0,1.5), ylim=c(0,.05),
     cex=rescale(means$vegshed.lakearea[1:11], c(2,15))[c(1,2,7:11)], xaxs='i', yaxs='i',
     bty='n', las=2, xaxt='n', xlab='', yaxt='n', ylab='',
     col=c(a[1],a[1],a[2],a[2],a[2],a[1],a[2]), lwd=4, type='n')
error.bars(x=means$DOC[c(1,2,7:11)], y=means$color_chlA[c(1,2,7:11)],
           upper=new[c(1,2,7:11)], cap.length=0, lwd=2, col='gray30', xpd=NA)
error.bars(x=means$DOC[c(1,2,7:11)], y=means$color_chlA[c(1,2,7:11)],
           upper=se$DOC[c(1,2,7:11)], cap.length=0, horiz=TRUE, lwd=2, col='gray30')
points(means$DOC[c(1,2,7:11)], means$color_chlA[c(1,2,7:11)], xlim=c(0,1.5), ylim=c(0,.05),
     cex=rescale(means$vegshed.lakearea[1:11], c(2,15))[c(1,2,7:11)], xaxs='i', yaxs='i',
     bty='n', las=2, xaxt='n', xlab='', yaxt='n', ylab='',
     col=c(a[1],a[1],a[2],a[2],a[2],a[1],a[2]), lwd=4)
axis(side=2, las=2, cex.axis=1.5)
axis(side=1, at=c(0,.4,.8,1.2), las=1, cex.axis=1.5)
mtext('a440nm (AU) : Chlorophyll-a (ug/L)', font=2, side=2, outer=TRUE,
      line=-2, cex=1.4)
par(mar=c(7,0,0,4))
plot(means$DOC[3:5], means$color_chlA[3:5], xlim=c(1.5,11.5),  bty='n', ylim=c(0,.05),
     cex=rescale(means$vegshed.lakearea[1:11], c(2,15))[3:5], xaxs='i', yaxs='i', axes=FALSE,
     col=a[3], lwd=4, xlab='', type='n')
error.bars(x=means$DOC[3:5], y=means$color_chlA[3:5],
           upper=new[3:5], cap.length=0, lwd=2, col='gray30')
error.bars(x=means$DOC[3:5], y=means$color_chlA[3:5],
           upper=se$DOC[3:5], cap.length=0, horiz=TRUE, lwd=2, col='gray30')
points(means$DOC[3:5], means$color_chlA[3:5], xlim=c(1.5,11.5),  bty='n', ylim=c(0,.05),
     cex=rescale(means$vegshed.lakearea[1:11], c(2,15))[3:5], xaxs='i', yaxs='i',
     col=a[3], lwd=4, xlab='')
axis(side=1, at=c(2,4,6,8,10), labels=c('2.0', '4.0', '6.0', '8.0', '10.0'),# pos=1250,
     las=1, cex.axis=1.5)
mtext('Dissolved Organic Carbon (mg/L)', side=1, outer=TRUE, line=-2, font=2, cex=2)

#see file "ordination covariates" for covariate PCA####

#NMDS####
FA_data <- read.csv("C:\\Users\\Mike\\Desktop\\Grad\\conferences\\FA_data.csv")
FA <- FA_data[1:60,12:42]
FA <- replacer(FA, '==NA', 0)
branched <- colSums(FA[which(FA_data$Branched == 1),])
long_safa <- colSums(FA[which(FA_data$Long.SAFA == 1),])
short_safa <- colSums(FA[which(FA_data$short.SAFA == 1),])
MOB <- colSums(FA[which(FA_data$MOB == 1),])
MUFA <- colSums(FA[which(FA_data$MUFA == 1),])
PUFA <- colSums(FA[which(FA_data$PUFA == 1),])

#gotta run a mixing model before this step if including MLEs
nmds_data <- cbind(t(rbind(branched,long_safa,short_safa,MOB,MUFA,PUFA)),
                   t(as.matrix(FA_data[61,12:42])))#,
                   # as.matrix(MLEs[3,]))
                   # t(MLEs))
# colnames(nmds_data)[7:8] <- c('w3-w6','terr_prop')
# colnames(nmds_data)[7:10] <- c('w3-w6','phyto_prop', 'peri_prop', 'terr_prop')
colnames(nmds_data)[7] <- c('w3-w6')

# nmds.scree(nmds_data, distance="bray", k=8, trymax=20)
# nmds.monte(nmds_data, distance="bray", k=2, trymax=100)

nmds <- metaMDS(nmds_data, distance="bray", k=2, trymax=100)
# stressplot(nmds)

#sample scores
scores <- nmds$points

#loadings
loadings<-envfit(nmds$points, nmds_data, perm=1000)

#plot
l <- a[1]
h <- a[2]
s <- a[3]
par(mar=c(5,5,4,4))
fig <- ordiplot(nmds, choices=c(1,2), type="none", cex.lab=1.6, cex.main=1.8,
                main="FA Ordination", xlim=c(-.6,.3), ylim=c(-.5,.5))
plot(loadings, p.max=.05, col="gray40", cex=1.3, arrow.mul=0.5)#, labels='')
points(fig, "sites", col=c(h,h,s,s,h,h,l,l,l,h,h,s,s,s,l,s,h,h,h,l,s,l,l,l,l,s,h,h,s,s,s),
       pch=c(3,3,17,17,17,17,17,17,17,17,17,17,16,16,17,16,16,16,16,4,16,15,16,16,15,15,16,4,15,15,15),
       cex=1.7, lwd=2)
legend('topleft', legend=c('Montane Small', 'Montane Large', 'Alpine', 'Calanoida', 'Cladocera',
                           'Trichoptera','Salvelinus','Chaoborus'), bty='n',
       border=NA, fill=c(a[3],a[1],a[2],rep(NA,5)), pch=c(NA,NA,NA,16,15,17,3,4),
       cex=1.2, pt.lwd=2)

#PCA original####
FA_data <- read.csv("D:\\Dropbox\\Grad\\presentations\\ASLO\\FA_data.csv")
FA <- FA_data[1:60,12:42]
FA <- replacer(FA, '==NA', 0)
branched <- colSums(FA[which(FA_data$Branched == 1),])
long_safa <- colSums(FA[which(FA_data$Long.SAFA == 1),])
short_safa <- colSums(FA[which(FA_data$short.SAFA == 1),])
MOB <- colSums(FA[which(FA_data$MOB == 1),])
MUFA <- colSums(FA[which(FA_data$MUFA == 1),])
PUFA <- colSums(FA[which(FA_data$PUFA == 1),])

pca_data <- cbind(t(rbind(branched,long_safa,short_safa,MOB,MUFA,PUFA)),
                   t(as.matrix(FA_data[61,12:42])))
colnames(pca_data)[7] <- c('w3-w6')

pca_data[,-7] <- replacer(pca_data[,-7], '==0', 0.00000001, asinsqrt=TRUE)

#format data for PCA
# tcal.asin.grouped<-t(as.matrix(cal.asin.grouped))

#perform PCA and associated tests
pca<-prcomp(pca_data, scale=T, scores=T, center=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which components are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,pca_data,dim=7,cutoff=0.2)
#sample scores
sample.scores<-pca$x[,1:7]

#PCA revamp####
FA_data <- read.csv("D:\\Dropbox\\Grad\\presentations\\ASLO\\FA_data_revamp.csv")
FA <- FA_data[1:60,13:43]
FA <- replacer(FA, '==NA', 0)
branched <- colSums(FA[which(FA_data$Branched == 1),])
long_safa <- colSums(FA[which(FA_data$Long.SAFA == 1),])
med_safa <- colSums(FA[which(FA_data$med.SAFA == 1),])
short_safa <- colSums(FA[which(FA_data$short.SAFA == 1),])
MOB <- colSums(FA[which(FA_data$MOB == 1),])
MUFA <- colSums(FA[which(FA_data$MUFA == 1),])
PUFA <- colSums(FA[which(FA_data$PUFA == 1),])

pca_data <- cbind(t(rbind(branched,long_safa,med_safa,short_safa,MOB,MUFA,PUFA)),
                   t(as.matrix(FA_data[61,13:43])))
colnames(pca_data)[8] <- c('w3-w6')

pca_data[,-8] <- replacer(pca_data[,-8], '==0', 0.00000001, asinsqrt=TRUE)

#format data for PCA
# tcal.asin.grouped<-t(as.matrix(cal.asin.grouped))

#perform PCA and associated tests
pca<-prcomp(pca_data, scale=T, center=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which components are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,pca_data,dim=7,cutoff=0.2)
#sample scores
sample.scores<-pca$x[,1:7]

#PCA revamp 4 ####
FA_data <- read.csv("D:\\Dropbox\\Grad\\presentations\\ASLO\\FA_data_revamp4.csv")
FA <- FA_data[1:60,13:43]
FA <- replacer(FA, '==NA', 0)
algPUFA <- colSums(FA[which(FA_data$algPUFA == 1),])
w36PUFA <- colSums(FA[which(FA_data$w36PUFA == 1),])
lSAFA <- colSums(FA[which(FA_data$lSAFA == 1),])
mSAFA <- colSums(FA[which(FA_data$mSAFA == 1),])
Bacteria <- colSums(FA[which(FA_data$Bacteria == 1),])
MOB <- colSums(FA[which(FA_data$MOB == 1),])
Other <- colSums(FA[which(FA_data$Other == 1),])

pca_data <- cbind(t(rbind(algPUFA,w36PUFA,lSAFA,mSAFA,Bacteria,MOB,Other)),
# pca_data <- cbind(t(rbind(algae,tOM,bact,MOB)),
                  t(as.matrix(FA_data[61,13:43])))
colnames(pca_data)[8] <- c('w3-w6')
# colnames(pca_data)[5] <- c('w3-w6')

pca_data[,-8] <- replacer(pca_data[,-8], '==0', 0.00000001, asinsqrt=TRUE)
# pca_data[,-5] <- replacer(pca_data[,-5], '==0', 0.00000001, asinsqrt=TRUE)

#format data for PCA
# tcal.asin.grouped<-t(as.matrix(cal.asin.grouped))

#perform PCA and associated tests
pca<-prcomp(pca_data, scale=T, center=T)
#determine eigenvalues (and proportion of variance)
eigenvalues<-pca.eigenval(pca)
#see which components are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,pca_data,dim=5,cutoff=0.0)
# apply(as.matrix(structure),2,as.numeric)^2
#sample scores
sample.scores<-pca$x[,1:4]

#PCA revamp 3 ####
FA_data <- read.csv("D:\\Dropbox\\Grad\\presentations\\ASLO\\FA_data_revamp3.csv")
FA <- FA_data[1:60,10:40]
FA <- replacer(FA, '==NA', 0)
algae <- colSums(FA[which(FA_data$Algae == 1),])
tOM <- colSums(FA[which(FA_data$t.OM == 1),])
bact <- colSums(FA[which(FA_data$Bacteria == 1),])
other <- colSums(FA[which(FA_data$Other == 1),])

pca_data <- cbind(t(rbind(algae,tOM,bact,other)),
# pca_data <- cbind(t(rbind(algae,tOM,bact)),
                  t(as.matrix(FA_data[61,10:40])))
colnames(pca_data)[5] <- c('w3-w6')
# colnames(pca_data)[4] <- c('w3-w6')

pca_data[,-5] <- replacer(pca_data[,-5], '==0', 0.00000001, asinsqrt=TRUE)
# pca_data[,-4] <- replacer(pca_data[,-4], '==0', 0.00000001, asinsqrt=TRUE)

#format data for PCA
# tcal.asin.grouped<-t(as.matrix(cal.asin.grouped))

#perform PCA and associated tests
pca<-prcomp(pca_data, scale=T, center=T)
#determine eigenvalues (and proportion of variance)
eigenvalues<-pca.eigenval(pca)
#see which components are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,pca_data,dim=5,cutoff=0.0)
# apply(as.matrix(structure),2,as.numeric)^2
#sample scores
sample.scores<-pca$x[,1:4]

##plots (individual) [used in paper] ####
# png(width=7, height=6, units='in', type='cairo', res=1200,
#     file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/fig5_FA.png')
pdf(width=7, height=6, compress=FALSE,
    file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/fig5_FA.pdf')
# library(viridis)
# cols <- plasma(3, begin=0, end=0.75)
# l <- cols[1]; h <- cols[2]; s <- cols[3]
l <- 'springgreen3'; h <- 'darkgoldenrod2'; s <- 'royalblue2'
n <- 1.3; o <- 1.55; p <- 2
par(mar=c(5,5,2,2))
fig<-ordiplot(pca, choices=c(1,2), type="none", main="", cex.main=2, cex.lab=1.7,
              xlab='PC1 (47.0%)', ylab='PC2 (23.3%)', las=1, cex.axis=1.3)#, xlim=c(-3,2), ylim=c(-2.8,2.6))
arrows(0,0,pca$rotation[,1]*3.7, pca$rotation[,2]*3.7, col="gray50", length=0.1, lwd=2,
       lty=c(1,1,2,1,1,1,1,1))
points(fig, "sites", lwd=3, cex=1.7,
       pch=NA,
       # pch=c(3,3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,21,21,21,4,NA,NA,21,21,NA,NA,21,4,NA,NA,NA),
       # pch=c(3,3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,21,21,NA,21,21,21,21,4,21,23,21,21,23,23,21,4,23,23,23),
       # pch=c(3,3,24,24,24,24,24,24,24,24,24,24,21,21,24,21,21,21,21,4,21,23,21,21,23,23,21,4,23,23,23),
       col=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s),
       bg=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s))
       #### cex=c(o,o, n, n, o, o, p, p, p, o, o, n, n, n, p, n, o, o, o, p, n, p,p, p, p, n, o, o, n, n, n))
# text(pca$rotation[,1]*2, pca$rotation[,2]*2, row.names(pca$rotation), col="purple")
text(c(-2.1,0.7,2.2,2.0,2.245,0.1,-2.5,-1.8),
     c(1.1,1.1,0.1,-0.3,-0.5,-2.6,-0.9,0.3),
     labels=c('Algae','t-OM (nonsig.)','mSAFA','sSAFA,','Other MUFA','(MO) Bacteria',
              'w3,w6 PUFA','w3:w6'), col='gray50',font=2)
legend('bottomleft', legend=c('Montane Small', 'Montane Large', 'Alpine'),
       fill=c(s,l,h), border='gray50', bty='n', cex=1.2)
legend('topleft', legend=c('Calanoida', 'Cladocera', 'Trichoptera','Salvelinus','Chaoborus'),
       pch=c(21,23,24,3,4), cex=1.2, pt.lwd=2, pt.bg='black', bty='n')
# legend('bottomleft', legend=c('Montane-Small', 'Montane-Large', 'Alpine', 'Calanoida', 'Cladocera',
#                               'Trichoptera','Salvelinus','Chaoborus'),
#        fill=c(s,l,h,rep(NA,5)), pch=c(NA,NA,NA,21,23,24,3,4),
#        cex=1.2, pt.lwd=3, border=NA, pt.bg='black', bty='n')
dev.off()
shell('D:/Dropbox/Grad/Projects/Thesis/Seven^ Lakes^ Project^ 2014/Manuscript/figures/fig5_FA.pdf',
      translate=TRUE)

#pcs 1 and 3
l <- 'springgreen3'; h <- 'darkgoldenrod2'; s <- 'royalblue2'
n <- 1.3; o <- 1.55; p <- 2
par(mar=c(5,5,2,2))
fig<-ordiplot(pca, choices=c(1,3), type="none", main="", cex.main=2, cex.lab=1.7,
              xlab='PC1 (47.0%)', ylab='PC3 (%)', las=1, cex.axis=1.3)#, xlim=c(-3,2), ylim=c(-2.8,2.6))
arrows(0,0,pca$rotation[,1]*3.7, pca$rotation[,3]*3.7, col="gray50", length=0.1, lwd=2)
points(fig, "sites", lwd=3, cex=1.7,
       pch=c(3,3,24,24,24,24,24,24,24,24,24,24,21,21,24,21,21,21,21,4,21,23,21,21,23,23,21,4,23,23,23),
       col=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s),
       bg=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s))
       # cex=c(o,o, n, n, o, o, p, p, p, o, o, n, n, n, p, n, o, o, o, p, n, p,p, p, p, n, o, o, n, n, n))

text(pca$rotation[,1]*2, pca$rotation[,3]*2, row.names(pca$rotation), col="purple")
legend('bottomleft', legend=c('Montane Small', 'Montane Large', 'Alpine'),
       fill=c(s,l,h), border='gray30', bty='n', cex=1.2)
legend('topleft', legend=c('Calanoida', 'Cladocera', 'Trichoptera','Salvelinus','Chaoborus'),
       pch=c(21,25,24,3,4), cex=1.2, pt.lwd=2, pt.bg='black', bty='n')
# legend('bottomleft', legend=c('Montane-Small', 'Montane-Large', 'Alpine', 'Calanoida', 'Cladocera',
#                               'Trichoptera','Salvelinus','Chaoborus'),
#        fill=c(s,l,h,rep(NA,5)), pch=c(NA,NA,NA,21,23,24,3,4),
#        cex=1.2, pt.lwd=3, border=NA, pt.bg='black', bty='n')


#first and third principal components
# l <- 'springgreen3'; h <- 'darkgoldenrod2'; s <- 'royalblue2'
# n <- 1.3; o <- 1.55; p <- 2
# par(mar=c(5,5,2,2))
# fig<-ordiplot(pca, choices=c(1,3), type="none", main="", cex.main=2, cex.lab=1.7,
#               xlab='PC1 (48.8%)', ylab='PC3 (9.7%)', las=1, cex.axis=1.3, xlim=c(-2.5,2.2), ylim=c(-2.4,2.2))
# arrows(0,0,pca$rotation[,1]*2.5, pca$rotation[,3]*2.5, col="gray50", length=0.1, lwd=2)
# points(fig, "sites", lwd=2, cex=1.7,
#        pch=c(3,3,24,24,24,24,24,24,24,24,24,24,21,21,24,21,21,21,21,4,21,25,21,21,25,25,21,4,25,25,25),
#        col=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s),
#        bg=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s))
# # cex=c(o,o, n, n, o, o, p, p, p, o, o, n, n, n, p, n, o, o, o, p, n, p,p, p, p, n, o, o, n, n, n))
#
# text(pca$rotation[,1]*2, pca$rotation[,3]*2, row.names(pca$rotation), col="purple")
# legend('topright', legend=c('Montane Small', 'Montane Large', 'Alpine', 'Calanoida', 'Cladocera',
#                               'Trichoptera','Salvelinus','Chaoborus'),
#        fill=c(s,l,h,rep(NA,5)), pch=c(NA,NA,NA,21,25,24,3,4),
#        cex=1.2, pt.lwd=2, border=NA, pt.bg='black', bty='n')


# plots (together) ####
#first and second principal components
l <- 'springgreen3'; h <- 'darkgoldenrod2'; s <- 'royalblue2'
n <- 1.3; o <- 1.55; p <- 2
par(mfrow=c(1,2),mar=c(5,5,1,1), oma=c(0,0,0,6))
fig<-ordiplot(pca, choices=c(1,2), type="none", main="", cex.main=2, cex.lab=1.7,
              xlab='', ylab='PC2 (32.9%)', las=1, cex.axis=1.3, xlim=c(-3.3,2.7), ylim=c(-3.4,2.9))
arrows(0,0,pca$rotation[,1]*3.6, pca$rotation[,2]*3.6, col="gray50", length=0.1, lwd=2)
points(fig, "sites", lwd=2, cex=1.7,
       pch=c(3,3,24,24,24,24,24,24,24,24,24,24,21,21,24,21,21,21,21,4,21,25,21,21,25,25,21,4,25,25,25),
       col=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s),
       bg=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s))
text(pca$rotation[,1]*2, pca$rotation[,2]*2, row.names(pca$rotation), col="purple")

#first and third principal components
l <- 'springgreen3'; h <- 'darkgoldenrod2'; s <- 'royalblue2'
n <- 1.3; o <- 1.55; p <- 2
fig<-ordiplot(pca, choices=c(1,3), type="none", main="", cex.main=2, cex.lab=1.7,
              xlab='', ylab='PC3 (9.7%)', las=1, cex.axis=1.3, xlim=c(-3.3,2.7), ylim=c(-1,1))
arrows(0,0,pca$rotation[,1]*2.5, pca$rotation[,3]*2.5, col="gray50", length=0.1, lwd=2)
points(fig, "sites", lwd=2, cex=1.7,
       pch=c(3,3,24,24,24,24,24,24,24,24,24,24,21,21,24,21,21,21,21,4,21,25,21,21,25,25,21,4,25,25,25),
       col=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s),
       bg=c(h,h, s, s, h, h, l, l, l, h, h, s, s, s, l, s, h, h, h,l, s, l, l, l, l, s, h,h, s, s, s))
# cex=c(o,o, n, n, o, o, p, p, p, o, o, n, n, n, p, n, o, o, o, p, n, p,p, p, p, n, o, o, n, n, n))
text(pca$rotation[,1]*2, pca$rotation[,3]*2, row.names(pca$rotation), col="purple")

mtext('PC1 (48.8%)', side=1, cex=1.7, outer=TRUE, line=-2)
legend('topright', legend=c('Montane Small', 'Montane Large', 'Alpine', 'Calanoida', 'Cladocera',
                            'Trichoptera','Salvelinus','Chaoborus'),
       fill=c(s,l,h,rep(NA,5)), pch=c(NA,NA,NA,21,25,24,3,4),
       cex=1.2, pt.lwd=2, border=NA, pt.bg='black', bty='n')

#allochth x production plot####
a <- brewer.pal(3, 'Dark2')
a[2:3] <- a[c(3,2)]
l <- a[1]
h <- a[2]
s <- a[3]
par(mar=c(5,5,4,4))
plot(covs$chlA, MLEs[3,], cex=1.5,
     xlab='Chlorophyll-a (ug/L)', ylab='Proportion Allochthony', cex.main=2, cex.lab=1.7,
     col=c(h,h,s,s,h,h,l,l,l,h,h,s,s,s,l,s,h,h,h,l,s,l,l,l,l,s,h,h,s,s,s), lwd=2,
     pch=c(3,3,17,17,17,17,17,17,17,17,17,17,16,16,17,16,16,16,16,4,16,15,16,16,15,15,16,4,15,15,15))
mod <- lm(MLEs[3,] ~ covs$chlA)
# mod2 <- lm(MLEs[3,] ~ covs$chlA + I(covs$chlA^-2))
summary(mod)
# summary(mod2)
abline(mod)
# plot(covs$chlA, fitted(mod2))

plot(covs$O2_Ar, MLEs[3,], cex=1.5,
     xlab=expression(paste(O[2], ':Ar')), ylab='Proportion Allochthony',
     col=c(h,h,s,s,h,h,l,l,l,h,h,s,s,s,l,s,h,h,h,l,s,l,l,l,l,s,h,h,s,s,s), lwd=2,
     pch=c(3,3,17,17,17,17,17,17,17,17,17,17,16,16,17,16,16,16,16,4,16,15,16,16,15,15,16,4,15,15,15))


#
