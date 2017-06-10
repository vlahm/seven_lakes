#Michael Vlah - University of Washington - Seven Lakes Project
#vlahm13@gmail.com
#created August 2015
#last edit: 9/15/15

#Notice: this script includes many trial methods that became obsolete.
#They have been retained for future reference and code-borrowing.
#To run the analysis in its final form, you'll need to run the material in
#part 0 - functions and packages.  Then start with part 8.

dev.off()
windows(record=T) #navigate plots in this separate window using PgUp and PgDn
setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")
data<-read.csv("ultimate_spreadsheet_POMandPeriCombinedLakes.csv")
########0 -functions and packages####
alldata<-read.csv("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Analysis/7lmeans.csv")
library(vioplot)
library(vegan)
library(plyr)
library(plotrix)
library(car)

#FISH560 functions
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

ordi.monte <-
  function(x,ord,dim=length(x),perm=1000,center=TRUE,
           scale=TRUE,digits=3,plot=TRUE,col.hist='blue',col.line='red',
           lty=2,las=1,lab=c(5,5,4),...){

    p<-length(x)
    if(dim>p){
      cat("Only",p,"axes available\n")
      dim<-p
    }

    if(ord=='pca'){
      z<-prcomp(x,center=center,scale=scale) #prcomp analysis
      z.eig<-z$sdev[1:dim]^2 #compute eigenvalues
      z.teig<-t(z.eig) #transpose eigenvalues
      z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
      write('',file='y.csv') #empty outfile if it exists
      for(i in 1:perm){
        y<-apply(x,2,sample) #permute data matrix
        y<-prcomp(y,center=center,scale=scale) #prcomp on permuted matrix
        y<-y$sdev[1:dim]^2 #compute eigenvalues
        y<-as.data.frame(t(y)) #coerce to data frame and transpose
        write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
      }
      y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
      p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
      p.value<-p.value/perm #compute p-value
      names<-colnames(z$rotation[,1:dim]) #add 'PC#' names
    }

    else if(ord=='ca'){
      library(vegan)
      z<-cca(x) #correspondence analysis
      z.eig<-z$CA$eig[1:dim] #get eigenvalues
      z.teig<-t(z.eig) #transpose eigenvalues
      z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
      write('',file='y.csv') #empty outfile if it exists
      for(i in 1:perm){
        y<-apply(x,2,sample) #permute data matrix
        y<-cca(y) #CA on permuted matrix
        y<-y$CA$eig[1:dim] #get eigenvalues
        y<-as.data.frame(t(y)) #coerce to data frame and transpose
        write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
      }
      y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
      p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
      p.value<-p.value/perm #compute p-value
      names<-names(z$CA$eig[1:dim]) #add 'CA#' names
    }

    else if(ord=='dca'){
      library(vegan)
      if(dim>4){
        cat("Only",4,"axes available\n")
        dim<-4
      }
      z<-decorana(x,...) #detrended correspondence analysis
      z.eig<-z$evals[1:dim] #get eigenvalues
      z.teig<-t(z.eig) #transpose eigenvalues
      z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
      write('',file='y.csv') #empty outfile if it exists
      for(i in 1:perm){
        y<-apply(x,2,sample) #permute data matrix
        y<-decorana(y,...) #DCA on permuted matrix
        y<-y$evals[1:dim] #get eigenvalues
        y<-as.data.frame(t(y)) #coerce to data frame and transpose
        write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
      }
      y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
      p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
      p.value<-p.value/perm #compute p-value
      names<-names(z$eval[1:dim]) #add 'CA#' names
    }

    if(plot==TRUE){
      for(i in 1:dim){
        xmin<-min(min(y[[i]],z.eig[i]))
        xmax<-max(max(y[[i]],z.eig[i]))
        hist(y[[i]],col=col.hist,las=las,lab=lab,
             xaxs='i',yaxs='i',xlim=c(xmin,xmax),
             xlab='Eigenvalue',
             main=paste('Random Permutation Distribution of Eigenvalues for',names[i],sep=' '),...)
        abline(v=z.eig[i],col=col.line,lty=lty,lwd=2,...)
        readline("Press return for next plot ")
      }
    }

    cat('Randomization Test of Eigenvalues:\n')
    z<-rbind('Eigenvalue'=z.eig,'P-value'=p.value)
    colnames(z)<-names
    z<-round(z,digits=digits)
    return(z)
  }

#1 - FA_PCAer####

############NOTICE
#this version of FA_PCAer should be set up to calculate percent_area correctly
#(based on the overall total from each sample, rather than just the subsample used in
#the PCA.  However, the spreadsheets have not been set up to run the analysis.
#This vein of the project was abandoned in favor of FA groupings.

# FA_PCAer<-function(FA.subset, all.FAs, ordi.monte=F)
# {
#   #read in the FAs and samples you want to test
#   subset<-read.csv(FA.subset)
#   all<-read.csv(all.FAs) #this is necessary for computing proportional areas
#
#   #find number and indices of columns containing Areal data
#   Nsamp<-length(grep("Area",colnames(subset)))
#   samp.ind<-grep("Area",colnames(subset))
#   proportional.area<-matrix(data=NA, nrow=length(subset[,1]), ncol=Nsamp)
#
#   #create matrix of proportional area
#   for(i in 1:Nsamp)
#   {
#     subset.col<-subset[,samp.ind[i]]
#     all.col<-all[,samp.ind[i]]
#     total.area<-sum(na.omit(all.col))
#     for(j in 1:length(subset.col))
#     {
#       proportional.area[j,i]<-(subset.col[j]*100)/total.area #find out why i and j are being created as objects
#     }
#     NAs<-which(is.na(proportional.area[,i]))
#     proportional.area[,i][NAs]<-0.0000001
#   }
#   proportional.area<-cbind(subset[,2:5],proportional.area)
#   colnames(proportional.area)[5:length(proportional.area[1,])]<-gsub("Reten_", "", colnames(csv)[grep("Reten", colnames(csv))])
#   row.names(proportional.area)<-subset[,1]
#
#   #the PCA...
#   transposed<-t(as.matrix(proportional.area[,-(1:4)]))
#   pca<-prcomp(transposed,scale=T, scores=T)
#   #determine eigenvalues
#   eigenvalues<-pca.eigenval(pca)
#   #see which eigenvalues are significant
#   screeplot(pca, bstick=T)
#   #do that a different way
#   if(ordi.monte==T)
#   {
#     ordimonte<-o
#   }
#   #see loadings.  square these to get percentage of variance in each original variable
#   #accounted for by each principal component
#   structure<-pca.structure(pca,transposed,dim=7,cutoff=0.2)
#   #sample.scores<-pca$x[,1:7]
#   #first plotting method
#   biplot(pca)
#   #second plotting method
#   ordiplot(pca,choices=c(1,2), type="text", display="sites")
#   arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="red")
#   text(pca$rotation[,1]*5.2, pca$rotation[,2]*5.2, row.names(pca$rotation))
#
#   if(ordi.monte==T)
#   {
#     details<-list(eigenvalues[,1:7], ordimonte, structure)
#   }
#   else
#   {
#     details<-list(eigenvalues[,1:7], structure)
#   }
#   return(details)
# }
#
# FA_PCAer("1_allConsumers_FAoverPointFive.csv", )
# FA_PCAer("2_allConsumers_FAoverPointSeven.csv")
# FA_PCAer("4_caddisonly_29FAs.csv")
# FA_PCAer("5_caddisonly_20FAs.csv")
# FA_PCAer("6_calonly_43FAs.csv")
# FA_PCAer("7_cladonly_28FAs.csv")
#2 - Inefficient Method (1 of 3) Grouped PCAs - probably ignore - Calanoid####

setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")

raw.calanoid<-read.csv("8_cal_grouped.csv")
calanoid<-t(raw.calanoid[,-(1:2)])
colnames(calanoid)<-raw.calanoid[,1]

cal.pca<-prcomp(calanoid,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(cal.pca)
#see which eigenvalues are significant
screeplot(cal.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(cal.pca,calanoid,dim=7,cutoff=0.2)
#sample.scores<-cal.pca$x[,1:7]
#first plotting method
biplot(cal.pca, main="calanoid")
#second plotting method
ordiplot(cal.pca,choices=c(1,2), type="text", display="sites", main="Calanoid")
arrows(0,0,cal.pca$rotation[,1]*5, cal.pca$rotation[,2]*5, col="blue")
text(cal.pca$rotation[,1]*5.2, cal.pca$rotation[,2]*5.2, row.names(cal.pca$rotation), col="blue")

#comparison of distributions:
vioplot(calanoid[1,], calanoid[2,], calanoid[3,], calanoid[4,], calanoid[5,],
        calanoid[6,], calanoid[7,], calanoid[8,], calanoid[9,], calanoid[10,], names=(rownames(calanoid)))
#3 - Caddis####

raw.caddis<-read.csv("9_caddis_grouped.csv")
caddis<-t(raw.caddis[,-(1:2)])
colnames(caddis)<-raw.caddis[,1]

caddis.pca<-prcomp(caddis,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(caddis.pca)
#see which eigenvalues are significant
screeplot(caddis.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(caddis.pca,caddis,dim=7,cutoff=0.2)
#sample.scores<-caddis.pca$x[,1:7]
#first plotting method
biplot(caddis.pca, main="caddis")
#second plotting method
ordiplot(caddis.pca,choices=c(1,2), type="text", display="sites", main="Caddisfly", ylim=c(-3,3))
arrows(0,0,caddis.pca$rotation[,1]*4, caddis.pca$rotation[,2]*4, col="blue")
text(caddis.pca$rotation[,1]*4.2, caddis.pca$rotation[,2]*4.2, row.names(caddis.pca$rotation), col="blue")

#comparison of distributions.  Do they seem to differ in the same way that the calanoid ones do?
vioplot(caddis[1,], caddis[2,], caddis[3,], caddis[4,], caddis[5,],
        caddis[6,], caddis[7,], caddis[8,], caddis[9,], caddis[10,],
        caddis[11,], names=(rownames(caddis)))
#Should run a k-sample K-S test here.  do this for cal and clad too? will implement later.

# rows<-rep(1:11, times=length(caddis[1,]))
# cols<-rep(1:7, each=length(caddis[,1]))
# counts<-round(stack(as.data.frame(caddis))[,1]*10000) #can only use integers in contingency tables, so here's
# #a "conversion" to integers that allows some retention of sigdigs.
# ctgcy<-as.data.frame(cbind(counts,rows,cols))
# ctgcy[,2]<-factor(ctgcy[,2])
# ctgcy[,3]<-factor(ctgcy[,3])
#4 - Cladocera####

raw.clado<-read.csv("10_clado_grouped.csv")
clado<-t(raw.clado[,-(1:2)])
colnames(clado)<-raw.clado[,1]

clado.pca<-prcomp(clado,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(clado.pca)
#see which eigenvalues are significant
screeplot(clado.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(clado.pca,clado,dim=7,cutoff=0.2)
#sample.scores<-clado.pca$x[,1:7]
#first plotting method
biplot(clado.pca, main="cladocera")
#second plotting method
ordiplot(clado.pca,choices=c(1,2), type="text", display="sites", main="Cladocera", ylim=c(-3,3))
arrows(0,0,clado.pca$rotation[,1]*5, clado.pca$rotation[,2]*5, col="blue")
text(clado.pca$rotation[,1]*5.2, clado.pca$rotation[,2]*5.2, row.names(clado.pca$rotation), col="blue")

#comparison of distributions
vioplot(clado[1,], clado[2,], clado[3,], clado[4,], clado[5,],
        clado[6,], names=(rownames(clado)))
#5 - Cladocera and chaoborus####

raw.cc<-read.csv("11_cladANDchaob_grouped.csv")
cc<-t(raw.cc[,-(1:2)])
colnames(cc)<-raw.cc[,1]

cc.pca<-prcomp(cc,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(cc.pca)
#see which eigenvalues are significant
screeplot(cc.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(cc.pca,cc,dim=7,cutoff=0.2)
#sample.scores<-cc.pca$x[,1:7]
#first plotting method
biplot(cc.pca, main="clad + chaob")
#second plotting method
ordiplot(cc.pca,choices=c(1,2), type="text", display="sites", main="Cladocera + Chaoborus", ylim=c(-3,3))
arrows(0,0,cc.pca$rotation[,1]*5, cc.pca$rotation[,2]*5, col="blue")
text(cc.pca$rotation[,1]*5.2, cc.pca$rotation[,2]*5.2, row.names(cc.pca$rotation), col="blue")

#comparison of distributions
vioplot(cc[1,], cc[2,], cc[3,], cc[4,], cc[5,],
        cc[6,], cc[7,], cc[8,], names=(rownames(cc)))
#6 - all together####
setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements/misc spreadsheets")
all<-read.csv("4a_all_samples_unprocessed.csv")

#store the first five columns (info columns)
info.cols<-all[,1:5]
#remove spacer column
#not necessary yet

#create matrix of relative area of each FA in each sample
area.cols<-grep("Area", colnames(all))
proportional.area<-matrix(data=NA, nrow=length(all[,1]), ncol=length(area.cols))
for(i in 1:length(area.cols))
{
  working.col<-all[,area.cols[i]]
  total.area<-sum(na.omit(working.col))
  for(j in 1:length(working.col))
  {
    proportional.area[j,i]<-(working.col[j]*100)/total.area
  }
}

#replace column names with IDs (ID key here: "C:\Users\Mike\Desktop\Grad\Projects\Thesis\Seven Lakes Project 2014\Data\FA\PCA_arrangements\ID_key.csv")
#first three digits = sample number, next one is lake type, next two are organism type, last two correspond to lake name
IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
       "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081318", "i082320", "i083420", "i084418",
       "i093118", "i105320", "i106317", "i097118", "i091118", "i095118", "i017111", "i032110", "i029410", "i034210", "i041108", "i044308", "i048308",
       "i051208", "i054408", "i064108", "i071109", "i07300f", "i07400b", "i07500h", "i07600s", "i07700c", "i109107")
colnames(proportional.area)<-IDs

#identify columns with no data / 1 datum
# rowSums(is.na(proportional.area)) == 48
# rowSums(is.na(proportional.area)) == 47

#combine FAs according to biomarker classes
#switching to Excel

write.csv(proportional.area, "C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements/misc spreadsheets/4b_converted_to_rel.csv")
#7 - all together - semi-efficient Method (2 of 3)####
#added sorting columns in excel
setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")
# data<-read.csv("ultimate_spreadsheet.csv")
# IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
#        "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081318", "i082320", "i083420", "i084418",
#        "i093118", "i105320", "i106317", "i097118", "i091118", "i095118", "i017111", "i032110", "i029410", "i034210", "i041108", "i044308", "i048308",
#        "i051208", "i054408", "i064108", "i071109", "i00000t", "i109107")
# colnames(data)[17:60]<-IDs
#
# #create groupings
# short_SAFA<-data[data$short_SAFA==1,c(1:5, 17:60)] #isolate rows
#   short_SAFA<-colSums(short_SAFA[,7:49], na.rm=T) #sum columns
# iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, 17:60)]
#   iSAFA_etc<-colSums(iSAFA_etc[,7:49], na.rm=T)
# other_MUFA<-data[data$other_MUFA==1,c(1:5, 17:60)]
#   other_MUFA<-colSums(other_MUFA[,7:49], na.rm=T)
# LIN<-data[data$LIN==1,c(1:5, 17:60)]
#   LIN<-colSums(LIN[,7:49], na.rm=T)
# other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, 17:60)]
#   other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,7:49], na.rm=T)
# OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, 17:60)]
#   OA_16.4n3<-colSums(OA_16.4n3[,7:49], na.rm=T)
# ALA_SDA<-data[data$ALA_SDA==1,c(1:5, 17:60)]
#   ALA_SDA<-colSums(ALA_SDA[,7:49], na.rm=T)
# long_SAFA<-data[data$long_SAFA==1,c(1:5, 17:60)]
#   long_SAFA<-colSums(long_SAFA[,7:49], na.rm=T)
# ARA<-data[data$ARA==1,c(1:5, 17:60)]
#   ARA<-colSums(ARA[,7:49], na.rm=T)
# EPA_DHA<-data[data$EPA_DHA==1,c(1:5, 17:60)]
#   EPA_DHA<-colSums(EPA_DHA[,7:49], na.rm=T)
# other_PUFA<-data[data$other_PUFA==1,c(1:5, 17:60)]
#   other_PUFA<-colSums(other_PUFA[,7:49], na.rm=T)
#
# grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
#       ALA_SDA, long_SAFA, ARA, EPA_DHA, other_PUFA)
#
# #need to combine samples for periphyton and POM.  too many zero values in those
# #categories.  Is it legal to combine, or are the distributions very dissimilar?
# #comparing probability distributions of periphyton
# qqplot(data[,49], data[,50]); abline(0,1)
# qqplot(data[,49], data[,51]); abline(0,1)
# qqplot(data[,51], data[,50]); abline(0,1)
# #probably not cool to combine periphyton.  Comparing POM.
# qqplot(data[,52], data[,53]); abline(0,1)
# qqplot(data[,52], data[,54]); abline(0,1)
# qqplot(data[,52], data[,55]); abline(0,1) #similar
# qqplot(data[,52], data[,56]); abline(0,1) #similar
# qqplot(data[,52], data[,57]); abline(0,1) #only two values
# qqplot(data[,53], data[,54]); abline(0,1)
# qqplot(data[,53], data[,55]); abline(0,1)
# qqplot(data[,53], data[,56]); abline(0,1)
# qqplot(data[,53], data[,57]); abline(0,1) #only two values
# qqplot(data[,54], data[,55]); abline(0,1)
# qqplot(data[,54], data[,56]); abline(0,1)
# qqplot(data[,54], data[,57]); abline(0,1) #only two values
# qqplot(data[,55], data[,56]); abline(0,1) #similar
# qqplot(data[,55], data[,57]); abline(0,1) #only two values
# qqplot(data[,56], data[,57]); abline(0,1) #only two values
# #probably not cool to combine.  Gonna do it anyway.

#combined POM and peri samples in Excel
data<-read.csv("ultimate_spreadsheet_POMandPeriCombinedLakes.csv")
IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
       "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081301", "i082320", "i083420", "i084401",
       "i093101", "i105320", "i106317", "i097101", "i091101", "i095101", "i017111", "i000504", "i000508", "i071109", "i000514", "i109107")
colnames(data)[16:52]<-IDs

#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, 16:52)] #isolate rows
short_SAFA<-colSums(short_SAFA[,6:42], na.rm=T) #sum columns
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, 16:52)]
iSAFA_etc<-colSums(iSAFA_etc[,6:42], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, 16:52)]
other_MUFA<-colSums(other_MUFA[,6:42], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, 16:52)]
LIN<-colSums(LIN[,6:42], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, 16:52)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:42], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, 16:52)]
OA_16.4n3<-colSums(OA_16.4n3[,6:42], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, 16:52)]
ALA_SDA<-colSums(ALA_SDA[,6:42], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, 16:52)]
long_SAFA<-colSums(long_SAFA[,6:42], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, 16:52)]
EPA_DHA<-colSums(EPA_DHA[,6:42], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, 16:52)]
other_PUFA<-colSums(other_PUFA[,6:42], na.rm=T)

#grouping 1
grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(grouped[1,]))
{
  for (j in 1:length(grouped[,1]))
  {
    if (grouped[j,i]==0)
    {
    grouped[j,i]<-0.000001
    }
  }
}

#PCA, ignoring "other_PUFA" (which has many zeros and is not very informative;
#note that it could be combined with "other_n3_n6_PUFA", which also has a lot of zeros.)

tgrouped<-t(as.matrix(grouped))
pca<-prcomp(tgrouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which eigenvalues are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,tgrouped,dim=7,cutoff=0.2)
#sample.scores<-pca$x[,1:7]
#first plotting method
biplot(pca)
#second plotting method
fig1<-ordiplot(pca, choices=c(1,2), type="none", xlim=c(-4,4), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))  #good example of argument matching for later use (match.arg)
arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="black")
text(pca$rotation[,1]*5.8, pca$rotation[,2]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

########8 - PCA and HCA functions, etc. (Method 3 of 3) - NEED THIS TO RUN ALL PCAS and HCAs####

#combined POM and peri samples in Excel

# IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010403", "i011203", "i012203", "i013103",
#        "i090120", "i092120", "i108303", "i094120", "i085220", "i086220", "i087220", "i080317", "i096120", "i081301", "i082320", "i083420", "i084401",
#        "i093101", "i105320", "i106317", "i097101", "i091101", "i095101", "i017111", "i000504", "i000508", "i071109", "i000514", "i109107")
# colnames(data)[16:52]<-IDs

#identify columns for each sample type
cal.cols<-grep("20", substr(colnames(data),6,7))
clad.cols<-grep("01", substr(colnames(data),6,7))
caddis.cols<-grep("03", substr(colnames(data),6,7))
chaob.cols<-grep("17", substr(colnames(data),6,7))
fish.cols<-grep("15", substr(colnames(data),6,7))
POM.cols<-grep("08", substr(colnames(data),6,7))
peri.cols<-grep("04", substr(colnames(data),6,7))
fil.cols<-grep("11", substr(colnames(data),6,7))
moss.cols<-grep("07", substr(colnames(data),6,7))
soil.cols<-grep("09", substr(colnames(data),6,7))
plant.cols<-grep("14", substr(colnames(data),6,7))

#replace NAs with small value (for frequency distribution plots)
dist.data<-data
for (i in 1:length(dist.data[1,16:52]))
{
  for (j in 1:length(dist.data[,16]))
  {
    if (is.na(dist.data[,16:52][j,i])==T)
    {
      dist.data[,16:52][j,i]<-0.000001
    }
  }
}
####9 - Calanoid + producers####
#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:21], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:21], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:21], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:21], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:21], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:21], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:21], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:21], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:21], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:21], na.rm=T)

#grouping 1
cal.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA")
# cal.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#modify the data
for (i in 1:length(cal.grouped[1,]))
{
  for (j in 1:length(cal.grouped[,1]))
  {
    #replace zeros with small values
    if (cal.grouped[j,i]==0)
    {
      cal.grouped[j,i]<-0.000001
    }
    #put data on proportional scale (0-1) and arc-sin square root transform them
    cal.asin.grouped<-cal.grouped
    cal.asin.grouped[j,i]<-asin(sqrt(cal.asin.grouped[j,i]/100))*(2/pi)
  }
}

tcal.asin.grouped<-t(as.matrix(cal.asin.grouped))
cal.pca<-prcomp(tcal.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(cal.pca)
#see which eigenvalues are significant
screeplot(cal.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(cal.pca,tcal.asin.grouped,dim=7,cutoff=0.2)
#sample.scores<-cal.pca$x[,1:7]
#first and second principal components
cal.fig<-ordiplot(cal.pca, choices=c(1,2), type="none", xlim=c(-6,3), ylim=c(-5,4), main="Cal PC1 and 2")
points(cal.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)))
arrows(0,0,cal.pca$rotation[,1]*5, cal.pca$rotation[,2]*5, col="black")
text(cal.pca$rotation[,1]*5.8, cal.pca$rotation[,2]*5.8, row.names(cal.pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


#plotting 3rd and 4th principal components
cal.fig<-ordiplot(cal.pca, choices=c(3,4), type="none", xlim=c(-4,4), ylim=c(-4,4), main="Cal PC3 and 4")
points(cal.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)))
arrows(0,0,cal.pca$rotation[,3]*5, cal.pca$rotation[,4]*5, col="black")
text(cal.pca$rotation[,3]*5.8, cal.pca$rotation[,4]*5.8, row.names(cal.pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

#are there obvious differences among the frequency distributions?

# #vio.data is the same as dist.data
# vioplot(vio.data[,cal.cols[1]], vio.data[,cal.cols[2]], vio.data[,cal.cols[3]], vio.data[,cal.cols[4]], vio.data[,cal.cols[5]],
#         vio.data[,cal.cols[6]], vio.data[,cal.cols[7]], vio.data[,cal.cols[8]], vio.data[,cal.cols[9]], vio.data[,cal.cols[10]],
#         names=(colnames(vio.data[,cal.cols])), h=3, wex=1)

par(mar=c(0,5,0,0))
par(mfrow=c(10,1))
for(i in 1:10)
{
  if (substr(colnames(dist.data[cal.cols[i]]),5,5)==1)
  {
    plot(dist.data[,cal.cols[i]], type="l", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="blue")
  }
  else if (substr(colnames(dist.data[cal.cols[i]]),5,5)==2)
  {
    plot(dist.data[,cal.cols[i]], type="l", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="red")
  }
  else if (substr(colnames(dist.data[cal.cols[i]]),5,5)==3)
  {
    plot(dist.data[,cal.cols[i]], type="l", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="black")
  }
  else if (substr(colnames(dist.data[cal.cols[i]]),5,5)==4)
  {
    plot(dist.data[,cal.cols[i]], type="l", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="green")
  }
}
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))
####10 - Caddisfly + producers####

#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:22], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:22], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:22], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:22], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:22], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:22], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:22], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:22], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:22], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:22], na.rm=T)

#grouping 1
caddis.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# caddis.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(caddis.grouped[1,]))
{
  for (j in 1:length(caddis.grouped[,1]))
  {
    if (caddis.grouped[j,i]==0)
    {
      caddis.grouped[j,i]<-0.000001
      #put data on proportional scale (0-1) and arc-sin square root transform them
      caddis.asin.grouped<-caddis.grouped
      caddis.asin.grouped[j,i]<-asin(sqrt(caddis.asin.grouped[j,i]/100))*(2/pi)
    }
  }
}

tcaddis.asin.grouped<-t(as.matrix(caddis.asin.grouped))
caddis.pca<-prcomp(tcaddis.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(caddis.pca)
#see which eigenvalues are significant
screeplot(caddis.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(caddis.pca,tcaddis.asin.grouped,dim=7,cutoff=0.2)
#sample.scores<-caddis.pca$x[,1:7]

#first and second principal components
caddis.fig<-ordiplot(caddis.pca, choices=c(1,2), type="none", xlim=c(-5,3), ylim=c(-6,3), main="Caddis PC1 and 2")
points(caddis.fig, "sites", pch=as.numeric(substr(row.names(tcaddis.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)))
arrows(0,0,caddis.pca$rotation[,1]*5, caddis.pca$rotation[,2]*5, col="black")
text(caddis.pca$rotation[,1]*5.8, caddis.pca$rotation[,2]*5.8, row.names(caddis.pca$rotation))
legend("bottomleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomright", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


#plotting 2nd and 3rd principal components
caddis.fig<-ordiplot(caddis.pca, choices=c(2,3), type="none", xlim=c(-5,3), ylim=c(-4,3), main="Caddis PC2 and 3")
points(caddis.fig, "sites", pch=as.numeric(substr(row.names(tcaddis.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)))
arrows(0,0,caddis.pca$rotation[,2]*5, caddis.pca$rotation[,3]*5, col="black")
text(caddis.pca$rotation[,2]*5.8, caddis.pca$rotation[,3]*5.8, row.names(caddis.pca$rotation))
legend("bottomleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomright", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

#plotting 1st and 3rd principal components
caddis.fig<-ordiplot(caddis.pca, choices=c(1,3), type="none", xlim=c(-5,3), ylim=c(-4,3), main="Caddis PC1 and 3")
points(caddis.fig, "sites", pch=as.numeric(substr(row.names(tcaddis.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)))
arrows(0,0,caddis.pca$rotation[,1]*5, caddis.pca$rotation[,3]*5, col="black")
text(caddis.pca$rotation[,1]*5.8, caddis.pca$rotation[,3]*5.8, row.names(caddis.pca$rotation))
legend("bottomleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                              "fish", "POM", "periphyton", "filamentous", "moss",
                              "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                               11, 7, 9, 14))
legend("bottomright", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

#are there obvious differences among the frequency distributions?
par(mar=c(0,5,0,0))
par(mfrow=c(11,1))
for(i in 1:11)
{
  if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==1)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="blue")
  }
  else if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==2)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="red")
  }
  else if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==3)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="black")
  }
  else if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==4)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="green")
  }
}
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))
####11 - Cladocera + producers####
#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:17], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:17], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:17], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:17], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:17], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:17], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:17], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:17], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:17], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:17], na.rm=T)

#grouping 1
clad.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# clad.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(clad.grouped[1,]))
{
  for (j in 1:length(clad.grouped[,1]))
  {
    if (clad.grouped[j,i]==0)
    {
      clad.grouped[j,i]<-0.000001
      #put data on proportional scale (0-1) and arc-sin square root transform them
      clad.asin.grouped<-clad.grouped
      clad.asin.grouped[j,i]<-asin(sqrt(clad.asin.grouped[j,i]/100))*(2/pi)
    }
  }
}

tclad.asin.grouped<-t(as.matrix(clad.asin.grouped))
clad.pca<-prcomp(tclad.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(clad.pca)
#see which eigenvalues are significant
screeplot(clad.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(clad.pca,tclad.asin.grouped,dim=7,cutoff=0.2)
#sample.scores<-clad.pca$x[,1:7]
#first and second principal components
clad.fig<-ordiplot(clad.pca, choices=c(1,2), type="none", xlim=c(-4,5), ylim=c(-4,5), main="Cladocera PC1 and 2")
points(clad.fig, "sites", pch=as.numeric(substr(row.names(tclad.asin.grouped),6,7)), col=as.numeric(substr(row.names(tclad.asin.grouped),5,5)))
arrows(0,0,clad.pca$rotation[,1]*5, clad.pca$rotation[,2]*5, col="black")
text(clad.pca$rotation[,1]*5.8, clad.pca$rotation[,2]*5.8, row.names(clad.pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


#plotting 3rd and 4th principal components
clad.fig<-ordiplot(clad.pca, choices=c(3,4), type="none", xlim=c(-4,5), ylim=c(-4,4), main="Cladocera PC3 and 4")
points(clad.fig, "sites", pch=as.numeric(substr(row.names(tclad.asin.grouped),6,7)), col=as.numeric(substr(row.names(tclad.asin.grouped),5,5)))
arrows(0,0,clad.pca$rotation[,3]*5, clad.pca$rotation[,4]*5, col="black")
text(clad.pca$rotation[,3]*5.8, clad.pca$rotation[,4]*5.8, row.names(clad.pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

#are there obvious differences among the frequency distributions?
par(mar=c(0,5,0,0))
par(mfrow=c(6,1))
for(i in 1:6)
{
  if (substr(colnames(dist.data[clad.cols[i]]),5,5)==1)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="blue")
  }
  else if (substr(colnames(dist.data[clad.cols[i]]),5,5)==2)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="red")
  }
  else if (substr(colnames(dist.data[clad.cols[i]]),5,5)==3)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="black")
  }
  else if (substr(colnames(dist.data[clad.cols[i]]),5,5)==4)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="green")
  }
}
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))
####12 - All####

short_SAFA<-data[data$short_SAFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:42], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:42], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:42], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:42], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:42], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:42], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:42], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:42], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:42], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:42], na.rm=T)

#grouping 1
all.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# all.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(all.grouped[1,]))
{
  for (j in 1:length(all.grouped[,1]))
  {
    if (all.grouped[j,i]==0)
    {
      all.grouped[j,i]<-0.000001
      #put data on proportional scale (0-1) and arc-sin square root transform them
      all.asin.grouped<-all.grouped
      all.asin.grouped[j,i]<-asin(sqrt(all.asin.grouped[j,i]/100))*(2/pi)
    }
  }
}

tall.asin.grouped<-t(as.matrix(all.asin.grouped))
all.pca<-prcomp(tall.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(all.pca)
#see which eigenvalues are significant
screeplot(all.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(all.pca,tall.asin.grouped,dim=7,cutoff=0.2)
#sample.scores<-all.pca$x[,1:7]
#first and second principal components
all.fig<-ordiplot(all.pca, choices=c(1,2), type="none", xlim=c(-4,5), ylim=c(-4,4), main="All PC1 and 2")
points(all.fig, "sites", pch=as.numeric(substr(row.names(tall.asin.grouped),6,7)), col=as.numeric(substr(row.names(tall.asin.grouped),5,5)))
arrows(0,0,all.pca$rotation[,1]*5, all.pca$rotation[,2]*5, col="black")
text(all.pca$rotation[,1]*5.8, all.pca$rotation[,2]*5.8, row.names(all.pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


# #plotting 3rd and 4th principal components
# all.fig<-ordiplot(all.pca, choices=c(3,4), type="none", xlim=c(-4,5), ylim=c(-4,4), main="All PC1 and 2")
# points(all.fig, "sites", pch=as.numeric(substr(row.names(tall.asin.grouped),6,7)), col=as.numeric(substr(row.names(tall.asin.grouped),5,5)))
# arrows(0,0,all.pca$rotation[,3]*5, all.pca$rotation[,4]*5, col="black")
# text(all.pca$rotation[,3]*5.8, all.pca$rotation[,4]*5.8, row.names(all.pca$rotation))
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                            "fish", "POM", "periphyton", "filamentous", "moss",
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
#                                                             11, 7, 9, 14))
# legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
#        fill=c(1, 2, 3, 4, 5))
#

########13 - NMDS functions - must run PCAs before corresponding NMDSs will work####
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
####14 - NMDS Calanoid ####

#this performs the NMDS
tcal.grouped<-t(as.matrix(cal.grouped))
cal.nmds<-metaMDS(tcal.grouped, distance="bray", k=2, trymax=100, autotransform=F) #see nmds.monte for dimension testing

#not sure how scree or stressplot are used here, but keep trying the monte until you have a
#managable number of dimensions and the last addition removes considerable stress.
nmds.scree(tcal.grouped, distance="bray", k=10, autotransform=F, trymax=20)
nmds.monte(tcal.grouped, distance="bray", k=2, autotransform=F, trymax=20) #stress is higher with 3 dimensions than with 2 or 4.  I went with 2 for the NMDS.
#based on the monte, 3 is the correct number of dimensions to consider
stressplot(cal.nmds)

#gangster plotting
# plot(cal.nmds, type="n")
# text(cal.nmds, labels=row.names(tcal.grouped))

#sample scores
cal.sample.scores<-cal.nmds$points

#calculate loadings (variable weights) on each NMDS axis
#these determine which FAs can be used to interpret the positions
#of the samples in ordination space.
cal.loadings<-envfit(cal.nmds$points, tcal.grouped, perm=1000)

#plot loadings on ordination plot
cal.fig2<-ordiplot(cal.nmds, choices=c(1,2), type="none", main="Cal NMDS1 and 2")
plot(cal.loadings, p.max=.05, col="black")
points(cal.fig2, "sites", pch=as.numeric(substr(row.names(tcal.grouped),6,7)), col=as.numeric(substr(row.names(tcal.grouped),5,5)))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

# #different method (not working)
# library(ggplot2)
# library(grid)
# cal.nmds<-metaMDS(tcal.grouped, distance="bray", k=3, trymax=20, autotransform=F)#copied from above
# cal.nmds2<-data.frame(MDS1 = cal.nmds$points[,1], MDS2 = cal.nmds$points[,2])
# cal.loadings<-envfit(cal.nmds$points, tcal.grouped, perm=1000)#copied from above
# cal.vec<-as.data.frame(cal.loadings$vectors$arrows*sqrt(cal.loadings$vectors$r))
# cal.vec$FAs<-rownames(cal.vec) #this step is unnecessary.  use rownames
# ggplot(data=cal.nmds2, aes(MDS1, MDS2))+
#   geom_point(aes(data=rownames(cal.nmds2), color=substr(rownames(cal.nmds2),5,5)))+
#   geom_segment(data=cal.vec, aes(c=0, xend=MDS1, y=0, yend=MDS2),
#                arrow = arrow(length=unit(0.5, "cm")), colour="gray", inherit_aes=F)+
#   geom_text(data=cal.vec, aes(x=MDS1, y=MDS2, label=FAs), size=5)+
#   coord_fixed()

# #third method (almost works as-is.  each legend is created in a different way and only lake_type
# #works right.  don't know how to label it though.  for the other one, could label within the scale_shape_identity call.
# library(ggplot2)
# library(grid)
# cal.nmds3<-metaMDS(tcal.grouped, distance="bray", k=3, trymax=20, autotransform=F)
# cal.scores<-as.data.frame(scores(cal.nmds3,display="sites"))
# cal.scores<-cbind(cal.scores, Lake_Type=substr(rownames(cal.scores),5,5), Organism=as.integer(substr(rownames(cal.scores),6,7)))
# #translate colname codes into lake and organism types for plotting
# # lake.types<-mapvalues(factor(substr(rownames(cal.scores),5,5)), from=c(1,2,3,4,5), to=c("Puddle", "High", "Clear + Milk", "Big", "NA"))
# # org.types<-mapvalues(factor(substr(rownames(cal.scores),6,7)), from=c("04","07","08","09","11","14","20"), to=c("Periphyton", "Moss", "POM", "Soil", "Filamentous", "Terr. Plant", "Calanoid"))
# cal.vf<-envfit(cal.nmds3, tcal.grouped, perm=1000)
# cal.spp<-as.data.frame(scores(cal.vf, display="vectors"))
# cal.spp<-cbind(cal.spp, FAs=rownames(cal.spp))
# cal.nmds.plot<-ggplot(cal.scores)+
#   geom_point(mapping=aes(x=NMDS1, y=NMDS2, colour=Lake_Type, shape=Organism))+
# #   scale_color_manual(values=c("darkslategray2", "darkseagreen2", "lightgoldenrod2", "darkorange2", "indianred2"),
# #                     name="Lake Type",
#   scale_shape_identity(guide="legend")+ #tells "shape" to accept integers as direct specifiers of pch
#   coord_fixed()+ #need aspect ratio of 1?
#   geom_segment(data=cal.spp,
#                aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
#                arrow=arrow(length=unit(0.25, "cm")), colour="grey")+
#   geom_text(data=cal.spp, aes(x=NMDS1, y=NMDS2, label=FAs), size=3)

# #combining clustering and ordination?
# cal.d<-vegdist(tcal.grouped, "bray")
# cal.ward<-hclust(cal.d, method="ward")
# cal.class<-cutree(cal.ward, k=3)
# cal.groups<-levels(factor(cal.class))
# cal.site<-scores(cal.nmds)
# ordiplot(cal.site, type="n", main="Cal NMDS + Clustering") #error here
# for (i in 1:length(cal.groups))
# {
#   points(cal.site[cal.class==i,], pch=(14+i), cex=2, col=i+1)
# }
# text(cal.site, row.names(#??????
####15 - NMDS Caddis ####

#this performs the NMDS
tcaddis.grouped<-t(as.matrix(caddis.grouped))
caddis.nmds<-metaMDS(tcaddis.grouped, distance="bray", k=2, trymax=100, autotransform=F) #again, 2 seems to be the best number of dimensions, as long as I'm interpreting my monte correctly

#not sure how scree or stressplot are used here, but keep trying the monte until you have a
#managable number of dimensions and the last addition removes considerable stress.
nmds.scree(tcaddis.grouped, distance="bray", k=10, autotransform=F, trymax=20)
nmds.monte(tcaddis.grouped, distance="bray", k=2, autotransform=F, trymax=20)
#based on the monte, 3 is the correct number of dimensions to consider
stressplot(caddis.nmds)

#gangster plotting
# plot(caddis.nmds, type="n")
# text(caddis.nmds, labels=row.names(tcaddis.grouped))

#sample scores
caddis.sample.scores<-caddis.nmds$points

#calculate loadings (variable weights) on each NMDS axis
#these determine which FAs can be used to interpret the positions
#of the samples in ordination space.
caddis.loadings<-envfit(caddis.nmds$points, tcaddis.grouped, perm=1000); caddis.loadings

#plot NMDS 1 and 2
caddis.fig2<-ordiplot(caddis.nmds, choices=c(1,2), type="none", main="Caddis NMDS1 and 2")
plot(caddis.loadings, p.max=.05, col="black")
points(caddis.fig2, "sites", pch=as.numeric(substr(row.names(tcaddis.grouped),6,7)), col=as.numeric(substr(row.names(tcaddis.grouped),5,5)))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))
####16 - NMDS Cladocera ####

#this performs the NMDS
tclad.grouped<-t(as.matrix(clad.grouped))
clad.nmds<-metaMDS(tclad.grouped, distance="bray", k=2, trymax=100, autotransform=F) #again, 2 seems to be the best number of dimensions, as long as I'm interpreting my monte correctly

#not sure how scree or stressplot are used here, but keep trying the monte until you have a
#managable number of dimensions and the last addition removes considerable stress.
nmds.scree(tclad.grouped, distance="bray", k=10, autotransform=F, trymax=20)
nmds.monte(tclad.grouped, distance="bray", k=2, autotransform=F, trymax=20)
#based on the monte, 3 is the correct number of dimensions to consider
stressplot(clad.nmds)

#gangster plotting
# plot(clad.nmds, type="n")
# text(clad.nmds, labels=row.names(tclad.grouped))

#sample scores
clad.sample.scores<-clad.nmds$points

#calculate loadings (variable weights) on each NMDS axis
#these determine which FAs can be used to interpret the positions
#of the samples in ordination space.
clad.loadings<-envfit(clad.nmds$points, tclad.grouped, perm=1000); clad.loadings

#plot NMDS 1 and 2
clad.fig2<-ordiplot(clad.nmds, choices=c(1,2), type="none", main="clad NMDS1 and 2")
plot(clad.loadings, p.max=.06, col="black")
points(clad.fig2, "sites", pch=as.numeric(substr(row.names(tclad.grouped),6,7)), col=as.numeric(substr(row.names(tclad.grouped),5,5)))
legend("topright", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomright", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))
####17 - NMDS all (gotta run PCA-all first)####

#this performs the NMDS
tall.grouped<-t(as.matrix(all.grouped))
all.nmds<-metaMDS(tall.grouped, distance="bray", k=3, trymax=100, autotransform=F) #three dimensions this time

#not sure how scree or stressplot are used here, but keep trying the monte until you have a
#managable number of dimensions and the last addition removes considerable stress.
nmds.scree(tall.grouped, distance="bray", k=10, autotransform=F, trymax=20)
nmds.monte(tall.grouped, distance="bray", k=3, autotransform=F, trymax=20)
#based on the monte, 3 is the correct number of dimensions to consider
stressplot(all.nmds)

#gangster plotting
# plot(all.nmds, type="n")
# text(all.nmds, labels=row.names(tall.grouped))

#sample scores
all.sample.scores<-all.nmds$points

#calculate loadings (variable weights) on each NMDS axis
#these determine which FAs can be used to interpret the positions
#of the samples in ordination space.
all.loadings<-envfit(all.nmds$points, tall.grouped, perm=1000); all.loadings

#plot NMDS 1 and 2
all.fig2<-ordiplot(all.nmds, choices=c(1,2), type="none", main="All NMDS1 and 2")
plot(all.loadings, p.max=.08, col="black")
points(all.fig2, "sites", pch=as.numeric(substr(row.names(tall.grouped),6,7)), col=as.numeric(substr(row.names(tall.grouped),5,5)))
legend("topright", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                            "fish", "POM", "periphyton", "filamentous", "moss",
                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                             11, 7, 9, 14))
legend("bottomright", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

####18 - Calanoid PCA etc. ####

#create FA groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, cal.cols)]
short_SAFA<-colSums(short_SAFA[,6:15], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, cal.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:15], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, cal.cols)]
other_MUFA<-colSums(other_MUFA[,6:15], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, cal.cols)]
LIN<-colSums(LIN[,6:15], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, cal.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:15], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, cal.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:15], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, cal.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:15], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, cal.cols)]
long_SAFA<-colSums(long_SAFA[,6:15], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, cal.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:15], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, cal.cols)]
other_PUFA<-colSums(other_PUFA[,6:15], na.rm=T)

#grouping combination 1 (all)
cal.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
                   ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA")
# cal.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#modify the grouped data for PCA (arcsin sqrt transform and missing data replacement)
for (i in 1:length(cal.grouped[1,]))
{
  for (j in 1:length(cal.grouped[,1]))
  {
    #replace zeros with small values
    if (cal.grouped[j,i]==0)
    {
      cal.grouped[j,i]<-0.000001
    }
    #put data on proportional scale (0-1) and arc-sin square root transform them
    cal.asin.grouped<-cal.grouped
    cal.asin.grouped[j,i]<-asin(sqrt(cal.asin.grouped[j,i]/100))*(2/pi)
  }
}

#format data for PCA
tcal.asin.grouped<-t(as.matrix(cal.asin.grouped))

#perform PCA and associated tests
cal.pca<-prcomp(tcal.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(cal.pca)
#see which eigenvalues are significant
screeplot(cal.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(cal.pca,tcal.asin.grouped,dim=7,cutoff=0.2)
#sample scores
sample.scores<-cal.pca$x[,1:7]

#create matrix of other covariates organized so that they can be used as cex values in the plot
cal.indexed.covariates<-matrix(NA, nrow=length(row.names(tcal.asin.grouped)), ncol=45)
for (i in 1:length(cal.indexed.covariates[1,]))
{
  for (j in 1:length(row.names(tcal.asin.grouped)))
  {
    #read lake_ID(but not class) off alldata and match it to lake ID in the indexed covariates
    cal.indexed.covariates[j,i]<-alldata[1:11,i][substr(alldata$lake_ID_andClass[1:11],1,2)==(as.numeric(substr(row.names(tcal.asin.grouped),8,9))[j])]
  }
}
cal.indexed.covariates<-as.data.frame(cal.indexed.covariates)
colnames(cal.indexed.covariates)<-colnames(alldata)
cal.indexed.covariates[,1]<-c("C", "O", "Z", "Mirror", "Y025", "Y015", "L", "Clear", "Morgenroth", "Milk")

#plot PCA
#first and second principal components (showing watershed area:lake area ratio)
cal.fig<-ordiplot(cal.pca, choices=c(1,2), type="none", xlim=c(-4,3), ylim=c(-3,4), main="Calanoid PCA; size = color:chlA")
points(cal.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)),
       col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)), cex=rescale(cal.indexed.covariates$color_chlA, c(2,15))) #rescales all variables to range: 2-15
# points(cal.fig, "sites", pch=3, bg="transparent", col="gray", cex=rescale(cal.indexed.covariates$lake_area, c(2,15)))
# points(cal.fig, "sites", pch=21, bg="transparent", col="gray", cex=rescale(cal.indexed.covariates$shed_lake_areaRatio, c(2,15)))
arrows(0,0,cal.pca$rotation[,1]*5, cal.pca$rotation[,2]*5, col="purple")
text(cal.pca$rotation[,1]*5.8, cal.pca$rotation[,2]*5.8, row.names(cal.pca$rotation), col="purple")
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                            "fish", "POM", "periphyton", "filamentous", "moss",
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4, 11, 7, 9, 14))
legend("bottomleft", legend=c("sm low forested", "med mid forested", "lg low forested", "high unforested"),
       fill=c(1, 3, 4, 2))
# legend("bottom", legend=c("gray cross = lake area", "gray circle = watershed:area"))


#to loop through with all covariates represented as point size:
#NOTE caddisd15n is misrepresented.  Also note issue with CH4
for (i in 3:length(cal.indexed.covariates[1,]))
{
  cal.fig<-ordiplot(cal.pca, choices=c(1,2), type="none", xlim=c(-6,3), ylim=c(-5,4), main=paste("Cal PC1 and 2, size = ", colnames(cal.indexed.covariates)[i]))
  points(cal.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)),
         col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)), cex=rescale(cal.indexed.covariates[,i], c(2,15))) #rescales all variables to range: 2-15
  points(cal.fig, "sites", pch=3, bg="transparent", col="gray", cex=rescale(cal.indexed.covariates$lake_area, c(2,15)))
  points(cal.fig, "sites", pch=21, bg="transparent", col="gray", cex=rescale(cal.indexed.covariates$shed_lake_areaRatio, c(2,15)))
  arrows(0,0,cal.pca$rotation[,1]*5, cal.pca$rotation[,2]*5, col="purple")
  text(cal.pca$rotation[,1]*5.8, cal.pca$rotation[,2]*5.8, row.names(cal.pca$rotation), col="purple")
  legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                             "fish", "POM", "periphyton", "filamentous", "moss",
                             "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                              11, 7, 9, 14))
  legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
         fill=c(1, 2, 3, 4, 5))
  legend("bottom", legend=c("gray cross = lake area", "gray circle = watershedA:lakeA"))
}

#notable covariates: CO2, CH4,

#plot cal d13c vs cal d15n (should show inverse relationship).  (then color by lake type) and size by other covariates
#until you find one that follows the trendline in either direction.

for(i in 3:length(cal.indexed.covariates[1,]))
{
  plot(cal.indexed.covariates$calanoidd13c, cal.indexed.covariates$calanoidd15n, cex=rescale(cal.indexed.covariates[,i], c(2,15)),
       main=paste(colnames(cal.indexed.covariates)[i]), col=substr(cal.indexed.covariates$lake_ID_andClass,3,3), pch=20, ylim = c(1,5),
       xlab = "Calanoid d13C", ylab = "Calanoid d15N")
  text(x = cal.indexed.covariates$calanoidd13c, y = cal.indexed.covariates$calanoidd15n + 0.4, labels = cal.indexed.covariates$lake)
  legend("topright", legend=c("puddles", "high", "milk+clear", "big"),
         fill=c(1, 2, 3, 4))
}


#plotting 3rd and 4th principal components
cal.fig<-ordiplot(cal.pca, choices=c(3,4), type="none", xlim=c(-4,4), ylim=c(-4,4), main="Cal PC3 and 4")
points(cal.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)))
arrows(0,0,cal.pca$rotation[,3]*5, cal.pca$rotation[,4]*5, col="black")
text(cal.pca$rotation[,3]*5.8, cal.pca$rotation[,4]*5.8, row.names(cal.pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

#are there obvious differences among the frequency distributions?

par(mar=c(0,5,0,0))
par(mfrow=c(10,1))
for(i in 1:10)
{
  if (substr(colnames(dist.data[cal.cols[i]]),5,5)==1)
  {
    plot(dist.data[,cal.cols[i]], type="o", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="black", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[cal.cols[i]]),5,5)==2)
  {
    plot(dist.data[,cal.cols[i]], type="o", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="red", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[cal.cols[i]]),5,5)==3)
  {
    plot(dist.data[,cal.cols[i]], type="o", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="green", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[cal.cols[i]]),5,5)==4)
  {
    plot(dist.data[,cal.cols[i]], type="o", ylab=colnames(dist.data[cal.cols[i]]), xlab="", col="blue", xaxt="n", ylim=c(0,25))
    axis(side=1, at=c(4,6,8,10,12,16,22,26,28,33,34,35,38,41,42,43,50,52,56,58,60),
         labels=dist.data$FA_name[c(4,6,8,10,12,16,22,26,28,33,34,35,38,41,42,43,50,52,56,58,60)], las=2)
  }
}
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))

#### plotting covariates with PC1
for(i in 3:length(cal.indexed.covariates[1,]))
{
  plot(sample.scores[,1], cal.indexed.covariates[,i], xlim=c(-4,4),
       main=paste(colnames(cal.indexed.covariates)[i]), col=substr(cal.indexed.covariates$lake_ID_andClass,3,3), pch=20, cex=3,
       xlab = "PC1", ylab = paste(colnames(cal.indexed.covariates)[i]))
  text(x = sample.scores[,1]+0.5, y = cal.indexed.covariates[,i], labels = cal.indexed.covariates$lake)
#   legend("topright", legend=c("puddles", "high", "milk+clear", "big"),
#          fill=c(1, 2, 3, 4))
}

for(i in 3:length(cal.indexed.covariates[1,]))
{
  plot(cal.indexed.covariates$color_chlA, cal.indexed.covariates[,i],
       main=paste("color:chl-a X ", colnames(cal.indexed.covariates)[i]), col=substr(cal.indexed.covariates$lake_ID_andClass,3,3), pch=20, cex=3,
       xlab = "color:chlorophyll-a", ylab = paste(colnames(cal.indexed.covariates)[i]))
  text(x = cal.indexed.covariates$color_chlA + 0.008, y = cal.indexed.covariates[,i], labels = cal.indexed.covariates$lake)
}

#test for correlations among the covariates that seem to diverge between the ponds and the other lakes when plotted against PC1
cor(cal.indexed.covariates[,c(4,6,11,14,24,29,34)])
scatterplotMatrix(cal.indexed.covariates[,c(4,6,11,14,24,29,34)], reg.line=F, smooth=F, spread=F, diagonal='density', id.n=0, span=0.5)
#see log for 9/22 to see how I resolved this

#### plotting covariates with PC1 JUST FOR THE NON-PUDDLES
for(i in 3:length(cal.indexed.covariates[1,]))
{
  plot(sample.scores[,1][-c(1,2,3,7)], cal.indexed.covariates[,i][-c(1,2,3,7)], xlim=c(-4,4),
       main=paste(colnames(cal.indexed.covariates)[i]), col=substr(cal.indexed.covariates$lake_ID_andClass,3,3)[-c(1,2,3,7)], pch=20, cex=3,
       xlab = "PC1", ylab = paste(colnames(cal.indexed.covariates)[i]))
  text(x = sample.scores[,1][-c(1,2,3,7)]+0.5, y = cal.indexed.covariates[,i][-c(1,2,3,7)], labels = cal.indexed.covariates$lake[-c(1,2,3,7)])
}

scatterplotMatrix(cal.indexed.covariates[-c(1,2,3,7),c(6,11,15,16,29,34)], reg.line=F, smooth=F, spread=F, diagonal='density', id.n=0, span=0.5)
####19 - Caddis PCA etc. ####

#create FA groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, caddis.cols)]
short_SAFA<-colSums(short_SAFA[,6:16], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, caddis.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:16], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, caddis.cols)]
other_MUFA<-colSums(other_MUFA[,6:16], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, caddis.cols)]
LIN<-colSums(LIN[,6:16], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, caddis.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:16], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, caddis.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:16], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, caddis.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:16], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, caddis.cols)]
long_SAFA<-colSums(long_SAFA[,6:16], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, caddis.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:16], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, caddis.cols)]
other_PUFA<-colSums(other_PUFA[,6:16], na.rm=T)

#grouping combination 1 (all)
caddis.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
                   ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA")
# caddis.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#modify the grouped data for PCA (arcsin sqrt transform and missing data replacement)
for (i in 1:length(caddis.grouped[1,]))
{
  for (j in 1:length(caddis.grouped[,1]))
  {
    #replace zeros with small values
    if (caddis.grouped[j,i]==0)
    {
      caddis.grouped[j,i]<-0.000001
    }
    #put data on proportional scale (0-1) and arc-sin square root transform them
    caddis.asin.grouped<-caddis.grouped
    caddis.asin.grouped[j,i]<-asin(sqrt(caddis.asin.grouped[j,i]/100))*(2/pi)
  }
}

#format data for PCA
tcaddis.asin.grouped<-t(as.matrix(caddis.asin.grouped))

#perform PCA and associated tests
caddis.pca<-prcomp(tcaddis.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(caddis.pca)
#see which eigenvalues are significant
screeplot(caddis.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(caddis.pca,tcaddis.asin.grouped,dim=7,cutoff=0.2)
#sample scores
sample.scores<-caddis.pca$x[,1:7]

#create matrix of other covariates organized so that they can be used as cex values in the plot
caddis.indexed.covariates<-matrix(NA, nrow=length(row.names(tcaddis.asin.grouped)), ncol=45)
for (i in 1:length(caddis.indexed.covariates[1,]))
{
  for (j in 1:length(row.names(tcaddis.asin.grouped)))
  {
    #read lake_ID(but not class) off alldata and match it to lake ID in the indexed covariates
    caddis.indexed.covariates[j,i]<-alldata[1:11,i][substr(alldata$lake_ID_andClass[1:11],1,2)==(as.numeric(substr(row.names(tcaddis.asin.grouped),8,9))[j])]
  }
}
caddis.indexed.covariates<-as.data.frame(caddis.indexed.covariates[-c(1,5),])
colnames(caddis.indexed.covariates)<-colnames(alldata)
caddis.indexed.covariates[,1]<-c("C", "Milk", "Mirror", "Morgenroth", "No_Name", "Y015", "Y025", "Z", "Clear")

#plot PC 1 and 2
#temp <- tcaddis.asin.grouped
caddis.ID<-c('C','C','k','m','M','M','N','y','Y','Z','c')
#IDcorrespond<-c('C','C','Milk','mirror','Morgenroth','Morgenroth','No_Name','Y015','Y025','Z','Clear')
#first and second principal components (showing watershed area:lake area ratio)
caddis.fig<-ordiplot(caddis.pca, choices=c(1,2), type="none", xlim=c(-4,3), ylim=c(-4,4), main="Caddis PC 1 and 2, sized by watershed:lake area ratio")
points(caddis.fig, "sites", pch=20,
       col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)), cex=rescale(caddis.indexed.covariates$shed_lake_areaRatio[c(1,1,2,3,4,4,5,6,7,8,9)], c(2,15))) #rescales all variables to range: 2-15
points(caddis.fig, "sites", pch=caddis.ID, col='yellow')
# points(caddis.fig, "sites", pch=3, bg="transparent", col="gray", cex=rescale(caddis.indexed.covariates$lake_area, c(2,15)))
# points(caddis.fig, "sites", pch=21, bg="transparent", col="gray", cex=rescale(caddis.indexed.covariates$shed_lake_areaRatio, c(2,15)))
arrows(0,0,caddis.pca$rotation[,1]*5, caddis.pca$rotation[,2]*5, col="purple")
text(caddis.pca$rotation[,1]*5.8, caddis.pca$rotation[,2]*5.8, row.names(caddis.pca$rotation), col="purple")
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                            "fish", "POM", "periphyton", "filamentous", "moss",
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4, 11, 7, 9, 14))
legend("topleft", legend=c("sm low forested", "med mid forested", "lg low forested", "high unforested"),
       fill=c(1, 3, 4, 2))
# legend("bottom", legend=c("gray cross = lake area", "gray circle = watershed:area"))

#plot PC 1 and 3
#first and second principal components (showing watershed area:lake area ratio)
caddis.fig<-ordiplot(caddis.pca, choices=c(1,3), type="none", xlim=c(-4,3), ylim=c(-4,5), main="Caddis PC 1 and 3, sized by watershed:lake area ratio")
points(caddis.fig, "sites", pch=20,
       col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)), cex=rescale(caddis.indexed.covariates$shed_lake_areaRatio[c(1,1,2,3,4,4,5,6,7,8,9)], c(2,15))) #rescales all variables to range: 2-15
points(caddis.fig, "sites", pch=ID, col='yellow')
arrows(0,0,caddis.pca$rotation[,1]*5, caddis.pca$rotation[,3]*5, col="purple")
text(caddis.pca$rotation[,1]*5.8, caddis.pca$rotation[,3]*5.8, row.names(caddis.pca$rotation), col="purple")
legend("bottomleft", legend=c("small forested", "med forested", "large forested", "unforested"),
       fill=c(1, 3, 4, 2))

#to loop through with all covariates represented as point size:
#NOTE caddisd15n is misrepresented.  Also note issue with CH4
# for (i in 3:length(caddis.indexed.covariates[1,]))
# {
#   caddis.fig<-ordiplot(caddis.pca, choices=c(1,2), type="none", xlim=c(-6,3), ylim=c(-5,4), main=paste("caddis PC1 and 2, size = ", colnames(caddis.indexed.covariates)[i]))
#   points(caddis.fig, "sites", pch=as.numeric(substr(row.names(tcaddis.asin.grouped),6,7)),
#          col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)), cex=rescale(caddis.indexed.covariates[,i], c(2,15))) #rescales all variables to range: 2-15
#   points(caddis.fig, "sites", pch=3, bg="transparent", col="gray", cex=rescale(caddis.indexed.covariates$lake_area, c(2,15)))
#   points(caddis.fig, "sites", pch=21, bg="transparent", col="gray", cex=rescale(caddis.indexed.covariates$shed_lake_areaRatio, c(2,15)))
#   arrows(0,0,caddis.pca$rotation[,1]*5, caddis.pca$rotation[,2]*5, col="purple")
#   text(caddis.pca$rotation[,1]*5.8, caddis.pca$rotation[,2]*5.8, row.names(caddis.pca$rotation), col="purple")
#   legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                              "fish", "POM", "periphyton", "filamentous", "moss",
#                              "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
#                                                               11, 7, 9, 14))
#   legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
#          fill=c(1, 2, 3, 4, 5))
#   legend("bottom", legend=c("gray cross = lake area", "gray circle = watershedA:lakeA"))
# }

#notable covariates: CO2, CH4,

#plot caddis d13c vs caddis d15n (should show inverse relationship).  (then color by lake type) and size by other covariates
#until you find one that follows the trendline in either direction.

for(i in 3:length(caddis.indexed.covariates[1,]))
{
  plot(caddis.indexed.covariates$caddisd13c, caddis.indexed.covariates$caddisd15n, cex=rescale(caddis.indexed.covariates[,i], c(2,15)),
       main=paste(colnames(caddis.indexed.covariates)[i]), col=substr(caddis.indexed.covariates$lake_ID_andClass,3,3), pch=20, ylim = c(1,5),
       xlab = "caddis d13C", ylab = "caddis d15N")
  text(x = caddis.indexed.covariates$caddisd13c, y = caddis.indexed.covariates$caddisd15n + 0.4, labels = caddis.indexed.covariates$lake)
  legend("topright", legend=c("puddles", "high", "milk+clear", "big"),
         fill=c(1, 2, 3, 4))
}


#plotting 3rd and 4th principal components
# caddis.fig<-ordiplot(caddis.pca, choices=c(3,4), type="none", xlim=c(-4,4), ylim=c(-4,4), main="caddis PC3 and 4")
# points(caddis.fig, "sites", pch=as.numeric(substr(row.names(tcaddis.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcaddis.asin.grouped),5,5)))
# arrows(0,0,caddis.pca$rotation[,3]*5, caddis.pca$rotation[,4]*5, col="black")
# text(caddis.pca$rotation[,3]*5.8, caddis.pca$rotation[,4]*5.8, row.names(caddis.pca$rotation))
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                            "fish", "POM", "periphyton", "filamentous", "moss",
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
#                                                             11, 7, 9, 14))
# legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
#        fill=c(1, 2, 3, 4, 5))
#
# #are there obvious differences among the frequency distributions?

par(mar=c(0,5,0,0))
par(mfrow=c(10,1))
for(i in 1:10)
{
  if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==1)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="black", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==2)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="red", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==3)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="green", xaxt="n", ylim=c(0,25))
    axis(side=1, at=c(4,6,10,12,16,20,22,26,28,30,33,34,35,38,41,44,46,50,52,53,59),
         labels=dist.data$FA_name[c(4,6,10,12,16,20,22,26,28,30,33,34,35,38,41,44,46,50,52,53,59)], las=2)
  }
  else if (substr(colnames(dist.data[caddis.cols[i]]),5,5)==4)
  {
    plot(dist.data[,caddis.cols[i]], type="l", ylab=colnames(dist.data[caddis.cols[i]]), xlab="", col="blue", xaxt="n", ylim=c(0,25))
  }
}
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))

#### plotting covariates with PC1
for(i in 3:length(caddis.indexed.covariates[1,]))
{
  plot(sample.scores[,1], caddis.indexed.covariates[,i], xlim=c(-4,4),
       main=paste(colnames(caddis.indexed.covariates)[i]), col=substr(caddis.indexed.covariates$lake_ID_andClass,3,3), pch=20, cex=3,
       xlab = "PC1", ylab = paste(colnames(caddis.indexed.covariates)[i]))
  text(x = sample.scores[,1]+0.5, y = caddis.indexed.covariates[,i], labels = caddis.indexed.covariates$lake)
  #   legend("topright", legend=c("puddles", "high", "milk+clear", "big"),
  #          fill=c(1, 2, 3, 4))
}

for(i in 3:length(caddis.indexed.covariates[1,]))
{
  plot(caddis.indexed.covariates$color_chlA, caddis.indexed.covariates[,i],
       main=paste("color:chl-a X ", colnames(caddis.indexed.covariates)[i]), col=substr(caddis.indexed.covariates$lake_ID_andClass,3,3), pch=20, cex=3,
       xlab = "color:chlorophyll-a", ylab = paste(colnames(caddis.indexed.covariates)[i]))
  text(x = caddis.indexed.covariates$color_chlA + 0.008, y = caddis.indexed.covariates[,i], labels = caddis.indexed.covariates$lake)
}

#test for correlations among the covariates that seem to diverge between the ponds and the other lakes when plotted against PC1
cor(caddis.indexed.covariates[,c(4,6,11,14,24,29,34)])
scatterplotMatrix(caddis.indexed.covariates[,c(4,6,11,14,24,29,34)], reg.line=F, smooth=F, spread=F, diagonal='density', id.n=0, span=0.5)
#see log for 9/22 to see how I resolved this

#### plotting covariates with PC1 JUST FOR THE NON-PUDDLES
for(i in 3:length(caddis.indexed.covariates[1,]))
{
  plot(sample.scores[,1][-c(1,2,3,7)], caddis.indexed.covariates[,i][-c(1,2,3,7)], xlim=c(-4,4),
       main=paste(colnames(caddis.indexed.covariates)[i]), col=substr(caddis.indexed.covariates$lake_ID_andClass,3,3)[-c(1,2,3,7)], pch=20, cex=3,
       xlab = "PC1", ylab = paste(colnames(caddis.indexed.covariates)[i]))
  text(x = sample.scores[,1][-c(1,2,3,7)]+0.5, y = caddis.indexed.covariates[,i][-c(1,2,3,7)], labels = caddis.indexed.covariates$lake[-c(1,2,3,7)])
}

scatterplotMatrix(caddis.indexed.covariates[-c(1,2,3,7),c(6,11,15,16,29,34)], reg.line=F, smooth=F, spread=F, diagonal='density', id.n=0, span=0.5)
####20 - Cladocera PCA etc. ####

#create FA groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, clad.cols)]
short_SAFA<-colSums(short_SAFA[,6:11], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, clad.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:11], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, clad.cols)]
other_MUFA<-colSums(other_MUFA[,6:11], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, clad.cols)]
LIN<-colSums(LIN[,6:11], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, clad.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:11], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, clad.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:11], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, clad.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:11], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, clad.cols)]
long_SAFA<-colSums(long_SAFA[,6:11], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, clad.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:11], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, clad.cols)]
other_PUFA<-colSums(other_PUFA[,6:11], na.rm=T)

#grouping combination 1 (all)
clad.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
                   ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA")
# clad.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#modify the grouped data for PCA (arcsin sqrt transform and missing data replacement)
for (i in 1:length(clad.grouped[1,]))
{
  for (j in 1:length(clad.grouped[,1]))
  {
    #replace zeros with small values
    if (clad.grouped[j,i]==0)
    {
      clad.grouped[j,i]<-0.000001
    }
    #put data on proportional scale (0-1) and arc-sin square root transform them
    clad.asin.grouped<-clad.grouped
    clad.asin.grouped[j,i]<-asin(sqrt(clad.asin.grouped[j,i]/100))*(2/pi)
  }
}

#format data for PCA
tclad.asin.grouped<-t(as.matrix(clad.asin.grouped))

#perform PCA and associated tests
clad.pca<-prcomp(tclad.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(clad.pca)
#see which eigenvalues are significant
screeplot(clad.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(clad.pca,tclad.asin.grouped,dim=6,cutoff=0.2)
#sample scores
sample.scores<-clad.pca$x[,1:6]

#create matrix of other covariates organized so that they can be used as cex values in the plot
clad.indexed.covariates<-matrix(NA, nrow=length(row.names(tclad.asin.grouped)), ncol=45)
for (i in 1:length(clad.indexed.covariates[1,]))
{
  for (j in 1:length(row.names(tclad.asin.grouped)))
  {
    #read lake_ID(but not class) off alldata and match it to lake ID in the indexed covariates
    clad.indexed.covariates[j,i]<-alldata[1:11,i][substr(alldata$lake_ID_andClass[1:11],1,2)==(as.numeric(substr(row.names(tclad.asin.grouped),8,9))[j])]
  }
}
clad.indexed.covariates<-as.data.frame(clad.indexed.covariates)
colnames(clad.indexed.covariates)<-colnames(alldata)
clad.indexed.covariates[,1]<-c("Clear", "Morgenroth", "O", "L", "C", "Z")

#plot PCA
#first and second principal components (showing watershed area:lake area ratio)
clad.ID <- c('c','M','O','L','C','Z')
clad.fig<-ordiplot(clad.pca, choices=c(1,2), type="none", xlim=c(-4,4), ylim=c(-3,3), main="Cladocera PCA, sized by watershed:lake area ratio")
points(clad.fig, "sites", pch=20,
       col=as.numeric(substr(row.names(tclad.asin.grouped),5,5)), cex=rescale(clad.indexed.covariates$shed_lake_areaRatio, c(2,15))) #rescales all variables to range: 2-15
points(clad.fig, "sites", pch=clad.ID, col='yellow')
# points(clad.fig, "sites", pch=3, bg="transparent", col="gray", cex=rescale(clad.indexed.covariates$lake_area, c(2,15)))
# points(clad.fig, "sites", pch=21, bg="transparent", col="gray", cex=rescale(clad.indexed.covariates$shed_lake_areaRatio, c(2,15)))
arrows(0,0,clad.pca$rotation[,1]*5, clad.pca$rotation[,2]*5, col="purple")
text(clad.pca$rotation[,1]*5.8, clad.pca$rotation[,2]*5.8, row.names(clad.pca$rotation), col="purple")
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                            "fish", "POM", "periphyton", "filamentous", "moss",
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4, 11, 7, 9, 14))
legend("bottomleft", legend=c("small forested", "med forested", "large forested", "unforested"),
       fill=c(1, 3, 4, 2))
# legend("bottom", legend=c("gray cross = lake area", "gray circle = watershed:area"))


#to loop through with all covariates represented as point size:
#NOTE caddisd15n is misrepresented.  Also note issue with CH4
# for (i in 3:length(clad.indexed.covariates[1,]))
# {
#   clad.fig<-ordiplot(clad.pca, choices=c(1,2), type="none", xlim=c(-6,3), ylim=c(-5,4), main=paste("clad PC1 and 2, size = ", colnames(clad.indexed.covariates)[i]))
#   points(clad.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)),
#          col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)), cex=rescale(clad.indexed.covariates[,i], c(2,15))) #rescales all variables to range: 2-15
#   points(clad.fig, "sites", pch=3, bg="transparent", col="gray", cex=rescale(clad.indexed.covariates$lake_area, c(2,15)))
#   points(clad.fig, "sites", pch=21, bg="transparent", col="gray", cex=rescale(clad.indexed.covariates$shed_lake_areaRatio, c(2,15)))
#   arrows(0,0,clad.pca$rotation[,1]*5, clad.pca$rotation[,2]*5, col="purple")
#   text(clad.pca$rotation[,1]*5.8, clad.pca$rotation[,2]*5.8, row.names(clad.pca$rotation), col="purple")
#   legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                              "fish", "POM", "periphyton", "filamentous", "moss",
#                              "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
#                                                               11, 7, 9, 14))
#   legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
#          fill=c(1, 2, 3, 4, 5))
#   legend("bottom", legend=c("gray cross = lake area", "gray circle = watershedA:lakeA"))
# }

#notable covariates: CO2, CH4,

#plot caddis d13c vs caddis d15n (should show inverse relationship).  (then color by lake type) and size by other covariates
#until you find one that follows the trendline in either direction.

for(i in 3:length(clad.indexed.covariates[1,]))
{
  plot(clad.indexed.covariates$calanoidd13c, clad.indexed.covariates$calanoidd15n, cex=rescale(clad.indexed.covariates[,i], c(2,15)),
       main=paste(colnames(clad.indexed.covariates)[i]), col=substr(clad.indexed.covariates$lake_ID_andClass,3,3), pch=20, ylim = c(1,5),
       xlab = "Calanoid d13C", ylab = "Calanoid d15N")
  text(x = clad.indexed.covariates$calanoidd13c, y = clad.indexed.covariates$calanoidd15n + 0.4, labels = clad.indexed.covariates$lake)
  legend("topright", legend=c("puddles", "high", "milk+clear", "big"),
         fill=c(1, 2, 3, 4))
}


#plotting 3rd and 4th principal components
# clad.fig<-ordiplot(clad.pca, choices=c(3,4), type="none", xlim=c(-4,4), ylim=c(-4,4), main="clad PC3 and 4")
# points(clad.fig, "sites", pch=as.numeric(substr(row.names(tcal.asin.grouped),6,7)), col=as.numeric(substr(row.names(tcal.asin.grouped),5,5)))
# arrows(0,0,clad.pca$rotation[,3]*5, clad.pca$rotation[,4]*5, col="black")
# text(clad.pca$rotation[,3]*5.8, clad.pca$rotation[,4]*5.8, row.names(clad.pca$rotation))
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
#                            "fish", "POM", "periphyton", "filamentous", "moss",
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
#                                                             11, 7, 9, 14))
# legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
#        fill=c(1, 2, 3, 4, 5))

#are there obvious differences among the frequency distributions?

par(mar=c(0,5,0,0))
par(mfrow=c(10,1))
for(i in 1:10)
{
  if (substr(colnames(dist.data[clad.cols[i]]),5,5)==1)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="black", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[clad.cols[i]]),5,5)==2)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="red", xaxt="n", ylim=c(0,25))
  }
  else if (substr(colnames(dist.data[clad.cols[i]]),5,5)==3)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="green", xaxt="n", ylim=c(0,25))
    axis(side=1, at=c(4,8,9,10,12,14,16,17,18,22,26,28,29,33,34,35,38,41,44,50,52,58),
         labels=dist.data$FA_name[c(4,8,9,10,12,14,16,17,18,22,26,28,29,33,34,35,38,41,44,50,52,58)], las=2)

  }
  else if (substr(colnames(dist.data[clad.cols[i]]),5,5)==4)
  {
    plot(dist.data[,clad.cols[i]], type="l", ylab=colnames(dist.data[clad.cols[i]]), xlab="", col="blue", xaxt="n", ylim=c(0,25))
  }
}
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))

#### plotting covariates with PC1
for(i in 3:length(clad.indexed.covariates[1,]))
{
  plot(sample.scores[,1], clad.indexed.covariates[,i], xlim=c(-4,4),
       main=paste(colnames(clad.indexed.covariates)[i]), col=substr(clad.indexed.covariates$lake_ID_andClass,3,3), pch=20, cex=3,
       xlab = "PC1", ylab = paste(colnames(clad.indexed.covariates)[i]))
  text(x = sample.scores[,1]+0.5, y = clad.indexed.covariates[,i], labels = clad.indexed.covariates$lake)
  #   legend("topright", legend=c("puddles", "high", "milk+clear", "big"),
  #          fill=c(1, 2, 3, 4))
}

for(i in 3:length(clad.indexed.covariates[1,]))
{
  plot(clad.indexed.covariates$color_chlA, clad.indexed.covariates[,i],
       main=paste("color:chl-a X ", colnames(clad.indexed.covariates)[i]), col=substr(clad.indexed.covariates$lake_ID_andClass,3,3), pch=20, cex=3,
       xlab = "color:chlorophyll-a", ylab = paste(colnames(clad.indexed.covariates)[i]))
  text(x = clad.indexed.covariates$color_chlA + 0.008, y = clad.indexed.covariates[,i], labels = clad.indexed.covariates$lake)
}

#test for correlations among the covariates that seem to diverge between the ponds and the other lakes when plotted against PC1
cor(clad.indexed.covariates[,c(4,6,11,14,24,29,34)])
scatterplotMatrix(clad.indexed.covariates[,c(4,6,11,14,24,29,34)], reg.line=F, smooth=F, spread=F, diagonal='density', id.n=0, span=0.5)
#see log for 9/22 to see how I resolved this

#### plotting covariates with PC1 JUST FOR THE NON-PUDDLES
for(i in 3:length(clad.indexed.covariates[1,]))
{
  plot(sample.scores[,1][-c(1,2,3,7)], clad.indexed.covariates[,i][-c(1,2,3,7)], xlim=c(-4,4),
       main=paste(colnames(clad.indexed.covariates)[i]), col=substr(clad.indexed.covariates$lake_ID_andClass,3,3)[-c(1,2,3,7)], pch=20, cex=3,
       xlab = "PC1", ylab = paste(colnames(clad.indexed.covariates)[i]))
  text(x = sample.scores[,1][-c(1,2,3,7)]+0.5, y = clad.indexed.covariates[,i][-c(1,2,3,7)], labels = clad.indexed.covariates$lake[-c(1,2,3,7)])
}

scatterplotMatrix(clad.indexed.covariates[-c(1,2,3,7),c(6,11,15,16,29,34)], reg.line=F, smooth=F, spread=F, diagonal='density', id.n=0, span=0.5)

####21 - Producer PCAs ####

#create FA groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:11], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:11], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:11], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:11], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:11], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:11], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:11], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:11], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:11], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:11], na.rm=T)

#grouping combination 1 (all)
prod.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
                   ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA")
# cal.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#modify the grouped data for PCA (arcsin sqrt transform and missing data replacement)
for (i in 1:length(prod.grouped[1,]))
{
  for (j in 1:length(prod.grouped[,1]))
  {
    #replace zeros with small values
    if (prod.grouped[j,i]==0)
    {
      prod.grouped[j,i]<-0.000001
    }
    #put data on proportional scale (0-1) and arc-sin square root transform them
    prod.asin.grouped<-prod.grouped
    prod.asin.grouped[j,i]<-asin(sqrt(prod.asin.grouped[j,i]/100))*(2/pi)
  }
}

#format data for PCA
tprod.asin.grouped<-t(as.matrix(prod.asin.grouped))

#perform PCA and associated tests
prod.pca<-prcomp(tprod.asin.grouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(prod.pca)
#see which eigenvalues are significant
screeplot(prod.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(prod.pca,tprod.asin.grouped,dim=7,cutoff=0.2)
#sample.scores<-prod.pca$x[,1:7]

#plot PCA
#first and second principal components (counterintuitive groupings)
prod.fig<-ordiplot(prod.pca, choices=c(1,2), type="none", xlim=c(-6,3), ylim=c(-5,4), main="prod PC1 and 2")
points(prod.fig, "sites", pch=as.numeric(substr(row.names(tprod.asin.grouped),6,7)),
       col=as.numeric(substr(row.names(tprod.asin.grouped),5,5)))
arrows(0,0,prod.pca$rotation[,1]*5, prod.pca$rotation[,2]*5, col="purple")
text(prod.pca$rotation[,1]*5.8, prod.pca$rotation[,2]*5.8, row.names(prod.pca$rotation), col="purple")
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

#first and second principal components (PC3 looks more like what I'd expect.)
prod.fig<-ordiplot(prod.pca, choices=c(2,3), type="none", xlim=c(-6,3), ylim=c(-5,4), main="prod PC2 and 3")
points(prod.fig, "sites", pch=as.numeric(substr(row.names(tprod.asin.grouped),6,7)),
       col=as.numeric(substr(row.names(tprod.asin.grouped),5,5)))
arrows(0,0,prod.pca$rotation[,2]*5, prod.pca$rotation[,3]*5, col="purple")
text(prod.pca$rotation[,2]*5.8, prod.pca$rotation[,3]*5.8, row.names(prod.pca$rotation), col="purple")
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly",
                           "fish", "POM", "periphyton", "filamentous", "moss",
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

########22 - Hierarchical cluster analysis functions####
library(cluster)
library(vegan)
library(pvclust)
library(NbClust)

hclus.table <-
  function(x){

    z1<-x$dist.method
    z2<-x$method
    z3<-seq(length(x$labels)-1,1)
    z3<-as.data.frame(cbind(z3,x$merge,x$height))
    colnames(z3)<-c('no. clusters','entity','entity','distance')
    z<-list(z1,z2,z3)
    names(z)<-c('dist.method','method','cluster.table')
    return(z)
  }

data.trans <-
  function(x,method,var='',exp=1,outfile='',
           plot=TRUE,save.plot=FALSE,col.hist='blue',col.line='black',
           las=1,lab=c(5,5,4),...){

    if(plot==TRUE){
      old.par<-par(no.readonly=TRUE)
    }

    if(!var==''){
      y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
      y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
      t1<-y1 #copy to work file for transformations
    }
    else{
      y1<-x #original variables
      t1<-x #copy to work file for transformations
    }
  }

hclus.cophenetic <-
  function(d,hclus,fit='lm',...){

    old.par<-par(no.readonly=TRUE)
    d.coph<-cophenetic(hclus)
    r<-round(cor(d,d.coph),2)
    plot(d,d.coph,xlab='Observed dissimilarity',
         ylab='Cophenetic dissimilarity',
         main=paste('Cophenetic Correlation ',
                    '(',hclus$dist.method,', ',hclus$method,')',sep=''),...)
    text(max(d),min(d.coph),paste('Cophenetic correlation = ',r,sep=''),col='red',pos=2)
    #	title(sub=paste('Cophenetic correlation = ',r,sep=''),col.sub='red',adj=0)
    if(fit=='lm'){
      abline(lm(d.coph~d),col='blue')
    }
    else if(fit=='rlm'){
      abline(rlm(d.coph~d),col='blue')
    }
    else if(fit=='qls'){
      abline(lqs(d.coph~d),col='blue')
    }

    par(old.par)
    return(r)
  }

hclus.scree <-
  function(x,...){

    old.par<-par(no.readonly=TRUE)
    z1<-seq(length(x$height),1)
    z<-as.data.frame(cbind(z1,sort(x$height)))
    plot(z[,1],z[,2],type='o',lwd=1.5,pch=19,col='blue',
         ylab='Dissimilarity',xlab='Number of Clusters',
         main=paste('Scree Plot of Hierarchical Clustering ',
                    '(',x$dist.method,', ',x$method,')',sep=''),...)
    par(old.par)
  }
####23 - HCA1 (all species)####
short_SAFA<-data[data$short_SAFA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:42], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:42], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:42], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:42], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:42], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:42], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:42], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:42], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:42], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, cal.cols, clad.cols, caddis.cols, chaob.cols, fish.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:42], na.rm=T)

#grouping 1
hca1.grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
                   ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA")
# hca1.grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#make more descriptive column names
colnames(hca1.grouped)<-c("cal C", "cal O", "cal Z", "cal Mirror", "cal Y025", "cal Y015", "cal L", "cal Clear",
                          "cal Morg", "cal Milk", "clad Clear", "clad Morg", "clad O", "clad L", "clad C", "clad Z",
                          "caddis C", "caddis C2", "caddis Milk", "caddis Mirror", "caddis Morg", "caddis Morg2",
                          "caddis NoName", "caddis Y015", "caddis Y025", "caddis Z", "caddis Clear", "chaob Clear",
                          "chaob Milk", "fish Mirror", "fish Mirror2", "POM", "peri", "filamentous", "moss", "soil",
                          "plant")

#modify the data
for (i in 1:length(hca1.grouped[1,]))
{
  for (j in 1:length(hca1.grouped[,1]))
  {
    #replace zeros with small values
    if (hca1.grouped[j,i]==0)
    {
      hca1.grouped[j,i]<-0.000001
    }
    #put data on proportional scale (0-1) and arc-sin square root transform them
    hca1.asin.grouped<-hca1.grouped
    hca1.asin.grouped[j,i]<-asin(sqrt(hca1.asin.grouped[j,i]/100))*(2/pi)
  }
}

#the HCA - seems like the order of these steps is wonky.  should be 1. find correct number of clusters with NbClust,
#2. plot, 3. highlight that number of clusters
#create dissimilarity matrix
hc1.distmat<-vegdist(t(hca1.asin.grouped), method="euclidean")
#perform clustering (many other methods can be used: ward, single, average, mcquitty, median, centroid); this method is agglomerative - can use divisive ones too
hc1.clust<-hclust(hc1.distmat, method="complete", members=NULL)
#summary results
  #hclus.table(hc1.clust)
#dendrogram
plot(hc1.clust, main="complete linkage dendrogram", xlab="species", ylab="euclidean distance", hang=-1)
#add rectangles around clusters by specifying desired number (how to choose the number of clusters?)
rect.hclust(hc1.clust, k=2)
#add rectangles by specifying a euclidean distance
  # rect.hclust(hc1.clust, h=22, border="green")
#cut the tree into clusters (when does this get used?)
  # hc1.classes<-cutree(hc1.clust, k=6)
#find the agglomerative coefficient - if obs quickly aggregate and only fully unite at much greater distance, coef = 1.  if they take forever to aggregate, coef = 0
coef.hclust(hc1.clust)
#cophenetic correlation coefficient: if high, dendrogram is appropriate summary of "some data." if not, it's only a description of the output of the clustering algorithm.
cor(hc1.distmat, cophenetic(hc1.clust))
# hclus.cophenetic(hc1.distmat, hc1.clust)
#if this screeplot has an "elbow", you should cut the tree at the number of clusters that corresponds to the elbow.  if it doesn't, try another clustering method.
#reading from right to left (says to do this in notes), take the first number you reach before the elbow occurs (my assumption)
hclus.scree(hc1.clust)
#test to make sure the groupings are stable through minor perturbations of the data (removal of one point at a time); should be transposed relative to the original data arrangement
hc1.bootstrap<-pvclust(hca1.asin.grouped, method.hclust="complete", method.dist="euclidean", nboot=200)
#plot to highlight significant clusters (red and green are different measurements of significance. use red)
plot(hc1.bootstrap, hang=-1)
pvrect(hc1.bootstrap, alpha=0.99)
#now you have to find the "optimal" number of clusters.  This function incorporates 30 different indices.  Majority rule is fine to go by, but
#CH, Duda, Cindex, Gamma, and Beale were shown to perform best by milligan and cooper 1985.
NbClust(t(hca1.asin.grouped), distance="euclidean", min.nc=2, max.nc=10, method="complete")
par(mfrow=c(1,1))

############two groups recommended.  peri, POM separate from plant, filamentous, soil.  only 1 fish, most cladocera, and chaoborus group with plant

#here it is again using method=average (slightly better cophanetic cor, otherwise same) -BEST
hc1.clust2<-hclust(hc1.distmat, method="average", members=NULL)
  coef.hclust(hc1.clust2)
  cor(hc1.distmat, cophenetic(hc1.clust2))
hclus.scree(hc1.clust2)
NbClust(t(hca1.asin.grouped), distance="euclidean", min.nc=2, max.nc=10, method="average")
  par(mfrow=c(1,1))
hc1.bootstrap2<-pvclust(hca1.asin.grouped, method.hclust="average", method.dist="euclidean", nboot=100)
  plot(hc1.bootstrap2, hang=-1)
  pvrect(hc1.bootstrap2, alpha=0.99)
plot(hc1.clust2, main="average linkage dendrogram", xlab="species", ylab="euclidean distance", hang=-1)
rect.hclust(hc1.clust, k=2)

#here it is again using method=centroid (meh)
hc1.clust3<-hclust(hc1.distmat, method="centroid", members=NULL)
coef.hclust(hc1.clust3)
cor(hc1.distmat, cophenetic(hc1.clust3))
hclus.scree(hc1.clust3)
NbClust(t(hca1.asin.grouped), distance="euclidean", min.nc=2, max.nc=10, method="centroid")
par(mfrow=c(1,1))
hc1.bootstrap3<-pvclust(hca1.asin.grouped, method.hclust="centroid", method.dist="euclidean", nboot=100)
plot(hc1.bootstrap3, hang=-1)
pvrect(hc1.bootstrap3, alpha=0.99)
plot(hc1.clust3, main="centroid linkage dendrogram", xlab="species", ylab="euclidean distance", hang=-1)
rect.hclust(hc1.clust, k=2)

#here it is again using method=single (weak)
hc1.clust4<-hclust(hc1.distmat, method="single", members=NULL)
coef.hclust(hc1.clust4)
cor(hc1.distmat, cophenetic(hc1.clust4))
hclus.scree(hc1.clust4)
NbClust(t(hca1.asin.grouped), distance="euclidean", min.nc=2, max.nc=10, method="single")
par(mfrow=c(1,1))
hc1.bootstrap4<-pvclust(hca1.asin.grouped, method.hclust="single", method.dist="euclidean", nboot=100)
plot(hc1.bootstrap4, hang=-1)
pvrect(hc1.bootstrap4, alpha=0.99)
plot(hc1.clust4, main="single linkage dendrogram", xlab="species", ylab="euclidean distance", hang=-1)
rect.hclust(hc1.clust, k=2)
####24 - cal.HCA####
