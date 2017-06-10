setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")
library(vegan)

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

FA_PCAer<-function(csv.name, ordi.monte=F)
{
  #read in the FAs and samples you want to test
  csv<-read.csv(csv.name)
  
  #find number and indices of columns containing Areal data
  Nsamp<-length(grep("Area",colnames(csv)))
  samp.ind<-grep("Area",colnames(csv))
  proportional.area<-matrix(data=NA, nrow=length(csv[,1]), ncol=Nsamp)
  
  #create matrix of proportional area
  for(i in 1:Nsamp)
  {
    col<-csv[,samp.ind[i]]
    total.area<-sum(na.omit(col))
    for(j in 1:length(col))
    {
      proportional.area[j,i]<-(col[j]*100)/total.area #find out why i and j are being created as objects
    }
    NAs<-which(is.na(proportional.area[,i]))
    proportional.area[,i][NAs]<-0.0000001
  }
  proportional.area<-cbind(csv[,2:5],proportional.area)
  colnames(proportional.area)[5:length(proportional.area[1,])]<-gsub("Reten_", "", colnames(csv)[grep("Reten", colnames(csv))])
  row.names(proportional.area)<-csv[,1]
  
  #the PCA...
  transposed<-t(as.matrix(proportional.area[,-(1:4)]))
  pca<-prcomp(transposed,scale=T, scores=T)
  #determine eigenvalues
  eigenvalues<-pca.eigenval(pca)
  #see which eigenvalues are significant
  screeplot(pca, bstick=T)
  #do that a different way
  if(ordi.monte==T)
  {
    ordimonte<-o
  }
  #see loadings.  square these to get percentage of variance in each original variable
  #accounted for by each principal component
  structure<-pca.structure(pca,transposed,dim=7,cutoff=0.2)
  #sample.scores<-pca$x[,1:7]
  #first plotting method
  biplot(pca)
  #second plotting method
  ordiplot(pca,choices=c(1,2), type="text", display="sites")
  arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="red")
  text(pca$rotation[,1]*5.2, pca$rotation[,2]*5.2, row.names(pca$rotation))
  
  if(ordi.monte==T)
  {
    details<-list(eigenvalues[,1:7], ordimonte, structure)
  }
  else
  {
    details<-list(eigenvalues[,1:7], structure)
  }
  return(details)
}

FA_PCAer("1_allConsumers_FAoverPointFive.csv")
FA_PCAer("2_allConsumers_FAoverPointSeven.csv")
FA_PCAer("4_caddisonly_29FAs.csv")
FA_PCAer("5_caddisonly_20FAs.csv")
FA_PCAer("6_calonly_43FAs.csv")
FA_PCAer("7_cladonly_28FAs.csv")



