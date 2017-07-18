
#FA plot is in conference folder, with ASLO stuff

#so far all attempts at fixing parameters seem to yield untrustworthy results
#truncated wH may be useful for showing how hard measuring allochthony becomes at its lower extreme
#next, try different truncations and also different fixed values


rm(list=ls()); cat('\014')

# funcs, packages, defaults ####
setwd('C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/mixing_models2/sensitivity') #this is needed for wH plots (and everything else, probs)
# setwd('C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/sensitivity_analysis') #this one is needed for (+/- 2 SD)
data_fetcher <- function(){

    conds <- c('allvary','fixedw','fixedeH','fixedCNA')
    for(i in conds){
        for(j in c('H','C','N')){
            filename <- paste(i,'_',j,'.rda', sep='')
            if(!filename %in% dir()) next
            load(filename)
        }
    }

    allvary_CN <<- allvary_C[[3]][,8]
    fixedwC_CN <<- fixedw_C[[3]][,8]
    fixedCN_CN <<- fixedCNA_C[[3]][,8]
    fixedeH_CN <<- allvary_C[[3]][,8]
    fixedwH_CN <<- allvary_C[[3]][,8]
    fixedwN_CN <<- fixedw_N[[3]][,8]

    allvary_wC <<- allvary_C[[3]][,2]
    fixedwC_wC <<- fixedw_C[[3]][,2]
    fixedCN_wC <<- fixedCNA_C[[3]][,2]
    fixedeH_wC <<- allvary_C[[3]][,2]
    fixedwH_wC <<- allvary_C[[3]][,2]
    fixedwN_wC <<- allvary_C[[3]][,2]

    allvary_wN <<- allvary_N[[3]][,2]
    fixedwC_wN <<- allvary_N[[3]][,2]
    fixedCN_wN <<- fixedCNA_N[[3]][,2]
    fixedeH_wN <<- allvary_N[[3]][,2]
    fixedwH_wN <<- allvary_N[[3]][,2]
    fixedwN_wN <<- fixedw_N[[3]][,2]

    allvary_wH <<- allvary_H[[3]][,5]
    fixedwC_wH <<- allvary_H[[3]][,5]
    fixedCN_wH <<- allvary_H[[3]][,5]
    fixedeH_wH <<- fixedeH_H[[3]][,5]
    fixedwH_wH <<- fixedw_H[[3]][,5]
    fixedwN_wH <<- allvary_H[[3]][,5]

    allvary_eH <<- allvary_H[[3]][,2]
    fixedwC_eH <<- allvary_H[[3]][,2]
    fixedCN_eH <<- allvary_H[[3]][,2]
    fixedeH_eH <<- fixedeH_H[[3]][,2]
    fixedwH_eH <<- fixedw_H[[3]][,2]
    fixedwN_eH <<- allvary_H[[3]][,2]
}
data_fetcher()

# posteriors for estimated parameters ####
par(mfrow=c(6,5), mar=c(0,0,0,0), oma=c(5,.5,3,.5))
titles <- list('None', expression(paste(Delta^13,'C',sep='')),
               expression(paste(Delta^15,'N',sep='')),
               expression(paste(omega,sep='')), 'C:N',
               expression(paste(epsilon,'H',sep='')))

conds <- c('allvary','fixedwC','fixedwN','fixedwH','fixedCN','fixedeH')
for(i in conds){

    #fixed wC
    if(!substr(i,nchar(i)-1,nchar(i)) == 'wC'){
        plot(density(eval(parse(text=paste(i,'wC',sep='_'))),na.rm=TRUE),
             ann=FALSE, xlim=c(-4.5,5.5), ylim=c(0,0.37), axes=FALSE, xaxs='i', yaxs='i')
        # axis(side=2,las=2, at=c(0,0.1,0.2,0.3), labels=c('0.0','0.1','0.2','0.3'))
        if(i=='fixedeH'){
            axis(side=1)
        }
        if(i=='allvary'){
            mtext(titles[[2]], side=3)
        }
        # mtext(titles[[which(conds==i)]],side=2)
    } else {
        plot(1,1,type='n',xlim=c(-4.5,6), ylim=c(0,0.4), axes=FALSE, xaxs='i', yaxs='i')
        abline(v=0.39, lty=2)
        abline(h=0, lty=1, col='gray')
        # mtext(titles[[which(conds==i)]],side=2)
    }
    #fixed wN
    if(!substr(i,nchar(i)-1,nchar(i)) == 'wN'){
        plot(density(eval(parse(text=paste(i,'wN',sep='_'))),na.rm=TRUE),
             ann=FALSE, xlim=c(-0.5,7.5), ylim=c(0,0.42), axes=FALSE, xaxs='i', yaxs='i')
        if(i=='fixedeH'){
            axis(side=1)
        }
        if(i=='allvary'){
            mtext(titles[[3]], side=3)
        }
    } else {
        plot(1,1,type='n',xlim=c(-0.5,7.5), ylim=c(0,0.42), axes=FALSE, xaxs='i', yaxs='i')
        abline(v=3.4, lty=2)
        abline(h=0, lty=1, col='gray')
    }
    #fixedwH
    if(!substr(i,nchar(i)-1,nchar(i)) == 'wH'){
        plot(density(eval(parse(text=paste(i,'wH',sep='_'))),na.rm=TRUE),
             ann=FALSE, xlim=c(-1.1,1.5), ylim=c(0,2.4), axes=FALSE, xaxs='i', yaxs='i')
        if(i=='fixedeH'){
            axis(side=1, at=c(-1,0,1), labels=c('-1', '0', '1'))
        }
        if(i=='allvary'){
            mtext(titles[[4]], side=3)
        }
    } else {
        plot(1,1,type='n',xlim=c(-1.1,1.5), ylim=c(0,2.4), axes=FALSE, xaxs='i', yaxs='i')
        abline(v=0.28, lty=2)
        abline(h=0, lty=1, col='gray')
    }
    #fixedCN
    if(!substr(i,nchar(i)-1,nchar(i)) == 'CN'){
        plot(density(eval(parse(text=paste(i,'CN',sep='_'))),na.rm=TRUE),
             ann=FALSE, xlim=c(-40,50), ylim=c(0,0.2), axes=FALSE, xaxs='i', yaxs='i')
        if(i=='fixedeH'){
            axis(side=1)
        }
        if(i=='allvary'){
            mtext(titles[[5]], side=3)
        }
    } else {
        plot(1,1,type='n',xlim=c(-40,50), ylim=c(0,0.2), axes=FALSE, xaxs='i', yaxs='i')
        abline(v=4.74, lty=2)
        abline(h=0, lty=1, col='gray')
    }
    #fixedeH
    if(!substr(i,nchar(i)-1,nchar(i)) == 'eH'){
        plot(density(eval(parse(text=paste(i,'eH',sep='_'))),na.rm=TRUE),
             ann=FALSE, xlim=c(-230,-70), ylim=c(0,0.03), axes=FALSE, xaxs='i', yaxs='i')
        if(i=='allvary'){
            mtext(titles[[6]], side=3)
        }
    } else {
        plot(1,1,type='n',xlim=c(-230,-70), ylim=c(0,0.03), axes=FALSE, xaxs='i', yaxs='i')
        abline(v=-150, lty=2)
        abline(h=0, lty=1, col='gray')
        if(i=='fixedeH'){
            axis(side=1)
        }
    }
    # mtext('Prior Parameter Density', side=2, line=2.5, outer=TRUE, cex=1.3)
    mtext('Posterior Parameter Density', side=1, line=3.5, outer=TRUE, cex=1.3)
    # mtext('', side=3, line=2.5, outer=TRUE, cex=1.3)
}

# grouped posterior plots (obsolete) ####
library(abind)
library(RColorBrewer)
library(plyr)
cols <- brewer.pal(3, 'Dark2')
reassembler_gen3 <- function(mcmc_list){ #return 'MLEs' and densities as arrays
    nchains <- length(mcmc_list)
    chain_length <- length(mcmc_list[[1]][,1])
    param_list <- gsub(" *\\[.*?\\] *", '', colnames(mcmc_list[[1]]))
    all_ps <- which(param_list == 'p' | param_list == 'mean_p_sample' | param_list == 'mean_p_predator')
    num_ps <- length(all_ps)
    first_p <- all_ps[1]
    total_vals <- num_ps * nchains
    sum_chain_lengths <- chain_length * nchains
    indices <- seq(1, npreds*npreytypes, npreytypes)
    MLE_arr <- array(data=NA, dim=c(npreytypes,npreds))
    dens_arr <- array(data=NA, dim=c(npreytypes,npreds,sum_chain_lengths))

    #create vector of MLEs across chains
    temp <- matrix(data=NA, nrow=chain_length, ncol=nchains)
    MLE_vec <- vector('numeric', length=num_ps)

    for(i in 1:num_ps){
        for(j in 1:nchains){
            temp[,j] <- mcmc_list[[j]][,i+first_p-1][1:chain_length]
        }
        chain <- as.vector(temp)
        MLE_vec[i] <- density(chain, n=10000)$x[which.max(density(chain, n=10000)$y)]
    }

    #assemble MLEs into array
    MLE_arr[1,1:31] <- MLE_vec[indices]
    MLE_arr[2,1:31] <- MLE_vec[indices+1]
    if(npreytypes == 3){
        MLE_arr[3,1:31] <- MLE_vec[indices+2]
    }

    #make separate array of combined chains
    temp2 <- matrix(data=NA, nrow=chain_length, ncol=nchains)
    temp3 <- matrix(data=NA, nrow=chain_length, ncol=nchains)
    temp4 <- matrix(data=NA, nrow=chain_length, ncol=nchains)

    for(i in 1:length(indices)){
        for(j in 1:nchains){
            temp2[,j] <- mcmc_list[[j]][,(first_p-1+indices[i])][1:length(mcmc_list[[j]][,1])]
            temp3[,j] <- mcmc_list[[j]][,(first_p-1+(indices+1)[i])][1:length(mcmc_list[[j]][,1])]
            if(npreytypes == 3){
                temp4[,j] <- mcmc_list[[j]][,(first_p-1+(indices+2)[i])][1:length(mcmc_list[[j]][,1])]
            }
        }
        dens_arr[1,i,] <- as.vector(temp2)
        dens_arr[2,i,] <- as.vector(temp3)
        if(npreytypes == 3){
            dens_arr[3,i,] <- as.vector(temp4)
        }
    }

    #make third array of combined source proportions across lake types
    lgphyto <- t(as.matrix(as.vector(dens_arr[1,which(dens_lktypes==1),])))
    smphyto <- t(as.matrix(as.vector(dens_arr[1,which(dens_lktypes==2),])))
    hiphyto <- t(as.matrix(as.vector(dens_arr[1,which(dens_lktypes==3),])))
    phy <- t(rbind.fill.matrix(lgphyto, smphyto, hiphyto))
    lgperi <- t(as.matrix(as.vector(dens_arr[2,which(dens_lktypes==1),])))
    smperi <- t(as.matrix(as.vector(dens_arr[2,which(dens_lktypes==2),])))
    hiperi <- t(as.matrix(as.vector(dens_arr[2,which(dens_lktypes==3),])))
    per <- t(rbind.fill.matrix(lgperi, smperi, hiperi))
    lgterr <- t(as.matrix(as.vector(dens_arr[3,which(dens_lktypes==1),])))
    smterr <- t(as.matrix(as.vector(dens_arr[3,which(dens_lktypes==2),])))
    hiterr <- t(as.matrix(as.vector(dens_arr[3,which(dens_lktypes==3),])))
    ter <- t(rbind.fill.matrix(lgterr, smterr, hiterr))
    comb <- abind(phy, per, ter, along=3) #rows are draws; cols are lg, sm, hi; depths are phyto, peri, terr

    #and a final array for laketype * predtype (cal, clad, caddis) * draws (allochth only)



    return(list(MLE_arr, dens_arr, comb))
}
box <- function(dens_chains){

    par(mar=c(4,4,2,9))#, oma=c(0,0,0,0))

    boxplot(dens_chains[npreytypes,23,], xlim=c(0,dim(dens_chains)[2]), at=1, ylim=c(0,1), outline=FALSE,
            pars=list(outpch=16, outcex=0.7, outcol=cols[1], boxcol=cols[1], whiskcol=cols[1],
                      staplecol=cols[1], medcol=cols[1]), axes=FALSE, frame=FALSE, lwd=2)

    loopcols <- cols[dens_lktypes]

    for(i in order[-1]){
        boxplot(dens_chains[npreytypes,i,], xlim=c(0,dim(dens_chains)[2]), at=which(order==i), add=TRUE, ylim=c(0,1),
                pars=list(outpch=16, outcex=0.7, outcol=loopcols[i], boxcol=loopcols[i],
                          whiskcol=loopcols[i], staplecol=loopcols[i], medcol=loopcols[i]),
                axes=FALSE, bty='n', frame=FALSE, lwd=2, outline=FALSE)
    }


    axis(side=1, at=1:31, labels=as.vector(pdfac[order]), lwd.ticks=-1, padj=-1.5)
    mtext('Proportion Allochthony', side=2, line=1, cex=1.2)
    axis(side=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1),
         labels=c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'), las=2, pos=0.5, font=2)
    mtext('Consumer Taxon', side=1, line=1.5, cex=1.2)
    # legend(x=1,y=1.032, legend=c('Montane-Large','Montane-Small','Alpine','Calanoida','Cladocera',
    #        'Trichoptera','Salvelinus','Chaoborus'), bty='n', border=NA, fill=c(cols,rep('transparent',5)),
    #        pch=c('','','','1','2','3','4','5'), cex=1.2)
    # lines(x=c(1.65,7),y=c(.843,.843), lwd=2)
    legend(x=31.5,y=.7, legend=c('Montane-Large','Montane-Small','Alpine','Calanoida','Cladocera',
                                 'Trichoptera','Salvelinus','Chaoborus'), bty='n', border=NA, fill=c(cols,rep('transparent',5)),
           pch=c('','','','1','2','3','4','5'), cex=1.2, xpd=NA)
    lines(x=c(32.2,39.4),y=c(.532,.532), lwd=2, xpd=NA)
}

SI <- read.csv("allvary_cor.csv", header=TRUE)
SI <- SI[c(13,23,21,27,17,24,14,19,18,16,28,30,22,29,25,26,31,1,2,3,4,15,5,6,7,8,9,10,11,12,20),] #converts to old order
SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)
npreds = 31
npreytypes = 3
nlaketypes = 3
dens_lktypes <- c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1)
dens_lknames <- SI[,1]
dens_pdtypes <- SI[,2]
pdfac <- factor(dens_pdtypes)
levels(pdfac) <- c(1,5,2,4,3)
order <- c(2,6,13,15,22,25,26,27,31,1,3,7,10,12,14,16,17,20,21,30,4,5,8,9,11,18,19,23,24,28,29)


par(mfrow=c(3,2))

filenames <- c('allvary','fixedeH','fixedcnA','fixedwC','fixedwH','fixedwN')
titles <- c('All Parameters Estimated', 'Fixed eH', 'Fixed C:N Phytoplankton',
            expression(paste('Fixed ', omega, ' C')), expression(paste('Fixed ', omega, ' H')),
            expression(paste('Fixed ', omega, ' N')))

for(i in 1:length(filenames)){
    load(file=paste('allochth_', filenames[i], '_raw.rda', sep=''))
    out <- reassembler_gen3(j1$samples)
    MLEs <- out[[1]]
    dens_chains <- out[[2]]
    comb_by_lake <- out[[3]] #rows are draws; cols are lg, sm, hi; depths are phyto, peri, terr

    plot(density(comb_by_lake[,1,3], n=5000, na.rm=TRUE), col=cols[1], lwd=2, main=titles[i])
    lines(density(comb_by_lake[,2,3], n=5000, na.rm=TRUE), col=cols[2], lwd=2)
    lines(density(comb_by_lake[,3,3], n=5000, na.rm=TRUE), col=cols[3], lwd=2)

    box(dens_chains)
}

# real sensitivity analysis (+/- 2 SD) (obsolete) ####


mat <- 1:4
counter <- 1
pars <- c('wC', 'wN', 'wH', 'eH', 'CNA')
parconds <- c('low', 'ave', 'hi')
ps <- c('phyto', 'peri', 'terr')
lakes <- c('sm','lg','hi')

for(par in 1:5){
    for(parcond in 1:3){
        load(file=paste0('sensitivity_2/allochth_bylake_', parconds[parcond], pars[par], '.rda'))

        for(i in 1:3){ #phyto, peri, terr
            for(j in 1:3){ #sm, lg, hi
                row <- 1:4
                counter <- counter + 1
                row[1] <- bylake_out$`credible intervals`[j,1,i]
                row[2] <- bylake_out$means[j,i]
                row[3] <- bylake_out$modes[j,i]
                row[4] <- bylake_out$`credible intervals`[j,2,i]

                mat <- rbind(mat, row)
                rownames(mat)[counter] <- paste(lakes[j], ps[i], parconds[parcond], pars[par], sep='_')
            }
        }
    }
}

mat <- mat[-1,]


lakenames <- c('Montane-Sm.', 'Montane-Lg', 'Alpine')
# consnames <- c('Calanoida', 'Trichoptera', 'Cladocera')
parnames <- list(expression(bold(paste(Delta^13,'C',sep=''))),
                 expression(bold(paste(Delta^15,'N',sep=''))),
                 expression(bold(paste(omega,sep=''))),
                 expression(bold(paste('C:N'['a']))),
                 expression(bold(paste(epsilon,'H',sep=''))))

par(mfrow=c(3,5), mar=c(.5, .5, .5, .5), oma=c(5, 7, 2, 6))

counter <- 0
for(lake in 1:3){
    for(par in 1:5){

        counter <- counter + 1

        low <- paste(lakes[lake], 'terr', 'low', pars[par], sep='_')
        ave <- paste(lakes[lake], 'terr', 'ave', pars[par], sep='_')
        hi <- paste(lakes[lake], 'terr', 'hi', pars[par], sep='_')

        low_low <- mat[,1][rownames(mat) == low]
        mean_low <- mat[,2][rownames(mat) == low]
        mode_low <- mat[,3][rownames(mat) == low]
        hi_low <- mat[,4][rownames(mat) == low]
        low_ave <- mat[,1][rownames(mat) == ave]
        mean_ave <- mat[,2][rownames(mat) == ave]
        mode_ave <- mat[,3][rownames(mat) == ave]
        hi_ave <- mat[,4][rownames(mat) == ave]
        low_hi <- mat[,1][rownames(mat) == hi]
        mean_hi <- mat[,2][rownames(mat) == hi]
        mode_hi <- mat[,3][rownames(mat) == hi]
        hi_hi <- mat[,4][rownames(mat) == hi]

        plot(1, 1, ylim=c(0,1), type='n', xlim=c(.8,3.2), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='')
        polygon(x=c(1,2,3,3,2,1), y=c(low_low, low_ave, low_hi, hi_hi, hi_ave, hi_low), col='gray90', border=FALSE)
        lines(x=1:3, y=c(mean_low, mean_ave, mean_hi), lty=2, lwd=2, col='red')
        lines(x=1:3, y=c(mode_low, mode_ave, mode_hi), lty=3, lwd=2, col='blue')

        if(counter %in% 1:5){
            mtext(parnames[[par]], 3, font=2)
        }
        if(counter %in% c(1,6,11)){
            mtext(lakenames[lake], 2, font=2, line=3)
            axis(2, at=c(0, .5, 1), labels=c(0.0, 0.5, 1.0), las=2)
        }
        if(counter %in% 11:15){
            axis(1, at=1:3, labels=c('-2 SD', 'mean', '+2 SD'))
            axis(1, col='white', tcl=0, at=1:3, labels=rep('',3))
        }
    }
}
mtext('Parameter Value', side=1, line=3, font=2, outer=TRUE, cex=1.15)
mtext('Allochthony', side=2, line=4.5, font=2, outer=TRUE, cex=1.15)
legend(x=3, y=1.58, xpd=NA, bty='n', lty=2:3, col=c('red', 'blue'), legend=c('Mean', 'Mode'), lwd=2, cex=1.3)


#wH biplots####
merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
    # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
    add.id.column.to.data <- function(DATA)
    {
        data.frame(DATA, id... = seq_len(nrow(DATA)))
    }
    # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
    order.by.id...and.remove.it <- function(DATA)
    {
        # gets in a data.frame with the "id..." column.  Orders by it and returns it
        if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")

        ss_r <- order(DATA$id...)
        ss_c <- colnames(DATA) != "id..."
        DATA[ss_r, ss_c]
    }

    # tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
    # tmp()

    if(!missing(keep_order))
    {
        if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
        if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
        # if you didn't get "return" by now - issue a warning.
        warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
    } else {return(merge(x=x,y=y,..., sort = sort))}
}
error.bars <- function(x, y, upper, lower=upper, cap.length=0.1, horiz=F,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("One or more vectors is not the same length")

    if(horiz==F) {
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
    } else if (horiz==T) {
        arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
    }
}
coda.samples.dic <- function (model, variable.names = NULL, n.iter, thin = 1, ...)
{
    load.module('dic') # necessary for pD and deviance monitor

    start <- model$iter() + thin
    varnames=c(variable.names, c('deviance', 'pD'))
    out <- jags.samples(model, varnames, n.iter, thin,
                        type = "trace", ...)
    deviance <- out$deviance
    pD <- out$pD
    out$deviance <- NULL
    out$pD <- NULL
    # ans <- vector("list", nchain(model))
    ans <- vector("list", model$nchain())
    for (ch in 1:model$nchain()) {
        ans.ch <- vector("list", length(out))
        vnames.ch <- NULL
        for (i in seq(along = out)) {
            varname <- names(out)[[i]]
            d <- dim(out[[i]])
            if (length(d) < 3) {
                stop("Invalid dimensions for sampled output")
            }
            vardim <- d[1:(length(d) - 2)]
            nvar <- prod(vardim)
            niter <- d[length(d) - 1]
            nchain <- d[length(d)]
            values <- as.vector(out[[i]])
            var.i <- matrix(NA, nrow = niter, ncol = nvar)
            for (j in 1:nvar) {
                var.i[, j] <- values[j + (0:(niter - 1)) * nvar +
                                         (ch - 1) * niter * nvar]
            }
            vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
            ans.ch[[i]] <- var.i
        }
        ans.ch <- do.call("cbind", ans.ch)
        colnames(ans.ch) <- vnames.ch
        ans[[ch]] <- mcmc(ans.ch, start = start, thin = thin)
    }

    dic <- list(deviance = mean(as.vector(deviance)), penalty = mean(as.vector(pD)), type = 'pD')
    class(dic) <- "dic"
    return(list(samples=mcmc.list(ans), dic=dic))
}
gmean <- function (x) {
    if (is.null(dim(x))) {
        exp(mean(log(x)))
    } else {
        t(apply(x, 1, function(y) {
            exp(mean(log(x)))
        }))
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
print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}
library(gtools)
library(rjags)
library(stringr)

# tiff(type='cairo', res=300, width=6, height=6, units='in',
#      filename='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/official_figures/figure_A1.tif')
par(mfcol=c(3,3), mar=c(3,2,2,2), oma=c(3,4.5,1,2))
SI <- read.csv("cor_loww.csv", header=TRUE)
# SI <- SI[c(13,23,21,27,17,24,14,19,18,16,28,30,22,29,25,26,31,1,2,3,4,15,5,6,7,8,9,10,11,12,20),] #converts to old order
SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)
ind_list2 <- list('phyto' = 1, 'peri' = 2, 'terr' = 3)
biplot_bylake <- function(tracer1, tracer2, label_seed=1, titles=TRUE, legend=TRUE){

    lakevec <- c(2,1,2,3,3,1,2,3,3,2,4,2,1,2,1,2,2,4,4,2,2,1,3,3,1,1,1,3,3,2,4)
    laketitles <- c('Montane-Lg.', 'Montane-Sm.', 'Alpine')
    agg <- aggregate(SI[,3:ncol(SI)], by=list(lakevec), FUN=mean)
    agg2 <- aggregate(SI[,3:ncol(SI)], by=list(lakevec,SI[,2]), FUN=mean)
    gg <- agg[1:3,]
    gg <- cbind(1:3, gg)

    inds <- setNames(c(0,4,8), c('N', 'C', 'H'))
    ind1 <- unname(inds[names(inds) == tracer1])
    ind2 <- unname(inds[names(inds) == tracer2])

    augs <- setNames(c(1,1,1), c('N', 'C', 'H'))
    aug1 <- unname(augs[names(augs) == tracer1])
    aug2 <- unname(augs[names(augs) == tracer2])

    # labs <- setNames(c("delta^15, 'N'", "delta^13, 'C'", "delta^2, 'H'"), c('N', 'C', 'H'))
    # lab1 <- unname(labs[names(labs) == tracer1])
    # lab2 <- unname(labs[names(labs) == tracer2])
    count <- label_seed - 1

    for(lake in c(2,1,3)){

        count <- count + 1

        row <- which(gg[,1] == lake)
        other_rows <- which(agg2[,1] == lake)
        means1 <- unlist(c(gg[row,(3:6)+ind1], agg2[other_rows,3+ind1]))
        means2 <- unlist(c(gg[row,(3:6)+ind2], agg2[other_rows,3+ind2]))
        sds1<- unlist(c(gg[row,(15:18)+ind1], agg2[other_rows,15+ind1]))
        sds2<- unlist(c(gg[row,(15:18)+ind2], agg2[other_rows,15+ind2]))


        xmin <- min(means1) - max(sds1)
        xmax <- max(means1) + max(sds1)
        ymin <- min(means2) - max(sds2)
        ymax <- max(means2) + max(sds2)


        x <- rbind(c(gg[row,4+ind1], gg[row,4+ind2] - gg[row,16+ind2]*aug2),
                   c(gg[row,4+ind1], gg[row,4+ind2] + gg[row,16+ind2]*aug2),
                   c(gg[row,4+ind1] - gg[row,16+ind1]*aug1, gg[row,4+ind2]),
                   c(gg[row,4+ind1] + gg[row,16+ind1]*aug1, gg[row,4+ind2]),
                   c(gg[row,5+ind1], gg[row,5+ind2] - gg[row,17+ind2]),
                   c(gg[row,5+ind1], gg[row,5+ind2] + gg[row,17+ind2]),
                   c(gg[row,5+ind1] - gg[row,17+ind1], gg[row,5+ind2]),
                   c(gg[row,5+ind1] + gg[row,17+ind1], gg[row,5+ind2]),
                   c(gg[row,6+ind1], gg[row,6+ind2] - gg[row,18+ind2]),
                   c(gg[row,6+ind1], gg[row,6+ind2] + gg[row,18+ind2]),
                   c(gg[row,6+ind1] - gg[row,18+ind1], gg[row,6+ind2]),
                   c(gg[row,6+ind1] + gg[row,18+ind1], gg[row,6+ind2]))

        hullind <- chull(x)

        err <- 'black'; pt <- 'black'

        if(titles){
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1,
                 main='')
        } else {
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1)
        }

        polygon(x[hullind,][,1], x[hullind,][,2], col='gray90', border=NA)


        # print.letter(letters[count], xy=c(0.9,0.9), cex=1.3, font=2)
        #peri
        error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind1], horiz=TRUE,
                   cap.length=0.0, col='darkgreen', lwd=3)
        error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind2], horiz=FALSE,
                   cap.length=0.0, col='darkgreen', lwd=3)
        #phyto
        error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind1]*aug1, horiz=TRUE,
                   cap.length=0.0, col='steelblue1', lwd=3)
        error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind2]*aug2, horiz=FALSE,
                   cap.length=0.0, col='steelblue1', lwd=3)
        #terr
        error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind1], horiz=TRUE,
                   cap.length=0.0, col='sienna4', lwd=3)
        error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind2], horiz=FALSE,
                   cap.length=0.0, col='sienna4', lwd=3)

        points(gg[row,5+ind1], gg[row,5+ind2], #peri
               pch=20, cex=2.5, col='darkgreen')
        points(gg[row,4+ind1], gg[row,4+ind2], #phyto
               pch=20, cex=2.5, col='steelblue1')
        points(gg[row,6+ind1], gg[row,6+ind2], #terr
               pch=20, cex=2.5, col='sienna4')

        #cons
        for(cons in c('Calanoida', 'Cladocera', 'Trichoptera')){
            if(lake == 3 & cons == 'Cladocera') next
            row2 <- which(agg2[,1] == lake & agg2[,2] == cons)
            error.bars(agg2[row2,3+ind1], agg2[row2,3+ind2], upper=agg2[row2,15+ind1], horiz=TRUE,
                       cap.length=0.015, col=err, lwd=1.4)
            error.bars(agg2[row2,3+ind1], agg2[row2,3+ind2], upper=agg2[row2,15+ind2], horiz=FALSE,
                       cap.length=0.015, col=err, lwd=1.4)
        }
        for(cons in c('Calanoida', 'Cladocera', 'Trichoptera')){
            if(lake == 3 & cons == 'Cladocera') next
            row2 <- which(agg2[,1] == lake & agg2[,2] == cons)
            if(cons == 'Calanoida'){
                points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=21, cex=1.5, col=pt, bg='white')
            } else {
                if(cons == 'Cladocera'){
                    points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=25, cex=1.5, col=pt, bg='white')
                } else {
                    points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=24, cex=1.5, col=pt, bg='white')
                }
            }
        }

        if(count == 1){
            mtext(expression(paste(omega, ': -2 SD')), 3, line=1.2, cex=1.2)
            # mtext('Mean', 3, font=2, line=1.2)
            # mtext('+2 SD', 3, font=2, line=1.2)
            mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5, cex=1.2)
        }
        if(count == 2){
            mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=4.5, side=2)
            mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5, cex=1.2)
        }
        if(count == 3){
            # mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=3, side=1)
            mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5, cex=1.2)
        }
    }

    if(legend==T){
        legend(x=-135, y=-48, legend=c('Phytoplank.', 'Periphyton', 'Terrestrial', 'Calanoida', 'Cladocera', 'Trichoptera'),
               pch=c(20,20,20,21,25,24), pt.bg='white', cex=1.5, bty='n', xpd=NA,
               col=c('steelblue1','darkgreen','sienna4','black','black','black'))
    }
}
biplot_bylake('H', 'C', label_seed=1, titles=T, legend=F)

tracer1='H';tracer2='C';label_seed=1;titles=T;legend=T

SI <- read.csv("cor_avew.csv", header=TRUE)
# SI <- SI[c(13,23,21,27,17,24,14,19,18,16,28,30,22,29,25,26,31,1,2,3,4,15,5,6,7,8,9,10,11,12,20),] #converts to old order
SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)
ind_list2 <- list('phyto' = 1, 'peri' = 2, 'terr' = 3)
biplot_bylake <- function(tracer1, tracer2, label_seed=1, titles=TRUE, legend=TRUE){

    lakevec <- c(2,1,2,3,3,1,2,3,3,2,4,2,1,2,1,2,2,4,4,2,2,1,3,3,1,1,1,3,3,2,4)
    laketitles <- c('Montane-Lg.', 'Montane-Sm.', 'Alpine')
    agg <- aggregate(SI[,3:ncol(SI)], by=list(lakevec), FUN=mean)
    agg2 <- aggregate(SI[,3:ncol(SI)], by=list(lakevec,SI[,2]), FUN=mean)
    gg <- agg[1:3,]
    gg <- cbind(1:3, gg)

    inds <- setNames(c(0,4,8), c('N', 'C', 'H'))
    ind1 <- unname(inds[names(inds) == tracer1])
    ind2 <- unname(inds[names(inds) == tracer2])

    augs <- setNames(c(1,1,1), c('N', 'C', 'H'))
    aug1 <- unname(augs[names(augs) == tracer1])
    aug2 <- unname(augs[names(augs) == tracer2])

    # labs <- setNames(c("delta^15, 'N'", "delta^13, 'C'", "delta^2, 'H'"), c('N', 'C', 'H'))
    # lab1 <- unname(labs[names(labs) == tracer1])
    # lab2 <- unname(labs[names(labs) == tracer2])
    count <- label_seed - 1

    for(lake in c(2,1,3)){

        count <- count + 1

        row <- which(gg[,1] == lake)
        other_rows <- which(agg2[,1] == lake)
        means1 <- unlist(c(gg[row,(3:6)+ind1], agg2[other_rows,3+ind1]))
        means2 <- unlist(c(gg[row,(3:6)+ind2], agg2[other_rows,3+ind2]))
        sds1<- unlist(c(gg[row,(15:18)+ind1], agg2[other_rows,15+ind1]))
        sds2<- unlist(c(gg[row,(15:18)+ind2], agg2[other_rows,15+ind2]))


        xmin <- min(means1) - max(sds1)
        xmax <- max(means1) + max(sds1)
        ymin <- min(means2) - max(sds2)
        ymax <- max(means2) + max(sds2)


        x <- rbind(c(gg[row,4+ind1], gg[row,4+ind2] - gg[row,16+ind2]*aug2),
                   c(gg[row,4+ind1], gg[row,4+ind2] + gg[row,16+ind2]*aug2),
                   c(gg[row,4+ind1] - gg[row,16+ind1]*aug1, gg[row,4+ind2]),
                   c(gg[row,4+ind1] + gg[row,16+ind1]*aug1, gg[row,4+ind2]),
                   c(gg[row,5+ind1], gg[row,5+ind2] - gg[row,17+ind2]),
                   c(gg[row,5+ind1], gg[row,5+ind2] + gg[row,17+ind2]),
                   c(gg[row,5+ind1] - gg[row,17+ind1], gg[row,5+ind2]),
                   c(gg[row,5+ind1] + gg[row,17+ind1], gg[row,5+ind2]),
                   c(gg[row,6+ind1], gg[row,6+ind2] - gg[row,18+ind2]),
                   c(gg[row,6+ind1], gg[row,6+ind2] + gg[row,18+ind2]),
                   c(gg[row,6+ind1] - gg[row,18+ind1], gg[row,6+ind2]),
                   c(gg[row,6+ind1] + gg[row,18+ind1], gg[row,6+ind2]))

        hullind <- chull(x)

        err <- 'black'; pt <- 'black'

        if(titles){
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1,
                 main='')
        } else {
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1)
        }

        polygon(x[hullind,][,1], x[hullind,][,2], col='gray90', border=NA)


        # print.letter(letters[count], xy=c(0.9,0.9), cex=1.3, font=2)
        #peri
        error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind1], horiz=TRUE,
                   cap.length=0.0, col='darkgreen', lwd=3)
        error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind2], horiz=FALSE,
                   cap.length=0.0, col='darkgreen', lwd=3)
        #phyto
        error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind1]*aug1, horiz=TRUE,
                   cap.length=0.0, col='steelblue1', lwd=3)
        error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind2]*aug2, horiz=FALSE,
                   cap.length=0.0, col='steelblue1', lwd=3)
        #terr
        error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind1], horiz=TRUE,
                   cap.length=0.0, col='sienna4', lwd=3)
        error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind2], horiz=FALSE,
                   cap.length=0.0, col='sienna4', lwd=3)

        points(gg[row,5+ind1], gg[row,5+ind2], #peri
               pch=20, cex=2.5, col='darkgreen')
        points(gg[row,4+ind1], gg[row,4+ind2], #phyto
               pch=20, cex=2.5, col='steelblue1')
        points(gg[row,6+ind1], gg[row,6+ind2], #terr
               pch=20, cex=2.5, col='sienna4')

        #cons
        for(cons in c('Calanoida', 'Cladocera', 'Trichoptera')){
            if(lake == 3 & cons == 'Cladocera') next
            row2 <- which(agg2[,1] == lake & agg2[,2] == cons)
            error.bars(agg2[row2,3+ind1], agg2[row2,3+ind2], upper=agg2[row2,15+ind1], horiz=TRUE,
                       cap.length=0.015, col=err, lwd=1.4)
            error.bars(agg2[row2,3+ind1], agg2[row2,3+ind2], upper=agg2[row2,15+ind2], horiz=FALSE,
                       cap.length=0.015, col=err, lwd=1.4)
        }
        for(cons in c('Calanoida', 'Cladocera', 'Trichoptera')){
            if(lake == 3 & cons == 'Cladocera') next
            row2 <- which(agg2[,1] == lake & agg2[,2] == cons)
            if(cons == 'Calanoida'){
                points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=21, cex=1.5, col=pt, bg='white')
            } else {
                if(cons == 'Cladocera'){
                    points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=25, cex=1.5, col=pt, bg='white')
                } else {
                    points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=24, cex=1.5, col=pt, bg='white')
                }
            }
        }

        if(count == 1){
            # mtext('-2 SD', 3, font=2, line=1.2)
            mtext(expression(paste(omega, ': Mean')), 3, font=2, line=1.2, cex=1.2)
            # mtext('+2 SD', 3, font=2, line=1.2)
            # mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5)
        }
        if(count == 2){
            # mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=4, side=2)
            # mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5)
        }
        if(count == 3){
            mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=3, side=1)
            # mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5)
        }
    }

    if(legend==T){
        legend(x=-135, y=-48, legend=c('Phytoplank.', 'Periphyton', 'Terrestrial', 'Calanoida', 'Cladocera', 'Trichoptera'),
               pch=c(20,20,20,21,25,24), pt.bg='white', cex=1.5, bty='n', xpd=NA,
               col=c('steelblue1','darkgreen','sienna4','black','black','black'))
    }
}
biplot_bylake('H', 'C', label_seed=1, titles=T, legend=F)

SI <- read.csv("cor_hiw.csv", header=TRUE)
# SI <- SI[c(13,23,21,27,17,24,14,19,18,16,28,30,22,29,25,26,31,1,2,3,4,15,5,6,7,8,9,10,11,12,20),] #converts to old order
SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)
ind_list2 <- list('phyto' = 1, 'peri' = 2, 'terr' = 3)
biplot_bylake <- function(tracer1, tracer2, label_seed=1, titles=TRUE, legend=TRUE){

    lakevec <- c(2,1,2,3,3,1,2,3,3,2,4,2,1,2,1,2,2,4,4,2,2,1,3,3,1,1,1,3,3,2,4)
    laketitles <- c('Montane-Lg.', 'Montane-Sm.', 'Alpine')
    agg <- aggregate(SI[,3:ncol(SI)], by=list(lakevec), FUN=mean)
    agg2 <- aggregate(SI[,3:ncol(SI)], by=list(lakevec,SI[,2]), FUN=mean)
    gg <- agg[1:3,]
    gg <- cbind(1:3, gg)

    inds <- setNames(c(0,4,8), c('N', 'C', 'H'))
    ind1 <- unname(inds[names(inds) == tracer1])
    ind2 <- unname(inds[names(inds) == tracer2])

    augs <- setNames(c(1,1,1), c('N', 'C', 'H'))
    aug1 <- unname(augs[names(augs) == tracer1])
    aug2 <- unname(augs[names(augs) == tracer2])

    # labs <- setNames(c("delta^15, 'N'", "delta^13, 'C'", "delta^2, 'H'"), c('N', 'C', 'H'))
    # lab1 <- unname(labs[names(labs) == tracer1])
    # lab2 <- unname(labs[names(labs) == tracer2])
    count <- label_seed - 1

    for(lake in c(2,1,3)){

        count <- count + 1

        row <- which(gg[,1] == lake)
        other_rows <- which(agg2[,1] == lake)
        means1 <- unlist(c(gg[row,(3:6)+ind1], agg2[other_rows,3+ind1]))
        means2 <- unlist(c(gg[row,(3:6)+ind2], agg2[other_rows,3+ind2]))
        sds1<- unlist(c(gg[row,(15:18)+ind1], agg2[other_rows,15+ind1]))
        sds2<- unlist(c(gg[row,(15:18)+ind2], agg2[other_rows,15+ind2]))


        xmin <- min(means1) - max(sds1)
        xmax <- max(means1) + max(sds1)
        ymin <- min(means2) - max(sds2)
        ymax <- max(means2) + max(sds2)


        x <- rbind(c(gg[row,4+ind1], gg[row,4+ind2] - gg[row,16+ind2]*aug2),
                   c(gg[row,4+ind1], gg[row,4+ind2] + gg[row,16+ind2]*aug2),
                   c(gg[row,4+ind1] - gg[row,16+ind1]*aug1, gg[row,4+ind2]),
                   c(gg[row,4+ind1] + gg[row,16+ind1]*aug1, gg[row,4+ind2]),
                   c(gg[row,5+ind1], gg[row,5+ind2] - gg[row,17+ind2]),
                   c(gg[row,5+ind1], gg[row,5+ind2] + gg[row,17+ind2]),
                   c(gg[row,5+ind1] - gg[row,17+ind1], gg[row,5+ind2]),
                   c(gg[row,5+ind1] + gg[row,17+ind1], gg[row,5+ind2]),
                   c(gg[row,6+ind1], gg[row,6+ind2] - gg[row,18+ind2]),
                   c(gg[row,6+ind1], gg[row,6+ind2] + gg[row,18+ind2]),
                   c(gg[row,6+ind1] - gg[row,18+ind1], gg[row,6+ind2]),
                   c(gg[row,6+ind1] + gg[row,18+ind1], gg[row,6+ind2]))

        hullind <- chull(x)

        err <- 'black'; pt <- 'black'

        if(titles){
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1,
                 main='')
        } else {
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1)
        }

        polygon(x[hullind,][,1], x[hullind,][,2], col='gray90', border=NA)


        # print.letter(letters[count], xy=c(0.9,0.9), cex=1.3, font=2)
        #peri
        error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind1], horiz=TRUE,
                   cap.length=0.0, col='darkgreen', lwd=3)
        error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind2], horiz=FALSE,
                   cap.length=0.0, col='darkgreen', lwd=3)
        #phyto
        error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind1]*aug1, horiz=TRUE,
                   cap.length=0.0, col='steelblue1', lwd=3)
        error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind2]*aug2, horiz=FALSE,
                   cap.length=0.0, col='steelblue1', lwd=3)
        #terr
        error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind1], horiz=TRUE,
                   cap.length=0.0, col='sienna4', lwd=3)
        error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind2], horiz=FALSE,
                   cap.length=0.0, col='sienna4', lwd=3)

        points(gg[row,5+ind1], gg[row,5+ind2], #peri
               pch=20, cex=2.5, col='darkgreen')
        points(gg[row,4+ind1], gg[row,4+ind2], #phyto
               pch=20, cex=2.5, col='steelblue1')
        points(gg[row,6+ind1], gg[row,6+ind2], #terr
               pch=20, cex=2.5, col='sienna4')

        #cons
        for(cons in c('Calanoida', 'Cladocera', 'Trichoptera')){
            if(lake == 3 & cons == 'Cladocera') next
            row2 <- which(agg2[,1] == lake & agg2[,2] == cons)
            error.bars(agg2[row2,3+ind1], agg2[row2,3+ind2], upper=agg2[row2,15+ind1], horiz=TRUE,
                       cap.length=0.015, col=err, lwd=1.4)
            error.bars(agg2[row2,3+ind1], agg2[row2,3+ind2], upper=agg2[row2,15+ind2], horiz=FALSE,
                       cap.length=0.015, col=err, lwd=1.4)
        }
        for(cons in c('Calanoida', 'Cladocera', 'Trichoptera')){
            if(lake == 3 & cons == 'Cladocera') next
            row2 <- which(agg2[,1] == lake & agg2[,2] == cons)
            if(cons == 'Calanoida'){
                points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=21, cex=1.5, col=pt, bg='white')
            } else {
                if(cons == 'Cladocera'){
                    points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=25, cex=1.5, col=pt, bg='white')
                } else {
                    points(agg2[row2,3+ind1], agg2[row2,3+ind2], pch=24, cex=1.5, col=pt, bg='white')
                }
            }
        }

        if(count == 1){
            # mtext('-2 SD', 3, font=2, line=1.2)
            # mtext('Mean', 3, font=2, line=1.2)
            mtext(expression(paste(omega, ': +2 SD')), 3, font=2, line=1.2, cex=1.2)
            # mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5)
        }
        if(count == 2){
            # mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=4, side=2)
            # mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5)
        }
        if(count == 3){
            # mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=3, side=1)
            # mtext(paste(laketitles[which(1:3==lake)]), side=2, line=2.5)
        }
    }

    if(legend==T){
        legend(x=-135, y=-48, legend=c('Phytoplank.', 'Periphyton', 'Terrestrial', 'Calanoida', 'Cladocera', 'Trichoptera'),
               pch=c(20,20,20,21,25,24), pt.bg='white', cex=1.5, bty='n', xpd=NA,
               col=c('steelblue1','darkgreen','sienna4','black','black','black'))
    }
}
biplot_bylake('H', 'C', label_seed=1, titles=T, legend=F)
# dev.off()

#correction param posterior plot####
load('allvary_C.rda')
load('allvary_N.rda')
load('allvary_H.rda')

colnames(allvary_C[[1]])
colnames(allvary_N[[1]])
colnames(allvary_H[[1]])

#exploratory
par(mfrow=c(3,5))
for(i in 1:13){ #decide which variables to plot here (line them up with the colnames just above)
    plot(density(allvary_C[[2]][1,i,], na.rm=T), main=colnames(allvary_C[[1]])[i])
    for(j in 2:sum(!(is.na(allvary_C[[1]][,i])))){
        lines(density(allvary_C[[2]][j,i,], na.rm=T))
    }
}
par(mfrow=c(3,5))
for(i in 1:13){ #decide which variables to plot here (line them up with the colnames just above)
    plot(density(allvary_N[[2]][1,i,], na.rm=T), main=colnames(allvary_N[[1]])[i])
    for(j in 2:sum(!(is.na(allvary_N[[1]][,i])))){
        lines(density(allvary_N[[2]][j,i,], na.rm=T))
    }
}
par(mfrow=c(3,5))
for(i in 1:9){ #decide which variables to plot here (line them up with the colnames just above)
    plot(density(allvary_H[[2]][1,i,], na.rm=T), main=colnames(allvary_H[[1]])[i])
    for(j in 2:sum(!(is.na(allvary_H[[1]][,i])))){
        lines(density(allvary_H[[2]][j,i,], na.rm=T))
    }
}


#for reals
groupvec1 <- c(1,1,2,2,2,2,3,3,3,1,3)
groupvec2 <- c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1)
specvec <- c(rep(1,10),2,rep(3,6),4,4,rep(5,11),2)
l <- 'springgreen3'; h <- 'darkgoldenrod2'; s <- 'royalblue2'
col_ind <- list('1'=l, '2'=s, '3'=h)
col_ind2 <- list('1'='green', '2'='orange', '3'='blue', '4' = 'black', '5'='purple')

par(mfrow=c(3,6), mar=c(4, 2, 3, 1), oma=c(2, 4, 0, 0))

#mean_d13Ca_lake C
column = 10 #set
vec = groupvec1 #set
plot(density(allvary_C[[2]][1,column,], na.rm=T), main='phyto C', xaxs='i', yaxs='i', #settitle
     xlim=c(-40,-18), ylim=c(0,.55), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
    lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#mean_d13Cperi_lake C
column = 11 #set
vec = groupvec1 #set
plot(density(allvary_C[[2]][1,column,], na.rm=T), main='peri C', xaxs='i', yaxs='i', #settitle
     xlim=c(-40,-5), ylim=c(0,.38), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
    lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#d13CT C
column = 5 #set
vec = groupvec1 #set
plot(density(allvary_C[[2]][1,column,], na.rm=T), main='terr C', xaxs='i', yaxs='i', #settitle
     ylim=c(0,.16), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
    lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#mean_cons_true C - byclass
column = 9 #set
vec = groupvec2 #set
plot(density(allvary_C[[2]][1,column,], na.rm=T), main='cons C', xaxs='i', yaxs='i', #settitle
     ylim=c(0,.10), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
    lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#d13C_enrich C
column = 2 #set
vec = groupvec2 #set
plot(density(allvary_C[[2]][1,column,], na.rm=T), main='troph enrich C', xaxs='i', yaxs='i',
     xlim=c(-4,5), ylim=c(0,.39), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
    lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

plot(1,1, type='n', ann=F, axes=FALSE, main='')
legend('center', legend=c('Montane-Small', 'Montane-large', 'Alpine'), fill=c(s,l,h))

#---#

#mean_d15Na_lake N
column = 10 #set
vec = groupvec1 #set
plot(density(allvary_N[[2]][1,column,], na.rm=T), main='phyto N', xaxs='i', yaxs='i', #settitle
     xlim=c(-3.5,5.3), ylim=c(0,.67), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_N[[1]][,column])))){
    lines(density(allvary_N[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#mean_d15Nperi_lake N
column = 11 #set
vec = groupvec1 #set
plot(density(allvary_N[[2]][1,column,], na.rm=T), main='peri N', xaxs='i', yaxs='i', #settitle
     xlim=c(-12,5), ylim=c(0,.55), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_N[[1]][,column])))){
    lines(density(allvary_N[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#d15NT N
column = 5 #set
vec = groupvec1 #set
plot(density(allvary_N[[2]][1,column,], na.rm=T), main='terr N', xaxs='i', yaxs='i', #settitle
     xlim=c(-23,15), ylim=c(0,.19), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_N[[1]][,column])))){
    lines(density(allvary_N[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#mean_cons_true N - byclass
column = 9 #set
vec = groupvec2 #set
plot(density(allvary_N[[2]][1,column,], na.rm=T), main='cons N', xaxs='i', yaxs='i', #settitle
     col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_N[[1]][,column])))){
    lines(density(allvary_N[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#d15N_enrich N
column = 2 #set
vec = groupvec2 #set
plot(density(allvary_N[[2]][1,column,], na.rm=T), main='troph enrich N', xaxs='i', yaxs='i',
     ylim=c(0,.45), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_N[[1]][,column])))){
    lines(density(allvary_N[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#mean_cnA_lake C
column = 8 #set
vec = groupvec1 #set
plot(density(allvary_C[[2]][1,column,], na.rm=T), main='C:N phyto', xaxs='i', yaxs='i',  #settitle
     xlim=c(-30,35), ylim=c(0,.14), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
    lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#---#

#mean_d2Ha_lake H
column = 4 #set
vec = groupvec1 #set
plot(density(allvary_H[[2]][1,column,], na.rm=T), main='phyto H', xaxs='i', yaxs='i', #settitle
     ylim=c(0,.02), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_H[[1]][,column])))){
    lines(density(allvary_H[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#mean_d2Hperi_lake H
column = 3 #set
vec = groupvec1 #set
plot(density(allvary_H[[2]][1,column,], na.rm=T), main='peri H', xaxs='i', yaxs='i', #settitle
     xlim=c(-300,-100), ylim=c(0,.06), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_H[[1]][,column])))){
    lines(density(allvary_H[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#d2HT H
load("withd2HT.rda") #alt H
column =  #set
vec = groupvec1 #set
plot(density(allvary_H[[2]][1,column,], na.rm=T), main='terr H', xaxs='i', yaxs='i', #settitle
     col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_H[[1]][,column])))){
    lines(density(allvary_H[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}
load('allvary_H.rda')

#mean_cons_true H - byclass
column = 6 #set
vec = groupvec2 #set
plot(density(allvary_H[[2]][1,column,], na.rm=T), main='cons H', xaxs='i', yaxs='i', #settitle
     col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_H[[1]][,column])))){
    lines(density(allvary_H[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#env water
column = 5 #set
vec = groupvec2 #set
plot(density(allvary_H[[2]][1,column,], na.rm=T), main='env water', xaxs='i', yaxs='i',
      ylim=c(0,2.8), col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_H[[1]][,column])))){
    lines(density(allvary_H[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

#photosynth discrim H
column = 2 #set
vec = groupvec1 #set
plot(density(allvary_H[[2]][1,column,], na.rm=T), main='photosynth discrim', xaxs='i', yaxs='i',
     col=col_ind[[paste(vec[1])]], ylab='') #setlim
for(j in 2:sum(!(is.na(allvary_H[[1]][,column])))){
    lines(density(allvary_H[[2]][j,column,], na.rm=T), col=col_ind[[paste(vec[j])]])
}

mtext('Posterior Density', side=2, outer=TRUE, line=1.5)



#mean_cons_true C - byspec
# column = 9 #set
# vec = specvec #set
# plot(density(allvary_C[[2]][1,column,], na.rm=T), main='cons C', xaxs='i', yaxs='i', #settitle
#      ylim=c(0,.10), col=col_ind2[[paste(vec[1])]]) #setlim
# for(j in 2:sum(!(is.na(allvary_C[[1]][,column])))){
#     lines(density(allvary_C[[2]][j,column,], na.rm=T), col=col_ind2[[paste(vec[j])]])
# }

# allvary_H[[1]]
