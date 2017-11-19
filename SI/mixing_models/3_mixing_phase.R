
rm(list=ls(all=TRUE))
# if (is.null(dev.list()) == FALSE){dev.off()}
# windows(record=T)
setwd("/home/mike/git/seven_lakes/SI/mixing_models")
#commented out inits, changed plot limits for posterior

# added functions ####
# library(RColorBrewer)
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
print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}

# SI - setup ####
# SI <- read.csv("cor_agg.csv", header=TRUE)
# SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)
# SI <- read.csv("source-as-cons_test.csv", header=TRUE)
SI <- read.csv("cor_all.csv", header=TRUE)
SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)

# library(gtools)
library(rjags)
# library(stringr)

# SI - set parameters (FA params set themselves for now) ####
isos <- c('C', 'N', 'H')
which_prey <- c('phyto', 'peri', 'terr')
npreytypes <- 3
npreds <- nrow(SI)
# nlake <- 11
# predstypes_SI <- c(rep(1,10),2,rep(3,6),4,4,rep(5,11),2)
# preytypes_SI <- rep(1:npreytypes, npreds)
# laketypes <- c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1)
# isovec <- preytypes_SI

#these ones will set themselves
niso <- length(isos)
# npredstypes <- length(levels(factor(predstypes_SI)))
# nlaketypes <- length(levels(factor(laketypes)))
alpha1 <- rep(1,npreytypes)

# SI - arrange data for models ####
prey_mean_data = array(NA, dim=c(3,niso,npreds))
prey_var_data = array(NA, dim=c(3,niso,npreds))

#arrays of producer values (prey X iso X preds)
index <- seq(3,11,4)
phytomean_temp <- t(SI[,1+index])
perimean_temp <- t(SI[,2+index])
terrmean_temp <- t(SI[,3+index])
phytosd_temp <- t(SI[,13+index])
perisd_temp <- t(SI[,14+index])
terrsd_temp <- t(SI[,15+index])

ind_list <- list('N' = 1, 'C' = 2, 'H' = 3)
isoindex <- as.vector(unlist(ind_list[isos]))
ind_list2 <- list('phyto' = 1, 'peri' = 2, 'terr' = 3)
preyindex <- as.vector(unlist(ind_list2[which_prey]))

if('phyto' %in% which_prey){
    prey_mean_data[1,1:niso,1:npreds] <- phytomean_temp[isoindex,]
    prey_var_data[1,1:niso,1:npreds] <- phytosd_temp[isoindex,]^2
}
if('peri' %in% which_prey){
    prey_mean_data[2,1:niso,1:npreds] <- perimean_temp[isoindex,]
    prey_var_data[2,1:niso,1:npreds] <- perisd_temp[isoindex,]^2
}
if('terr' %in% which_prey){
    prey_mean_data[3,1:niso,1:npreds] <- terrmean_temp[isoindex,]
    prey_var_data[3,1:niso,1:npreds] <- terrsd_temp[isoindex,]^2
}

if(npreytypes == 2){
    prey_mean_data <- prey_mean_data[-which(is.na(prey_mean_data[,1,1])),,]
    prey_var_data <- prey_var_data[-which(is.na(prey_var_data[,1,1])),,]
}

#matrix of predator values
preds_mean_data = as.matrix(SI[,index][,isoindex])

# #prey aggregates (means by lake group and prey types)
# #aggregate SI dataframe columns by lake
# if(npreytypes == 3){
#     agg_ind <- c(1:4, 13:16)
# } else {
#     if('phyto' %in% which_prey & 'peri' %in% which_prey){
#         agg_ind <- c(1:3, 13:15)
#     } else {
#         if('phyto' %in% which_prey & 'terr' %in% which_prey){
#             agg_ind <- c(1,2,4,13,14,16)
#         } else {
#             if('terr' %in% which_prey & 'peri' %in% which_prey){
#                 agg_ind <- c(1,3,4,13,15,16)
#             }
#         }
#     }
# }
# agg_list <- list('N' = agg_ind, 'C' = agg_ind+4, 'H' = agg_ind+8)
# agg_cols <- sort(as.vector(unlist(agg_list[isos])))
#
# temp <- matrix(nrow=nlake, ncol=length(agg_cols))
# for(col in 1:length(agg_cols)){
#     temp[,col] <- aggregate(SI[,agg_cols[col]+2], by=list(SI$lake_name), mean)[,2]
# }
#
# #then isolate prey cols and aggregate by lake type
# if(npreytypes == 3){
#     predscols <- seq(1,(length(agg_cols)-3),4)
# } else {
#     predscols <- seq(1,(ncol(temp)-2),3)
# }
#
# prey_pre_agg <- temp[,-predscols]
# temp2 <- matrix(nrow=nlaketypes, ncol=ncol(prey_pre_agg))
# for(col in 1:ncol(temp2)){
#     temp2[,col] <- aggregate(prey_pre_agg[,col], by=list(c(2,1,2,3,3,1,1,2,3,3,2)), mean)[,2]
# }
# # colnames(temp2) <- colnames(SI[,agg_cols+2])[-predscols]
# # rownames(temp2) <- c('Large', 'Small', 'High')
# #then reassemble into arrays for mean and var
# temp2mean <- temp2[,1:(length(temp2[1,])/2)]
# temp2sd <- temp2[,-(1:(length(temp2[1,])/2))]
#
# prey_mean_agg <- array(data=NA, dim=c(npreytypes,3,niso))
# prey_var_agg <- array(data=NA, dim=c(npreytypes,3,niso))
# if(npreytypes == 3){
#     for(i in 0:(niso-1)){
#         prey_mean_agg[,,i+1] <- t(temp2mean[,(1:3)+i*3])
#         prey_var_agg[,,i+1] <- t(temp2sd[,(1:3)+i*3])^2
#     }
# } else {
#     for(i in 0:(niso-1)){
#         prey_mean_agg[,,i+1] <- t(temp2mean[,(1:2)+i*2])
#         prey_var_agg[,,i+1] <- t(temp2sd[,(1:2)+i*2])^2
#     }
# }

#run model ####
jags_data = list('niso'=niso, 'npreytypes'=npreytypes, 'npreds'=npreds,
                 'prey_mean_data1'=prey_mean_data, 'prey_mean_data2'=prey_mean_data,
                 'prey_var'=prey_var_data, 'preds_mean_data'=preds_mean_data,
                 'alpha'=alpha1)
inits1 <- list('p'=matrix(c(0.1, 0.1, 0.8), nrow=npreytypes, ncol=npreds, byrow=TRUE))
               # 'e'=rep(0, niso)), 'preytau'=rep(0.01, niso))
inits2 <- list('p'=matrix(c(0.8, 0.1, 0.1), nrow=npreytypes, ncol=npreds, byrow=TRUE))
# 'e'=rep(100, niso)), 'preytau'=rep(1, niso))
inits3 <- list('p'=matrix(c(0.1, 0.8, 0.1), nrow=npreytypes, ncol=npreds, byrow=TRUE))
               # 'e'=rep(10, niso)), 'preytau'=rep(0.1, niso))
inits = list(inits1, inits2, inits3)
model = 'mod_mix.txt'
jags_params = c('p', 'preds_mean', 'preds_mean_data','prey_mean_data2', 'prey_mean')

# out <- jags.model(model, data=jags_data, inits=inits, n.chains=3, n.adapt=3000)
mod <- jags.model(model, data=jags_data, inits=inits, n.chains=3, n.adapt=3000)
update(mod, n.iter=5e3)
j1 <- coda.samples.dic(model=mod, variable.names=jags_params, n.iter=1e5, thin=100)
DIC <- j1$dic

#convergence
# gelman.plot(cons_cor) #shouldn't vary from 1 by more than .05
g <- matrix(NA, nrow=nvar(j1$samples), ncol=2)
for (v in 1:nvar(j1$samples)) {
    g[v,] <- gelman.diag(j1$samples[,v])$psrf
}; min(g, na.rm=T); max(g, na.rm=T)
# which(g[,1] > 1.1 | g[,2] > 1.1)i
# g[g[,1] > 1.1 | g[,2] > 1.1,]
# windows(record=T)
# plot(cons_cor)

#extract
reassembler_MCMClist <- function(mcmc_list,var){ #return 'MLEs' and densities as arrays
    nchains <- length(mcmc_list)
    chain_length <- length(mcmc_list[[1]][,1])
    param_list <- gsub(" *\\[.*?\\] *", '', colnames(mcmc_list[[1]]))
    all_ps <- which(param_list == var)
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
    MLE_arr[1,1:npreds] <- MLE_vec[indices]
    MLE_arr[2,1:npreds] <- MLE_vec[indices+1]
    if(npreytypes == 3){
        MLE_arr[3,1:npreds] <- MLE_vec[indices+2]
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

    return(list(MLE_arr, dens_arr))
}
out <- reassembler_MCMClist(j1$samples, 'p')
DIC <- j1$dic
# save(j1, file='allochth.rda')
MLEs <- out[[1]]
dens_chains <- out[[2]]

#testing####
# predrows <- 1:31
# SI3 <- SI[predrows,]
# preds_mean_data <- cbind(SI3$cons_N_u,SI3$cons_H_u)
# preds_mean_data[,1] <- scale(preds_mean_data[,1])
# preds_mean_data[,2] <- scale(preds_mean_data[,2])
# prey_mean_data <- prey_mean_data[,,predrows]
# prey_mean_data <- apply(prey_mean_data, c(2,3), scale)
# prey_var_data <- prey_var_data[,,predrows]
# prey_var_data <- array(0.1, dim=c(3,2,length(predrows)))
# jags_data = list('niso'=2, 'npreytypes'=3, 'npreds'=length(predrows),
#                  'prey_mean_data'=prey_mean_data, 'prey_var'=prey_var_data,
#                  'preds_mean_data'=preds_mean_data,
#                  'alpha'=alpha1)
# pvals = c(0.33333, 0.33333, 0.33334)
# inits = function() { list('p' =  matrix(pvals, nrow=3, ncol=3, byrow=T)) }
# model = 'mod_alloch.txt'
# jags_params = c('p', 'preds_mean', 'preds_mean_data','prey_mean_data')

#plot ####
biplot_bylake <- function(tracer1, tracer2, label_seed=1, titles=TRUE, legend=TRUE){
    lakevec <- c(2,1,2,3,3,1,2,3,3,2,4,2,1,2,1,2,2,4,4,2,2,1,3,3,1,1,1,3,3,2,4)
    laketitles <- c('Montane-Lg.', 'Montane-Sm.', 'Alpine')
    agg <- aggregate(SI[,3:ncol(SI)], by=list(lakevec), FUN=mean)
    agg2 <- aggregate(SI[,3:ncol(SI)], by=list(lakevec,SI[,2]), FUN=mean)
    gg <- agg[1:3,]
    gg <- cbind(1:3, gg)

    #patches
    # ggcolname <- colnames(gg)
    # agg2colname <- colnames(agg2)
    # trueSDs <- agg_gen()
    # gg <- cbind(gg[,1:14], trueSDs[1])
    # colnames(gg) <- ggcolname
    # agg2 <- cbind(agg2[,1:14], trueSDs[2])
    # colnames(agg2) <- agg2colname

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


        xmin <- min(means1, na.rm=T) - max(sds1, na.rm=T)
        xmax <- max(means1, na.rm=T) + max(sds1, na.rm=T)
        ymin <- min(means2, na.rm=T) - max(sds2, na.rm=T)
        ymax <- max(means2, na.rm=T) + max(sds2, na.rm=T)


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
                 main=paste(laketitles[which(1:3==lake)]))
        } else {
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1)
        }

        polygon(x[hullind,][,1], x[hullind,][,2], col='gray90', border=NA)


        print.letter(letters[count], xy=c(0.9,0.9), cex=1.3, font=2)
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

        if(count %in% c(2,5,8)){
            if(tracer1=='N'){
                mtext(expression(bold(paste(delta^15,'N'))), outer=FALSE, line=2.7, side=1)
            } else {
                if (tracer1=='C'){
                    mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=2.7, side=1)
                } else {
                    mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=2.7, side=1)
                }
            }
        }

        if(count %in% c(1,4,7)){
            if(tracer2=='N'){
                mtext(expression(bold(paste(delta^15,'N'))), outer=FALSE, line=3.1, side=2)
            } else {
                if (tracer2=='C'){
                    mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=3.1, side=2)
                } else {
                    mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=3.1, side=2)
                }
            }
        }

    }
    if(legend==TRUE){
        legend(x=-124, y=-60, legend=c('Phytoplank.', 'Periphyton', 'Terrestrial', 'Calanoida', 'Cladocera', 'Trichoptera'),
               pch=c(20,20,20,21,25,24), pt.bg='white', cex=1.5, bty='n', xpd=NA,
               col=c('steelblue1','darkgreen','sienna4','black','black','black'))

        return(list(gg, agg2))
    }

}
par(mfrow=c(3,3), mar=c(3,2,2,2), oma=c(4,4,1,7.5))
SD_data <- biplot_bylake('H', 'C', label_seed=1, titles=T, legend=T)
biplot_bylake('H', 'N', label_seed=4, titles=F, legend=F)
biplot_bylake('C', 'N', label_seed=7, titles=F, legend=F)

byspeclake <- function(){

    cal_sm <- c(1,3,7,10);  cal_lg <- c(2,6);  cal_hi <- c(4,5,8,9)
    caddis_sm <- c(20,21,30);  caddis_lg <- c(22,25,26,27);  caddis_hi <- c(23,24,28,29)
    clad_sm <- c(12,14,16,17);  clad_lg <- c(13,15)
    # cal_sm <- 2;  cal_lg <- 1;  cal_hi <- 3
    # caddis_sm <- 7;  caddis_lg <- 6;  caddis_hi <- 8
    # clad_sm <- 5;  clad_lg <- 4

    par(mfrow=c(3,3), mar=c(0,2,2,2), oma=c(5,5,2,1))
    lakes <- c('sm','lg','hi')
    laketitles <- c('Montane-Sm.', 'Montane-Lg.', 'Alpine')
    consumers <- c('cal','caddis','clad')
    constitles <- c('Calanoida','Trichoptera','Cladocera')
    count <- 0
    modes <- array(NA, dim=c(length(lakes), length(consumers), npreytypes),
                   list(lakes, consumers, c('phyto','peri','terr')))
    means <- modes
    phytocreds <- array(NA, dim=c(length(lakes), length(consumers), 2),
                        list(lakes, consumers, c('cred_mins', 'cred_maxs')))
    pericreds <- phytocreds
    terrcreds <- phytocreds

    for(cons in consumers){
        for(lake in lakes){
            count <- count + 1

            if(cons=='clad' & lake=='hi'){break}
            ind <- eval(parse(text=paste(cons,'_',lake,sep='')))

            phyto_draws <- as.vector(out[[2]][1,ind,])
            phyto <- density(phyto_draws)
            phyto_quant <- quantile(phyto_draws, probs=c(0.025,0.975))
            peri_draws <- as.vector(out[[2]][2,ind,])
            peri <- density(peri_draws)
            peri_quant <- quantile(peri_draws, probs=c(0.025,0.975))
            terr_draws <- as.vector(out[[2]][3,ind,])
            terr <- density(terr_draws)
            terr_quant <- quantile(terr_draws, probs=c(0.025,0.975))

            phyto_left_shade_x <- phyto$x[which(phyto$x <= phyto_quant[1])]
            phyto_left_shade_y <- phyto$y[which(phyto$x <= phyto_quant[1])]
            phyto_right_shade_x <- phyto$x[which(phyto$x >= phyto_quant[2])]
            phyto_right_shade_y <- phyto$y[which(phyto$x >= phyto_quant[2])]
            peri_left_shade_x <- peri$x[which(peri$x <= peri_quant[1])]
            peri_left_shade_y <- peri$y[which(peri$x <= peri_quant[1])]
            peri_right_shade_x <- peri$x[which(peri$x >= peri_quant[2])]
            peri_right_shade_y <- peri$y[which(peri$x >= peri_quant[2])]
            terr_left_shade_x <- terr$x[which(terr$x <= terr_quant[1])]
            terr_left_shade_y <- terr$y[which(terr$x <= terr_quant[1])]
            terr_right_shade_x <- terr$x[which(terr$x >= terr_quant[2])]
            terr_right_shade_y <- terr$y[which(terr$x >= terr_quant[2])]
            ymax <- max(phyto$y, peri$y, terr$y)

            plot(phyto, xlim=c(0,1), ylim=c(0,5 + ymax^(1/1.5)),
                 las=1, bty='l',
                 main='', lty=3, xaxs='i', yaxs='i', xaxt='n', col='steelblue1', lwd=2)
            # abline(v=c(0,1), lty=2)
            polygon(x=c(phyto_left_shade_x, rev(phyto_left_shade_x)),
                    y=c(phyto_left_shade_y, rep(0,length(phyto_left_shade_y))), col='steelblue1',
                    border=NA)
            polygon(x=c(phyto_right_shade_x, rev(phyto_right_shade_x)),
                    y=c(phyto_right_shade_y, rep(0,length(phyto_right_shade_y))), col='steelblue1',
                    border=NA)
            lines(peri, lty=2, col='darkgreen', lwd=2)
            polygon(x=c(peri_left_shade_x, rev(peri_left_shade_x)),
                    y=c(peri_left_shade_y, rep(0,length(peri_left_shade_y))), col='darkgreen',
                    border=NA)
            polygon(x=c(peri_right_shade_x, rev(peri_right_shade_x)),
                    y=c(peri_right_shade_y, rep(0,length(peri_right_shade_y))), col='darkgreen',
                    border=NA)
            lines(terr, lwd=3, lty=1, col='sienna4')
            polygon(x=c(terr_left_shade_x, rev(terr_left_shade_x)),
                    y=c(terr_left_shade_y, rep(0,length(terr_left_shade_y))), col='sienna4',
                    border=NA)
            polygon(x=c(terr_right_shade_x, rev(terr_right_shade_x)),
                    y=c(terr_right_shade_y, rep(0,length(terr_right_shade_y))), col='sienna4',
                    border=NA)
            if(cons=='clad' | (lake=='hi' & cons=='caddis')){
                axis(side=1)
            }
            print.letter(letters[count], xy=c(0.5,0.9), cex=1.3, font=2)
            abline(h=0)

            if(count==7){
                mtext(paste(constitles[which(consumers==cons)]), side=2, line=2.7)
            } else {
                if(count==2){
                    mtext(paste(laketitles[which(lakes==lake)]), side=3, line=1)
                } else {
                    if (count == 3){
                        mtext(paste(laketitles[which(lakes==lake)]), side=3, line=1)
                    } else {
                        if(count == 1){
                            mtext(paste(constitles[which(consumers==cons)]), side=2, line=2.7)
                            mtext(paste(laketitles[which(lakes==lake)]), side=3, line=1)
                        } else {
                            if(count==4){
                                mtext(paste(constitles[which(consumers==cons)]), side=2, line=2.7)
                            }
                        }
                    }
                }
            }

            modes[lake,cons,1] <- phyto$x[which.max(phyto$y)]
            modes[lake,cons,2] <- peri$x[which.max(peri$y)]
            modes[lake,cons,3] <- terr$x[which.max(terr$y)]
            means[lake,cons,1] <- mean(phyto_draws)
            means[lake,cons,2] <- mean(peri_draws)
            means[lake,cons,3] <- mean(terr_draws)
            terrcreds[lake,cons,1] <- terr_quant[1]
            terrcreds[lake,cons,2] <- terr_quant[2]
            pericreds[lake,cons,1] <- peri_quant[1]
            pericreds[lake,cons,2] <- peri_quant[2]
            phytocreds[lake,cons,1] <- phyto_quant[1]
            phytocreds[lake,cons,2] <- phyto_quant[2]
        }
    }
    plot(1, 1, type='n', ann=FALSE, bty='n', axes=F)
    legend(x=.5,y=1.1, legend=c('Phytoplankton', 'Periphyton', 'Terrestrial'), cex=1.5,
           lty=c(3,2,1), bty='n', lwd=c(2,2,3), col=c('steelblue1', 'darkgreen', 'sienna4'))
    mtext('Source Contribution', side=1, outer=TRUE, line=3.2, font=2)
    mtext('Density', side=2, outer=TRUE, line=3, font=2)

    return(setNames(list(modes, means, phytocreds, pericreds, terrcreds), c('modes', 'means',
                                                                            'phyto credible intervals',
                                                                            'peri credible intervals',
                                                                            'terr credible intervals')))
}
byspeclake_out <- byspeclake()

#save/load ####
# save(out, file='allochth.rda') #use this one? (it was uncommented when i revisited this)
# save(j1, file=paste0('C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/mixing_models2/',
#                 'allochth_', isos[1], isos[2], '.rda'))
load('allochth.rda')

# more plots ####
#this load and following gelman diag are broken. just skip down to the main plots.
load('allochth_allvary3.rda') #allochth_allfree 1, 2, or sensitivity_2 are for the new C:Na output
g <- matrix(NA, nrow=nvar(j1$samples), ncol=2)
for (v in 1:nvar(j1$samples)) {
    g[v,] <- gelman.diag(j1$samples[,v])$psrf
}

npreds=31
npreytypes=3


dens_lktypes <- c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1)
# dens_lktypes <- c(3,3,2,2,3,3,1,1,1,3,3,2,2,2,1,2,3,3,3,1,2,1,1,1,1,2,3,3,2,2,2)
dens_lknames <- SI[,1]
dens_pdtypes <- SI[,2]
pdfac <- factor(dens_pdtypes)
levels(pdfac) <- c(1,5,2,4,3) #this produces the correct x-axis numbers, but does not correspond to the true
#numbering scheme for preds as of the great SI-FA data-synchronization of 5/10/16.  to get those numbers, use
#the following:
# levels(pdfac) <- c(3,4,5,1,2)
order <- c(2,6,13,15,22,25,26,27,31,1,3,7,10,12,14,16,17,20,21,30,4,5,8,9,11,18,19,23,24,28,29)
#this is the numbering scheme that worked with the old SI data order
# order <- c(23,24,22,25,15,7,8,9,20,13,21,14,16,30,29,26,31,3,4,12,27,17,19,18,28,1,2,5,6,10,11)

#boxplot of allochthony (this seems to be broken; not worth much anyway)
box <- function(dens_chains){

    par(mar=c(4,4,2,9), oma=c(0,0,0,0))

    cols <- brewer.pal(3, 'Dark2')

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


    axis(side=1, at=1:31, labels=as.vector(pdfac[order]), lwd.ticks=-1, padj=-1.5, font=2)
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
box(dens_chains)

#regression stuff
reg <- function(){
    a <- brewer.pal(3, 'Dark2')
    a[2:3] <- a[c(3,2)]
    l <- a[1]
    h <- a[2]
    s <- a[3]


    for(i in 1:npreytypes){
        if(npreytypes == 3){
            maintitle <- names(ind_list2)[i]
        } else {
            maintitle <- c(which_prey[which(which_prey != 'terr')], 'terr')
        }
        plot(covs$PC1, MLEs[i,], col=c(h,h,s,s,h,h,l,l,l,h,h,s,s,s,l,s,h,h,h,l,s,l,l,l,l,s,h,h,s,s,s),
             pch=c(3,3,17,17,17,17,17,17,17,17,17,17,16,16,17,16,16,16,16,9,16,15,16,16,15,15,16,9,15,15,15),
             cex=1, lwd=2, main=maintitle[i], xlab='PC1', cex.lab=1.5, cex.main=2)
    }

    return(maintitle)
}
# plot(covs$PC1, MLEs[3,], col=c(h,h,s,s,h,h,l,l,l,h,h,s,s,s,l,s,h,h,h,l,s,l,l,l,l,s,h,h,s,s,s),
#      pch=c(3,3,17,17,17,17,17,17,17,17,17,17,16,16,17,16,16,16,16,9,16,15,16,16,15,15,16,9,15,15,15),
#      cex=1.7, lwd=2, main='Allochthony of Consumers', xlab='PC1', ylab='Terrestrial Proportion', cex.lab=1.5,
#      cex.main=2)

byspeclake <- function(){

    cal_sm <- c(1,3,7,10);  cal_lg <- c(2,6);  cal_hi <- c(4,5,8,9)
    caddis_sm <- c(20,21,30);  caddis_lg <- c(22,25,26,27);  caddis_hi <- c(23,24,28,29)
    clad_sm <- c(12,14,16,17);  clad_lg <- c(13,15)

    par(mfrow=c(3,3), mar=c(0,2,2,2), oma=c(5,5,2,1))
    lakes <- c('sm','lg','hi')
    laketitles <- c('Montane-Sm.', 'Montane-Lg.', 'Alpine')
    consumers <- c('cal','caddis','clad')
    constitles <- c('Calanoida','Trichoptera','Cladocera')
    count <- 0
    modes <- array(NA, dim=c(length(lakes), length(consumers), npreytypes),
                   list(lakes, consumers, c('phyto','peri','terr')))
    means <- modes
    phytocreds <- array(NA, dim=c(length(lakes), length(consumers), 2),
                        list(lakes, consumers, c('cred_mins', 'cred_maxs')))
    pericreds <- phytocreds
    terrcreds <- phytocreds

    for(cons in consumers){
        for(lake in lakes){
            count <- count + 1

            if(cons=='clad' & lake=='hi'){break}
            ind <- eval(parse(text=paste(cons,'_',lake,sep='')))

            phyto_draws <- as.vector(out[[2]][1,ind,])
            phyto <- density(phyto_draws)
            phyto_quant <- quantile(phyto_draws, probs=c(0.025,0.975))
            peri_draws <- as.vector(out[[2]][2,ind,])
            peri <- density(peri_draws)
            peri_quant <- quantile(peri_draws, probs=c(0.025,0.975))
            terr_draws <- as.vector(out[[2]][3,ind,])
            terr <- density(terr_draws)
            terr_quant <- quantile(terr_draws, probs=c(0.025,0.975))

            phyto_left_shade_x <- phyto$x[which(phyto$x <= phyto_quant[1])]
            phyto_left_shade_y <- phyto$y[which(phyto$x <= phyto_quant[1])]
            phyto_right_shade_x <- phyto$x[which(phyto$x >= phyto_quant[2])]
            phyto_right_shade_y <- phyto$y[which(phyto$x >= phyto_quant[2])]
            peri_left_shade_x <- peri$x[which(peri$x <= peri_quant[1])]
            peri_left_shade_y <- peri$y[which(peri$x <= peri_quant[1])]
            peri_right_shade_x <- peri$x[which(peri$x >= peri_quant[2])]
            peri_right_shade_y <- peri$y[which(peri$x >= peri_quant[2])]
            terr_left_shade_x <- terr$x[which(terr$x <= terr_quant[1])]
            terr_left_shade_y <- terr$y[which(terr$x <= terr_quant[1])]
            terr_right_shade_x <- terr$x[which(terr$x >= terr_quant[2])]
            terr_right_shade_y <- terr$y[which(terr$x >= terr_quant[2])]
            ymax <- max(phyto$y, peri$y, terr$y)

            plot(phyto, xlim=c(0,1), ylim=c(0,(5 + ymax^(1/1.5))), las=1, bty='l',
                 main='', lty=3, xaxs='i', yaxs='i', xaxt='n', col='steelblue1', lwd=2)
            polygon(x=c(phyto_left_shade_x, rev(phyto_left_shade_x)),
                    y=c(phyto_left_shade_y, rep(0,length(phyto_left_shade_y))), col='steelblue1',
                    border=NA)
            polygon(x=c(phyto_right_shade_x, rev(phyto_right_shade_x)),
                    y=c(phyto_right_shade_y, rep(0,length(phyto_right_shade_y))), col='steelblue1',
                    border=NA)
            lines(peri, lty=2, col='darkgreen', lwd=2)
            polygon(x=c(peri_left_shade_x, rev(peri_left_shade_x)),
                    y=c(peri_left_shade_y, rep(0,length(peri_left_shade_y))), col='darkgreen',
                    border=NA)
            polygon(x=c(peri_right_shade_x, rev(peri_right_shade_x)),
                    y=c(peri_right_shade_y, rep(0,length(peri_right_shade_y))), col='darkgreen',
                    border=NA)
            lines(terr, lwd=3, lty=1, col='sienna4')
            polygon(x=c(terr_left_shade_x, rev(terr_left_shade_x)),
                    y=c(terr_left_shade_y, rep(0,length(terr_left_shade_y))), col='sienna4',
                    border=NA)
            polygon(x=c(terr_right_shade_x, rev(terr_right_shade_x)),
                    y=c(terr_right_shade_y, rep(0,length(terr_right_shade_y))), col='sienna4',
                    border=NA)
            if(cons=='clad' | (lake=='hi' & cons=='caddis')){
                axis(side=1)
            }
            print.letter(letters[count], xy=c(0.5,0.9), cex=1.3, font=2)
            abline(h=0)

            if(count==7){
                mtext(paste(constitles[which(consumers==cons)]), side=2, line=2.7)
            } else {
                if(count==2){
                    mtext(paste(laketitles[which(lakes==lake)]), side=3, line=1)
                } else {
                    if (count == 3){
                        mtext(paste(laketitles[which(lakes==lake)]), side=3, line=1)
                    } else {
                        if(count == 1){
                            mtext(paste(constitles[which(consumers==cons)]), side=2, line=2.7)
                            mtext(paste(laketitles[which(lakes==lake)]), side=3, line=1)
                        } else {
                            if(count==4){
                                mtext(paste(constitles[which(consumers==cons)]), side=2, line=2.7)
                            }
                        }
                    }
                }
            }

            modes[lake,cons,1] <- phyto$x[which.max(phyto$y)]
            modes[lake,cons,2] <- peri$x[which.max(peri$y)]
            modes[lake,cons,3] <- terr$x[which.max(terr$y)]
            means[lake,cons,1] <- mean(phyto_draws)
            means[lake,cons,2] <- mean(peri_draws)
            means[lake,cons,3] <- mean(terr_draws)
            terrcreds[lake,cons,1] <- terr_quant[1]
            terrcreds[lake,cons,2] <- terr_quant[2]
            pericreds[lake,cons,1] <- peri_quant[1]
            pericreds[lake,cons,2] <- peri_quant[2]
            phytocreds[lake,cons,1] <- phyto_quant[1]
            phytocreds[lake,cons,2] <- phyto_quant[2]
        }
    }
    plot(1, 1, type='n', ann=FALSE, bty='n', axes=F)
    legend(x=.5,y=1.1, legend=c('Phytoplankton', 'Periphyton', 'Terrestrial'), cex=1.5,
           lty=c(3,2,1), bty='n', lwd=c(2,2,3), col=c('steelblue1', 'darkgreen', 'sienna4'))
    mtext('Source Contribution', side=1, outer=TRUE, line=3.2, font=2)
    mtext('Density', side=2, outer=TRUE, line=3, font=2)

    return(setNames(list(modes, means, phytocreds, pericreds, terrcreds), c('modes', 'means',
                                                                            'phyto credible intervals',
                                                                            'peri credible intervals',
                                                                            'terr credible intervals')))
}
# tiff(type='cairo', res=300, width=6, height=6, units='in',
#      filename='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/official_figures/figure_3.tif')
byspeclake_out <- byspeclake()
# dev.off()

# Mode <- function(x) {
#     ux <- unique(x)
#     ux[which.max(tabulate(match(x, ux)))]
# }
means=modes=medians=sds=1:31
for(i in 1:31){
    means[i] <- mean(out[[2]][3,i,])
    # modes[i] <- Mode(out[[2]][3,i,])
    medians[i] <- median(out[[2]][3,i,])
    sds[i] <- sd(out[[2]][3,i,])
}
write.csv(cbind(means, medians, sds), 'allochth_all_out.csv')


bylake <- function(){

    sm <- c(1,3,7,10,12,14,16,17,20,21,30)
    lg <- c(2,6,13,15,22,25,26,27,31)
    hi <- c(4,5,8,9,11,18,19,23,24,28,29)
    lakes <- c('sm','lg','hi')
    titles <- c('Small - Montane', 'Large - Montane', 'Alpine')

    modes <- array(NA, dim=c(length(lakes), npreytypes),
                   list(lakes, c('phyto','peri','terr')))
    means <- array(NA, dim=c(length(lakes), npreytypes),
                   list(lakes, c('phyto','peri','terr')))
    creds <- array(NA, dim=c(length(lakes), 2, 3),
                   list(lakes, c('cred_mins','cred_maxs'), c('phyto','peri','terr')))

    par(mfrow=c(1,3), mar=c(2,2,0,1), oma=c(3,3,3,0))
    for(lake in lakes){

        ind <- eval(parse(text=paste(lake)))

        terr_draws <- as.vector(out[[2]][3,ind,])
        terr <- density(terr_draws)
        terr_quant <- quantile(terr_draws, probs=c(0.025,0.975))
        phyto_draws <- as.vector(out[[2]][1,ind,])
        phyto <- density(phyto_draws)
        phyto_quant <- quantile(phyto_draws, probs=c(0.025,0.975))
        peri_draws <- as.vector(out[[2]][2,ind,])
        peri <- density(peri_draws)
        peri_quant <- quantile(peri_draws, probs=c(0.025,0.975))

        left_shade_x <- terr$x[which(terr$x <= terr_quant[1])]
        left_shade_y <- terr$y[which(terr$x <= terr_quant[1])]
        right_shade_x <- terr$x[which(terr$x >= terr_quant[2])]
        right_shade_y <- terr$y[which(terr$x >= terr_quant[2])]
        ymax <- max(phyto$y, peri$y, terr$y)

        plot(phyto, xlim=c(0,1), ylim=c(0,14), lwd=1, lty=3,
             main='', las=1, xlab='', ylab='', bty='l', xaxs='i', yaxs='i')
        polygon(x=c(left_shade_x, rev(left_shade_x)),
                y=c(left_shade_y, rep(0,length(left_shade_y))), col='gray50',
                border=NA)
        polygon(x=c(right_shade_x, rev(right_shade_x)),
                y=c(right_shade_y, rep(0,length(right_shade_y))), col='gray50',
                border=NA)

        lines(peri, lwd=1, lty=2)
        lines(terr, lwd=2, lty=1)
        mtext('Source Contribution', side=1, outer=TRUE, line=1)
        mtext('Density', side=2, outer=TRUE, line=1)
        # mtext(titles[which(lakes==lake)], side=3, outer=FALSE, line=-2)
        print.letter(letters[which(lakes==lake)], xy=c(0.85,0.9), cex=1.3, font=2)
        abline(h=0)

        modes[lake,1] <- phyto$x[which.max(phyto$y)]
        modes[lake,2] <- peri$x[which.max(peri$y)]
        modes[lake,3] <- terr$x[which.max(terr$y)]
        means[lake,1] <- mean(phyto_draws)
        means[lake,2] <- mean(peri_draws)
        means[lake,3] <- mean(terr_draws)
        creds[lake,1,1] <- phyto_quant[1]
        creds[lake,2,1] <- phyto_quant[2]
        creds[lake,1,2] <- peri_quant[1]
        creds[lake,2,2] <- peri_quant[2]
        creds[lake,1,3] <- terr_quant[1]
        creds[lake,2,3] <- terr_quant[2]
    }
    legend(x=-1.6, y=16.7, legend=c('Phytoplankton', 'Periphyton', 'Terrestrial'),
           lty=c(3,2,1), xpd=NA, bty='n', horiz=TRUE, cex=1.2, text.width=c(.4,.4,.35),
           lwd=c(1,1,2))

    return(setNames(list(modes, means, creds), c('modes', 'means', 'credible intervals')))
}
bylake_out <- bylake()

biplot_byspeclake <- function(tracer1, tracer2){

    lakevec <- c(2,1,2,3,3,1,2,3,3,2,4,2,1,2,1,2,2,4,4,2,2,1,3,3,1,1,1,3,3,2,4)
    agg <- aggregate(SI[,3:ncol(SI)], by=list(lakevec,SI[,2]), FUN=mean)
    gg <- agg[c(1:3,5,6,8:10),]

    inds <- setNames(c(0,4,8), c('N', 'C', 'H'))
    ind1 <- unname(inds[names(inds) == tracer1])
    ind2 <- unname(inds[names(inds) == tracer2])

    augs <- setNames(c(1,1,3), c('N', 'C', 'H'))
    aug1 <- unname(augs[names(augs) == tracer1])
    aug2 <- unname(augs[names(augs) == tracer2])

    # labs <- setNames(c("delta^15, 'N'", "delta^13, 'C'", "delta^2, 'H'"), c('N', 'C', 'H'))
    # lab1 <- unname(labs[names(labs) == tracer1])
    # lab2 <- unname(labs[names(labs) == tracer2])
    count <- 0

    par(mfrow=c(3,3), mar=c(2,2,2,2), oma=c(4,5,1,1))
    for(lake in 1:3){
        for(cons in c('Calanoida', 'Trichoptera', 'Cladocera')){

            count <- count + 1
            if(cons=='Cladocera' & lake=='3'){break}

            row <- which(gg[,1] == lake & gg[,2] == cons)
            xmin <- min(gg[row,(3:6)+ind1]) - max(gg[row,(15:18)+ind1])
            xmax <- max(gg[row,(3:6)+ind1]) + max(gg[row,(15:18)+ind1])
            ymin <- min(gg[row,(3:6)+ind2]) - max(gg[row,(15:18)+ind2])
            ymax <- max(gg[row,(3:6)+ind2]) + max(gg[row,(15:18)+ind2])

            err <- 'gray30'; pt <- 'black'

            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1)
            polygon(x=c(gg[row,4+ind1], gg[row,5+ind1], gg[row,6+ind1]), border=NA,
                    y=c(gg[row,4+ind2], gg[row,5+ind2], gg[row,6+ind2]), col='gray75')


            # polygon(x=c(gg[row,4+ind1], gg[row,5+ind1], gg[row,6+ind1]), border=NA,
            #         y=c(gg[row,4+ind2], gg[row,5+ind2], gg[row,6+ind2]), col='gray90')

            #cons
            error.bars(gg[row,3+ind1], gg[row,3+ind2], upper=gg[row,15+ind1], horiz=TRUE,
                       cap.length=0.01, col=err)
            error.bars(gg[row,3+ind1], gg[row,3+ind2], upper=gg[row,15+ind2], horiz=FALSE,
                       cap.length=0.01, col=err)
            points(gg[row,3+ind1], gg[row,3+ind2], pch=20, cex=2, col=pt)
            print.letter(letters[count], xy=c(0.95,0.95), cex=1.3, font=2)
            #phyto
            error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind1]*aug1, horiz=TRUE,
                       cap.length=0.01, col=err)
            error.bars(gg[row,4+ind1], gg[row,4+ind2], upper=gg[row,16+ind2]*aug2, horiz=FALSE,
                       cap.length=0.01, col=err)
            points(gg[row,4+ind1], gg[row,4+ind2],
                   pch=22, cex=2, col=pt, bg=pt)
            #peri
            error.bars(gg[row,5+ind1], gg[row,5+ind2], upper=gg[row,17+ind1], horiz=TRUE,
                       cap.length=0.01, col=err)
            error.bars(gg[row,5+ind1], gg[row,4+ind2], upper=gg[row,17+ind2], horiz=FALSE,
                       cap.length=0.01, col=err)
            points(gg[row,5+ind1], gg[row,5+ind2],
                   pch=23, cex=2, col=pt, bg=pt)
            #terr
            error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind1], horiz=TRUE,
                       cap.length=0.01, col=err)
            error.bars(gg[row,6+ind1], gg[row,6+ind2], upper=gg[row,18+ind2], horiz=FALSE,
                       cap.length=0.01, col=err)
            points(gg[row,6+ind1], gg[row,6+ind2],
                   pch=24, cex=2, col=pt, bg=pt)

        }

        if(tracer1=='N'){
            mtext(expression(bold(paste(delta^15,'N'))), outer=TRUE, line=1.5, side=1)
        } else {
            if (tracer1=='C'){
                mtext(expression(bold(paste(delta^13,'C'))), outer=TRUE, line=1.5, side=1)
            } else {
                mtext(expression(bold(paste(delta^2,'H'))), outer=TRUE, line=1.5, side=1)
            }
        }

        if(tracer2=='N'){
            mtext(expression(bold(paste(delta^15,'N'))), outer=TRUE, line=1.5, side=2)
        } else {
            if (tracer2=='C'){
                mtext(expression(bold(paste(delta^13,'C'))), outer=TRUE, line=1.5, side=2)
            } else {
                mtext(expression(bold(paste(delta^2,'H'))), outer=TRUE, line=1.5, side=2)
            }
        }

    }
    plot(1, 1, type='n', ann=FALSE, bty='n', axes=F)
    legend(x='center', legend=c('Consumer', 'Phytoplankton', 'Periphyton', 'Terrestrial'),
           pch=c(20,22,23,24), pt.bg='black', cex=c(2,2,2,2), bty='n')
}
biplot_byspeclake('C', 'N')


# tiff(type='cairo', res=300, width=7, height=6, units='in',
#      filename='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/official_figures/figure_2.tif')
agg_gen <- function(){
    #to get SDs
    load("C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/sensitivity_analysis/real_C.rda")
    load("C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/sensitivity_analysis/real_N.rda")
    load("C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/sensitivity_analysis/real_H.rda")

    phyto_sm_c <- sd(as.vector(allvary_C[[2]][3:6,10,]), na.rm=TRUE) * 1
    peri_sm_c <- sd(as.vector(allvary_C[[2]][3:6,11,]), na.rm=TRUE) * 1
    terr_sm_c <- sd(as.vector(allvary_C[[2]][3:6,5,]), na.rm=TRUE) * 1
    phyto_lg_c <- sd(as.vector(allvary_C[[2]][c(1,2,10),10,]), na.rm=TRUE) * 1
    peri_lg_c <- sd(as.vector(allvary_C[[2]][c(1,2,10),11,]), na.rm=TRUE) * 1
    terr_lg_c <- sd(as.vector(allvary_C[[2]][c(1,2,10),5,]), na.rm=TRUE) * 1
    phyto_hi_c <- sd(as.vector(allvary_C[[2]][c(7:9,11),10,]), na.rm=TRUE) * 1
    peri_hi_c <- sd(as.vector(allvary_C[[2]][c(7:9,11),11,]), na.rm=TRUE) * 1
    terr_hi_c <- sd(as.vector(allvary_C[[2]][c(7:9,11),5,]), na.rm=TRUE) * 1
    phyto_sm_n <- sd(as.vector(allvary_N[[2]][3:6,10,]), na.rm=TRUE) * 1
    peri_sm_n <- sd(as.vector(allvary_N[[2]][3:6,11,]), na.rm=TRUE) * 1
    terr_sm_n <- sd(as.vector(allvary_N[[2]][3:6,5,]), na.rm=TRUE) * 1
    phyto_lg_n <- sd(as.vector(allvary_N[[2]][c(1,2,10),10,]), na.rm=TRUE) * 1
    peri_lg_n <- sd(as.vector(allvary_N[[2]][c(1,2,10),11,]), na.rm=TRUE) * 1
    terr_lg_n <- sd(as.vector(allvary_N[[2]][c(1,2,10),5,]), na.rm=TRUE) * 1
    phyto_hi_n <- sd(as.vector(allvary_N[[2]][c(7:9,11),10,]), na.rm=TRUE) * 1
    peri_hi_n <- sd(as.vector(allvary_N[[2]][c(7:9,11),11,]), na.rm=TRUE) * 1
    terr_hi_n <- sd(as.vector(allvary_N[[2]][c(7:9,11),5,]), na.rm=TRUE) * 1
    phyto_sm_h <- sd(as.vector(allvary_H[[2]][3:6,4,]), na.rm=TRUE) * 1
    peri_sm_h <- sd(as.vector(allvary_H[[2]][3:6,3,]), na.rm=TRUE) * 1
    phyto_lg_h <- sd(as.vector(allvary_H[[2]][c(1,2,10),4,]), na.rm=TRUE) * 1
    peri_lg_h <- sd(as.vector(allvary_H[[2]][c(1,2,10),3,]), na.rm=TRUE) * 1
    phyto_hi_h <- sd(as.vector(allvary_H[[2]][c(7:9,11),4,]), na.rm=TRUE) * 1
    peri_hi_h <- sd(as.vector(allvary_H[[2]][c(7:9,11),3,]), na.rm=TRUE) * 1
    load("C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/sensitivity_analysis/withd2HT.rda")
    terr_sm_h <- sd(as.vector(allvary_H[[2]][3:6,2,]), na.rm=TRUE) * 1
    terr_lg_h <- sd(as.vector(allvary_H[[2]][c(1,2,10),2,]), na.rm=TRUE) * 1
    terr_hi_h <- sd(as.vector(allvary_H[[2]][c(7:9,11),2,]), na.rm=TRUE) * 1
    load("C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/sensitivity_analysis/real_H.rda")

    cal_sm_c <- sd(as.vector(allvary_C[[2]][c(1,3,7,10),9,]), na.rm=TRUE) * 1
    caddis_sm_c <- sd(as.vector(allvary_C[[2]][c(20,21,30),9,]), na.rm=TRUE) * 1
    clad_sm_c <- sd(as.vector(allvary_C[[2]][c(12,14,16,17),9,]), na.rm=TRUE) * 1
    cal_lg_c <- sd(as.vector(allvary_C[[2]][c(2,6),9,]), na.rm=TRUE) * 1
    caddis_lg_c <- sd(as.vector(allvary_C[[2]][c(25:27,31),9,]), na.rm=TRUE) * 1
    clad_lg_c <- sd(as.vector(allvary_C[[2]][c(13,15),9,]), na.rm=TRUE) * 1
    cal_hi_c <- sd(as.vector(allvary_C[[2]][c(4,5,8,9),9,]), na.rm=TRUE) * 1
    caddis_hi_c <- sd(as.vector(allvary_C[[2]][c(23,24,28,29),9,]), na.rm=TRUE) * 1
    cal_sm_n <- sd(as.vector(allvary_N[[2]][c(1,3,7,10),9,]), na.rm=TRUE) * 1
    caddis_sm_n <- sd(as.vector(allvary_N[[2]][c(20,21,30),9,]), na.rm=TRUE) * 1
    clad_sm_n <- sd(as.vector(allvary_N[[2]][c(12,14,16,17),9,]), na.rm=TRUE) * 1
    cal_lg_n <- sd(as.vector(allvary_N[[2]][c(2,6),9,]), na.rm=TRUE) * 1
    caddis_lg_n <- sd(as.vector(allvary_N[[2]][c(25:27,31),9,]), na.rm=TRUE) * 1
    clad_lg_n <- sd(as.vector(allvary_N[[2]][c(13,15),9,]), na.rm=TRUE) * 1
    cal_hi_n <- sd(as.vector(allvary_N[[2]][c(4,5,8,9),9,]), na.rm=TRUE) * 1
    caddis_hi_n <- sd(as.vector(allvary_N[[2]][c(23,24,28,29),9,]), na.rm=TRUE) * 1
    cal_sm_h <- sd(as.vector(allvary_H[[2]][c(1,3,7,10),6,]), na.rm=TRUE) * 1
    caddis_sm_h <- sd(as.vector(allvary_H[[2]][c(20,21,30),6,]), na.rm=TRUE) * 1
    clad_sm_h <- sd(as.vector(allvary_H[[2]][c(12,14,16,17),6,]), na.rm=TRUE) * 1
    cal_lg_h <- sd(as.vector(allvary_H[[2]][c(2,6),6,]), na.rm=TRUE) * 1
    caddis_lg_h <- sd(as.vector(allvary_H[[2]][c(25:27,31),6,]), na.rm=TRUE) * 1
    clad_lg_h <- sd(as.vector(allvary_H[[2]][c(13,15),6,]), na.rm=TRUE) * 1
    cal_hi_h <- sd(as.vector(allvary_H[[2]][c(4,5,8,9),6,]), na.rm=TRUE) * 1
    caddis_hi_h <- sd(as.vector(allvary_H[[2]][c(23,24,28,29),6,]), na.rm=TRUE) * 1

    out <- as.data.frame(matrix(NA, nrow=3, ncol=12))
    out[,1] <- rep(NA,3)
    out[,2] <- c(phyto_lg_n, phyto_sm_n, phyto_hi_n)
    out[,3] <- c(peri_lg_n, peri_sm_n, peri_hi_n)
    out[,4] <- c(terr_lg_n, terr_sm_n, terr_hi_n)
    out[,5] <- rep(NA,3)
    out[,6] <- c(phyto_lg_c, phyto_sm_c, phyto_hi_c)
    out[,7] <- c(peri_lg_c, peri_sm_c, peri_hi_c)
    out[,8] <- c(terr_lg_c, terr_sm_c, terr_hi_c)
    out[,9] <- rep(NA,3)
    out[,10] <- c(phyto_lg_h, phyto_sm_h, phyto_hi_h)
    out[,11] <- c(peri_lg_h, peri_sm_h, peri_hi_h)
    out[,12] <- c(terr_lg_h, terr_sm_h, terr_hi_h)

    out2 <- as.data.frame(matrix(NA, nrow=10, ncol=12))
    # colnames(out2) <- c('Group.1','Group.2','cons_N_u','phyto_N_u','peri_N_u','terr_N_u',
    #                     'cons_C_u','phyto_C_u','peri_C_u','terr_C_u','cons_H_u','phyto_H_u',
    #                     'peri_H_u','terr_H_u','cons_N_sd','phyto_N_sd','peri_N_sd','terr_N_sd',
    #                     'cons_C_sd','phyto_C_sd','peri_C_sd','terr_C_sd','cons_H_sd','phyto_H_sd',
    #                     'peri_H_sd','terr_H_sd')
    # out2[,1] <- 1:10; out2[,2] <- c('Calanoida','Calanoida','Calanoida','Chaoborus','Cladocera',
    #                                 'Cladocera','Salvelinus','Trichoptera','Trichoptera','Trichoptera')
    out2[,1] <- c(cal_lg_n,cal_sm_n,cal_hi_n,NA,clad_lg_n,clad_sm_n,NA,caddis_lg_n,caddis_sm_n,caddis_hi_n)
    out2[,2] <- c(phyto_lg_n,phyto_sm_n,phyto_hi_n,NA,phyto_lg_n,phyto_sm_n,NA,phyto_lg_n,phyto_sm_n,phyto_hi_n)
    out2[,3] <- c(peri_lg_n,peri_sm_n,peri_hi_n,NA,peri_lg_n,peri_sm_n,NA,peri_lg_n,peri_sm_n,peri_hi_n)
    out2[,4] <- c(terr_lg_n,terr_sm_n,terr_hi_n,NA,terr_lg_n,terr_sm_n,NA,terr_lg_n,terr_sm_n,terr_hi_n)
    out2[,5] <- c(cal_lg_c,cal_sm_c,cal_hi_c,NA,clad_lg_c,clad_sm_c,NA,caddis_lg_c,caddis_sm_c,caddis_hi_c)
    out2[,6] <- c(phyto_lg_c,phyto_sm_c,phyto_hi_c,NA,phyto_lg_c,phyto_sm_c,NA,phyto_lg_c,phyto_sm_c,phyto_hi_c)
    out2[,7] <- c(peri_lg_c,peri_sm_c,peri_hi_c,NA,peri_lg_c,peri_sm_c,NA,peri_lg_c,peri_sm_c,peri_hi_c)
    out2[,8] <- c(terr_lg_c,terr_sm_c,terr_hi_c,NA,terr_lg_c,terr_sm_c,NA,terr_lg_c,terr_sm_c,terr_hi_c)
    out2[,9] <- c(cal_lg_h,cal_sm_h,cal_hi_h,NA,clad_lg_h,clad_sm_h,NA,caddis_lg_h,caddis_sm_h,caddis_hi_h)
    out2[,10] <- c(phyto_lg_h,phyto_sm_h,phyto_hi_h,NA,phyto_lg_h,phyto_sm_h,NA,phyto_lg_h,phyto_sm_h,phyto_hi_h)
    out2[,11] <- c(peri_lg_h,peri_sm_h,peri_hi_h,NA,peri_lg_h,peri_sm_h,NA,peri_lg_h,peri_sm_h,peri_hi_h)
    out2[,12] <- c(terr_lg_h,terr_sm_h,terr_hi_h,NA,terr_lg_h,terr_sm_h,NA,terr_lg_h,terr_sm_h,terr_hi_h)

    return(list(out, out2))
}
# tracer1='H'; tracer2='C'; label_seed=1; titles=T; legend=T
biplot_bylake <- function(tracer1, tracer2, label_seed=1, titles=TRUE, legend=TRUE){

    lakevec <- c(2,1,2,3,3,1,2,3,3,2,4,2,1,2,1,2,2,4,4,2,2,1,3,3,1,1,1,3,3,2,4)
    laketitles <- c('Montane-Lg.', 'Montane-Sm.', 'Alpine')
    agg <- aggregate(SI[,3:ncol(SI)], by=list(lakevec), FUN=mean)
    agg2 <- aggregate(SI[,3:ncol(SI)], by=list(lakevec,SI[,2]), FUN=mean)
    gg <- agg[1:3,]
    gg <- cbind(1:3, gg)

    #patches
    # ggcolname <- colnames(gg)
    # agg2colname <- colnames(agg2)
    # trueSDs <- agg_gen()
    # gg <- cbind(gg[,1:14], trueSDs[1])
    # colnames(gg) <- ggcolname
    # agg2 <- cbind(agg2[,1:14], trueSDs[2])
    # colnames(agg2) <- agg2colname

    inds <- setNames(c(0,4,8), c('N', 'C', 'H'))
    ind1 <- unname(inds[names(inds) == tracer1])
    ind2 <- unname(inds[names(inds) == tracer2])

    augs <- setNames(c(1,1,1), c('N', 'C', 'H'))
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


        xmin <- min(means1, na.rm=T) - max(sds1, na.rm=T)
        xmax <- max(means1, na.rm=T) + max(sds1, na.rm=T)
        ymin <- min(means2, na.rm=T) - max(sds2, na.rm=T)
        ymax <- max(means2, na.rm=T) + max(sds2, na.rm=T)


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
                 main=paste(laketitles[which(1:3==lake)]))
        } else {
            plot(1,1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type='n', bty='l',las=1)
        }

        polygon(x[hullind,][,1], x[hullind,][,2], col='gray90', border=NA)


        print.letter(letters[count], xy=c(0.9,0.9), cex=1.3, font=2)
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

        if(count %in% c(2,5,8)){
            if(tracer1=='N'){
                mtext(expression(bold(paste(delta^15,'N'))), outer=FALSE, line=2.7, side=1)
            } else {
                if (tracer1=='C'){
                    mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=2.7, side=1)
                } else {
                    mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=2.7, side=1)
                }
            }
        }

        if(count %in% c(1,4,7)){
            if(tracer2=='N'){
                mtext(expression(bold(paste(delta^15,'N'))), outer=FALSE, line=3.1, side=2)
            } else {
                if (tracer2=='C'){
                    mtext(expression(bold(paste(delta^13,'C'))), outer=FALSE, line=3.1, side=2)
                } else {
                    mtext(expression(bold(paste(delta^2,'H'))), outer=FALSE, line=3.1, side=2)
                }
            }
        }
    }
    if(legend==TRUE){
        legend(x=-124, y=-48, legend=c('Phytoplank.', 'Periphyton', 'Terrestrial', 'Calanoida', 'Cladocera', 'Trichoptera'),
               pch=c(20,20,20,21,25,24), pt.bg='white', cex=1.5, bty='n', xpd=NA,
               col=c('steelblue1','darkgreen','sienna4','black','black','black'))

        return(list(gg, agg2))
    }

}
par(mfrow=c(3,3), mar=c(3,2,2,2), oma=c(4,4,1,7.5))
SD_data <- biplot_bylake('H', 'C', label_seed=1, titles=T, legend=T)
biplot_bylake('H', 'N', label_seed=4, titles=F, legend=F)
biplot_bylake('C', 'N', label_seed=7, titles=F, legend=F)
# write.csv(SD_data[[1]], "C:/Users/Mike/Desktop/correctionPhase_byclass.csv")
# write.csv(SD_data[[2]], "C:/Users/Mike/Desktop/correctionPhase.csv")
# dev.off()

#wonky density plot
wonky <- function() {
    par(mar=c(4,4,2,2),oma=c(1,1,1,1))
    #most of this code is unnecessary.  just borrowed it from above
    cal_sm <- c(1,3,7,10);  cal_lg <- c(2,6);  cal_hi <- c(4,5,8,9)
    caddis_sm <- c(20,21,30);  caddis_lg <- c(22,25,26,27);  caddis_hi <- c(23,24,28,29)
    clad_sm <- c(12,14,16,17);  clad_lg <- c(13,15)

    par(mfrow=c(1,1))#, mar=c(0,2,2,2), oma=c(5,5,2,1))
    lakes <- c('sm','lg','hi')
    laketitles <- c('Montane-Sm.', 'Montane-Lg.', 'Alpine')
    consumers <- c('cal','caddis','clad')
    constitles <- c('Calanoida','Trichoptera','Cladocera')
    count <- 0
    modes <- array(NA, dim=c(length(lakes), length(consumers), npreytypes),
                   list(lakes, consumers, c('phyto','peri','terr')))
    creds <- array(NA, dim=c(length(lakes), length(consumers), 2),
                   list(lakes, consumers, c('cred_mins', 'cred_maxs')))

    cons = 'caddis'
    lake = 'lg'

    ind <- eval(parse(text=paste(cons,'_',lake,sep='')))

    phyto_draws <- as.vector(out[[2]][1,ind,])
    phyto <- density(phyto_draws)
    phyto_quant <- quantile(phyto_draws, probs=c(0.025,0.975))
    peri_draws <- as.vector(out[[2]][2,ind,])
    peri <- density(peri_draws)
    peri_quant <- quantile(peri_draws, probs=c(0.025,0.975))
    terr_draws <- as.vector(out[[2]][3,ind,])
    terr <- density(terr_draws)
    terr_quant <- quantile(terr_draws, probs=c(0.025,0.975))

    plot(phyto_draws, peri_draws, xlab='Phytoplankton Proportion', ylab='Periphyton Proportion')
    # plot(phyto_draws, terr_draws, xlab='Phytoplankton Proportion', ylab='Terrestrial Proportion')
    # plot(peri_draws, terr_draws, xlab='Periphyton Proportion', ylab='Terrestrial Proportion')
}
wonky()

# N_means <- SI[,3:6]; C_means <- SI[,7:10]; H_means <- SI[,11:14]
# N_sds <- SI[,15:18]; C_sds <- SI[,19:22]; H_sds <- SI[,23:26]
# for(i in 1:npreds){
#
#     plot(C_means[i,1], H_means[i,1], xlim=c(min(C_means)-max(C_sds), max(C_means)+max(C_sds)),
#          ylim=c(min(H_means)-max(H_sds), max(H_means)+max(H_sds)), pch=16, cex=2,
#          main=paste(SI$lake_name[i], ' - ', SI$sample_type[i]),
#          xlab='C', ylab='H') #consumer
#     error.bars(x=C_means[i,1], y=H_means[i,1], upper=C_sds[i,1], cap.length=0.01, horiz=TRUE, lwd=1.5)
#     error.bars(x=C_means[i,1], y=H_means[i,1], upper=H_sds[i,1], cap.length=0.01, horiz=FALSE, lwd=1.5)
#     error.bars(x=C_means[i,2], y=H_means[i,2], upper=C_sds[i,2], cap.length=0.01, horiz=TRUE, lwd=1.5)
#     error.bars(x=C_means[i,2], y=H_means[i,2], upper=H_sds[i,2], cap.length=0.01, horiz=FALSE, lwd=1.5)
#     points(C_means[i,2], H_means[i,2], pch=16, cex=2, col='lightgreen') #phyto
#     error.bars(x=C_means[i,3], y=H_means[i,3], upper=C_sds[i,3], cap.length=0.01, horiz=TRUE, lwd=1.5)
#     error.bars(x=C_means[i,3], y=H_means[i,3], upper=H_sds[i,3], cap.length=0.01, horiz=FALSE, lwd=1.5)
#     points(C_means[i,3], H_means[i,3], pch=16, cex=2, col='brown') #peri
#     error.bars(x=C_means[i,4], y=H_means[i,4], upper=C_sds[i,4], cap.length=0.01, horiz=TRUE, lwd=1.5)
#     error.bars(x=C_means[i,4], y=H_means[i,4], upper=H_sds[i,4], cap.length=0.01, horiz=FALSE, lwd=1.5)
#     points(C_means[i,4], H_means[i,4], pch=16, cex=2, col='darkgreen') #terr
# }
# plot(1, 1, type='n')
# legend(x='center', legend=c('cons', 'phyto', 'peri', 'terr'),
#        pch=16, cex=1, col=c('black','lightgreen','brown','darkgreen'))


#posteriors and DIC + box and regression
tracer='SI'
plotter <- function(MLE_arr, dens_arr, tracer){
    #combined posteriors for consumer species x lake
    byspeclake()
    #and by lake alone
    bylake()

    #regression plots
    preynames <- reg()

    #posteriors
    par(mfrow=c(4,4))
    for(i in 1:npreds){

        if(npreytypes == 3){
            max_y <- max(c(density(dens_arr[1,i,], n=10000)$y, density(dens_arr[2,i,], n=10000)$y,
                           density(dens_arr[3,i,], n=10000)$y))
        } else {
            max_y <- max(c(density(dens_arr[1,i,], n=10000)$y, density(dens_arr[2,i,], n=10000)$y))
        }

        if(tracer == 'SI'){
            title <- paste(SI$lake_name[i], ' - ', SI$sample_type[i])
        } else {
            if(tracer == 'FA'){
                title <- paste(lakesByPred[i], ' - ', predsvec[i])
            } else {
                print('oh, shitballs. gotta work out SI-FA combo')
            }
        }

        if(npreytypes == 3){
            plot(density(dens_arr[1,i,], n=10000), col='darkred', type='l', main=title,
                 xlim=c(0,1), ylim=c(0,max_y * 1.1))
            lines(density(dens_arr[2,i,], n=10000), col='darkblue')
            lines(density(dens_arr[3,i,], n=10000), col='darkgreen')
            text(x=MLE_arr[1,i], y=max_y * 1.05, labels=round(MLE_arr[1,i],2), col='darkred')
            text(x=MLE_arr[2,i], y=max_y * 1.05, labels=round(MLE_arr[2,i],2), col='darkblue')
            text(x=MLE_arr[3,i], y=max_y * 1.05, labels=round(MLE_arr[3,i],2), col='darkgreen')
        } else {
            if(preynames[1] == 'phyto'){
                plot(density(dens_arr[1,i,], n=10000), col='darkred', type='l', main=title,
                     xlim=c(0,1), ylim=c(0,max_y * 1.1))
                lines(density(dens_arr[2,i,], n=10000), col='darkgreen')
                text(x=MLE_arr[1,i], y=max_y * 1.05, labels=round(MLE_arr[1,i],2), col='darkred')
                text(x=MLE_arr[2,i], y=max_y * 1.05, labels=round(MLE_arr[2,i],2), col='darkgreen')
            } else {
                plot(density(dens_arr[1,i,], n=10000), col='darkblue', type='l', main=title,
                     xlim=c(0,1), ylim=c(0,max_y * 1.1))
                lines(density(dens_arr[2,i,], n=10000), col='darkgreen')
                text(x=MLE_arr[1,i], y=max_y * 1.05, labels=round(MLE_arr[1,i],2), col='darkblue')
                text(x=MLE_arr[2,i], y=max_y * 1.05, labels=round(MLE_arr[2,i],2), col='darkgreen')
            }
        }
    }
    plot(1, 1, type='n')
    legend(x='center', legend=c('phyto', 'peri', 'terr'),
           lty=1,
           col=c('darkred', 'darkblue', 'darkgreen'))
    par(ask=F,mfrow=c(4,4))

    #then does biplots (C and H can be changed - just haven't updated the variable names)
    N_means <- SI[,3:6]; C_means <- SI[,7:10]; H_means <- SI[,11:14]
    N_sds <- SI[,15:18]; C_sds <- SI[,19:22]; H_sds <- SI[,23:26]
    for(i in 1:npreds){

        plot(C_means[i,1], H_means[i,1], xlim=c(min(C_means)-max(C_sds), max(C_means)+max(C_sds)),
             ylim=c(min(H_means)-max(H_sds), max(H_means)+max(H_sds)), pch=16, cex=2,
             main=paste(SI$lake_name[i], ' - ', SI$sample_type[i]),
             xlab='C', ylab='H') #consumer
        error.bars(x=C_means[i,1], y=H_means[i,1], upper=C_sds[i,1], cap.length=0.01, horiz=TRUE, lwd=1.5)
        error.bars(x=C_means[i,1], y=H_means[i,1], upper=H_sds[i,1], cap.length=0.01, horiz=FALSE, lwd=1.5)
        error.bars(x=C_means[i,2], y=H_means[i,2], upper=C_sds[i,2], cap.length=0.01, horiz=TRUE, lwd=1.5)
        error.bars(x=C_means[i,2], y=H_means[i,2], upper=H_sds[i,2], cap.length=0.01, horiz=FALSE, lwd=1.5)
        points(C_means[i,2], H_means[i,2], pch=16, cex=2, col='lightgreen') #phyto
        error.bars(x=C_means[i,3], y=H_means[i,3], upper=C_sds[i,3], cap.length=0.01, horiz=TRUE, lwd=1.5)
        error.bars(x=C_means[i,3], y=H_means[i,3], upper=H_sds[i,3], cap.length=0.01, horiz=FALSE, lwd=1.5)
        points(C_means[i,3], H_means[i,3], pch=16, cex=2, col='brown') #peri
        error.bars(x=C_means[i,4], y=H_means[i,4], upper=C_sds[i,4], cap.length=0.01, horiz=TRUE, lwd=1.5)
        error.bars(x=C_means[i,4], y=H_means[i,4], upper=H_sds[i,4], cap.length=0.01, horiz=FALSE, lwd=1.5)
        points(C_means[i,4], H_means[i,4], pch=16, cex=2, col='darkgreen') #terr
    }
    plot(1, 1, type='n')
    legend(x='center', legend=c('cons', 'phyto', 'peri', 'terr'),
           pch=16, cex=1, col=c('black','lightgreen','brown','darkgreen'))

    #then plots the DIC stats
    par(mfrow=c(1,1))
    plot(1, 1, type='n')
    text(1, 1.15, labels=(paste('mean deviance = ', round(DIC$deviance, 2))))
    text(1, 1, labels=(paste('penalty = ', round(DIC$penalty, 2))))
    text(1, 0.85, labels=(paste('penalized deviance = ',
                                round(DIC$deviance + DIC$penalty, 2))))
    #then does regression plot (visually imprecise)
    #     par(mfrow=c(1,1))
    #     temp <- j1[[1]][,1][1:length(j1[[1]][,1])]
    #     int <- density(temp, n=10000)$x[which.max(density(temp, n=10000)$y)]
    #     temp <- j1[[1]][,2][1:length(j1[[1]][,2])]
    #     slope <- density(temp, n=10000)$x[which.max(density(temp, n=10000)$y)]
    #     plot(MLE_arr[npreytypes,], y, xlab='clr allochthony', ylab='abs 440')
    #     abline(a=int, b=slope)

    #boxplot of terrestrial fraction
    box(dens_chains)

    #traces and densities
    plot(j1[[1]])
}

addendum <- 'pericor_allvary_regIndiv_PC1test'

model='mod_allochth.txt'
pdf(file=paste("sensitivity_2/", str_match(model, "mod_(.*)\\.txt")[2],
               "_", addendum, "_out.pdf", sep=''), onefile=TRUE)
# pdf(file=paste("C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\",
#                "Seven Lake
#                "model_outputs\\", str_match(model, "mod_(.*)\\.txt")[2],
#                "_", addendum, "_out.pdf", sep=''), onefile=TRUE)
# pdf(file=paste("~/Desktop/grad/Projects/Thesis/",
#                "Seven Lakes Project 2014/Analysis/mixing_models/",
#                "model_outputs/", str_match(model, "mod_(.*)/.txt")[2],
#                "_", addendum, "_out.pdf", sep=''), onefile=TRUE)

plotter(MLEs, dens_chains, 'SI')
dev.off()
# write.csv(cbind(MLEs[3,], SI$lake_name, SI$sample_type), file="C:\\Users\\Mike\\Desktop\\allochth_fixedvar_bypred.csv")
