#Michael Vlah
#University of Washington
#vlahm13@gmail.com
#last edit: 7/15/16

#check data list, inits list, and output names before running!

rm(list=ls(all=TRUE))

# setup####
nlakes <- 11
nsamp <- 31
group_vector_1 <- c(1,1,2,2,2,2,3,3,3,1,3)
group_vector_2 <- c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1)
species_vector <- c(rep(1,10),2,rep(3,6),4,4,rep(5,11),2)
nspecies <- length(levels(factor(species_vector)))
ngroups <- length(levels(factor(group_vector_1)))
npreds=5
npreytypes=3

# if (is.null(dev.list()) == FALSE){dev.off()}
# windows(record=T) #open plot window (on Windows only) --use quartz() or x11() for Mac
library(rjags)
library(RColorBrewer)
library(car)
setwd("C:\\Users\\Mike\\git\\seven_lakes\\Analysis\\mixing_models\\mixing_models2")
# setwd("~/git/seven_lakes/stable_isotope_mixing_models/mixing_models2")
H_data <- read.csv("raw_H.csv")
C_data <- read.csv("raw_C.csv")
N_data <- read.csv("raw_N.csv")
#covariates <- read.csv("C:/Users/Mike/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Analysis/7lmeans.csv")
#iso <- read.csv("C:/Users/Mike/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Analysis/Isomeans.csv")
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

merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  order.by.id...and.remove.it <- function(DATA)
  {
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")

    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }
  if(!missing(keep_order))
  {
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
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

mcmc_extractor <- function(mcmc_list){ #return posteriors, modes, and combined posteriors across samples
    nchains <- length(mcmc_list)
    chain_length <- length(mcmc_list[[1]][,1])
    param_list <- gsub(" *\\[.*?\\] *", '', colnames(mcmc_list[[1]]))
    params <- levels(factor(param_list))
    nparams <- length(params)
    param_lengths <- vector('numeric', length=nparams)
    for(i in 1:nparams){
        param_lengths[i] <- length(which(param_list==paste(params[i])))
    }
    longest_param <- max(param_lengths)

    #create arrays of mode, density, and combined density (for each lake) across all chains
    temp <- matrix(data=NA, nrow=chain_length, ncol=nchains)
    mode <- matrix(data=NA, nrow=longest_param, ncol=nparams)
    dens <- array(data=NA, dim=c(longest_param, nparams, chain_length*3))
    combdens <- matrix(data=NA, nrow=chain_length*3*longest_param, ncol=nparams)
    for(j in 1:nparams){
        indices <- which(param_list==params[j])
        for(i in 1:length(indices)){
            for(k in 1:nchains){
                temp[,k] <- mcmc_list[[k]][,indices[i]][1:length(mcmc_list[[k]][,1])]
            }
            dens[i,j,] <- as.vector(temp)
            mode[i,j] <- density(dens[i,j,], n=10000)$x[which.max(density(dens[i,j,], n=10000)$y)]
        }
        combdens[,j] <- as.vector(dens[,j,])
    }
    colnames(mode) <- params
    colnames(combdens) <- params

    return(list(mode, dens, combdens))
}

# source correction ####

#initialize the model
inits1s <- list('mean_cnA_group'=rep(6.7, ngroups), 'mean_d13CA_group'=rep(-10, ngroups),
               'mean_d13Cperi_group'=rep(-22.1, ngroups), 'sd_d13Ca'=rep(0.1, nlakes),'sd_d13Cperi'=rep(.1, nlakes), 'cnT'=rep(60,nlakes),
               'd13CT'=rep(-30,nlakes), 'mean_d15NA_group'=rep(-1,ngroups), 'mean_d15Nperi_group'=rep(-1,ngroups),
               'sd_d15Na'=rep(0.1, nlakes),'sd_d15Nperi'=rep(0.1, nlakes), 'd15NT'=rep(-1,nlakes), 'mean_eH_lake'=rep(-150,nlakes),
               'mean_d2HA_group'=rep(-300, ngroups), 'mean_d2Hperi_group'=rep(-250,ngroups),
               'd2HT'=rep(-150,nlakes))
               # 'xi.1' = rnorm(nlakes, 0, 1/4^2), 'tau.eta.1' = rgamma(nlakes, shape=0.5, scale=1/0.5),
               # 'xi.2' = rnorm(nlakes, 0, 1/4^2), 'tau.eta.2' = rgamma(nlakes, shape=0.5, scale=1/0.5),
               # 'xi.3' = rnorm(nlakes, 0, 1/30^2), 'tau.eta.3' = rgamma(nlakes, shape=0.5, scale=1/0.5))
inits2s <- list('mean_cnA_group'=rep(3, ngroups), 'mean_d13CA_group'=rep(-20, ngroups),
               'mean_d13Cperi_group'=rep(-5, ngroups), 'sd_d13Ca'=rep(0.1, nlakes),'sd_d13Cperi'=rep(1, nlakes), 'cnT'=rep(25,nlakes),
               'd13CT'=rep(-25,nlakes), 'mean_d15NA_group'=rep(-5,ngroups), 'mean_d15Nperi_group'=rep(-5,ngroups),
               'sd_d15Na'=rep(0.1, nlakes),'sd_d15Nperi'=rep(1, nlakes), 'd15NT'=rep(-5,nlakes), 'mean_eH_lake'=rep(-200,nlakes),
               'mean_d2HA_group'=rep(-200, ngroups), 'mean_d2Hperi_group'=rep(-100,ngroups),
               'd2HT'=rep(0,nlakes))
               # 'xi.1' = rnorm(nlakes, 0, 1/4^2), 'tau.eta.1' = rgamma(nlakes, shape=0.5, scale=1/0.5),
               # 'xi.2' = rnorm(nlakes, 0, 1/4^2), 'tau.eta.2' = rgamma(nlakes, shape=0.5, scale=1/0.5),
               # 'xi.3' = rnorm(nlakes, 0, 1/30^2), 'tau.eta.3' = rgamma(nlakes, shape=0.5, scale=1/0.5))
inits3s <- list('mean_cnA_group'=rep(8, ngroups), 'mean_d13CA_group'=rep(-35, ngroups),
               'mean_d13Cperi_group'=rep(-40, ngroups), 'sd_d13Ca'=rep(0.1, nlakes),'sd_d13Cperi'=rep(10, nlakes), 'cnT'=rep(95,nlakes),
               'd13CT'=rep(-20,nlakes), 'mean_d15NA_group'=rep(-10,ngroups), 'mean_d15Nperi_group'=rep(-10,ngroups),
               'sd_d15Na'=rep(0.1, nlakes),'sd_d15Nperi'=rep(5, nlakes), 'd15NT'=rep(-10,nlakes), 'mean_eH_lake'=rep(-100,nlakes),
               'mean_d2HA_group'=rep(-50, ngroups), 'mean_d2Hperi_group'=rep(0,ngroups),
               'd2HT'=rep(-300,nlakes))
               # 'xi.1' = rnorm(nlakes, 0, 1/4^2), 'tau.eta.1' = rgamma(nlakes, shape=0.5, scale=1/0.5),
               # 'xi.2' = rnorm(nlakes, 0, 1/4^2), 'tau.eta.2' = rgamma(nlakes, shape=0.5, scale=1/0.5),
               # 'xi.3' = rnorm(nlakes, 0, 1/30^2), 'tau.eta.3' = rgamma(nlakes, shape=0.5, scale=1/0.5))

mod_source_cor <- jags.model('mod_correct_sources.txt',
                              data = list('nlakes' = nlakes, 'cn_terr' = C_data[1:11,3],
                                          'cnPOM' = C_data[1:11,2], 'd13C_periraw_data' = C_data[1:11,9],
                                          'cnperi' = C_data[1:11,8], 'd13C_POM_data' = C_data[1:11,4],
                                          'grps1' = group_vector_1, 'ngroups' = ngroups,
                                          'cn_terr_sd' = C_data[1:11,10], 'd13CT_sd' = C_data[1:11,11],

                                          'prec_cnA_grp' = C_data[1:3,13], 'prec_d13CA_grp' = C_data[1:3,15],
                                          'prec_d13Cperi_grp' = C_data[1:3,16],

                                          'd15N_periraw_data' = N_data[1:11,9],
                                          'd15N_POM_data' = N_data[1:11,4], 'd15NT_sd' = N_data[1:11,11],
                                          'prec_d15NA_grp' = N_data[1:3,15], 'prec_d15Nperi_grp' = N_data[1:3,16],

                                          'prec_d2Hperi_grp' = H_data[1:3,20], 'd2H_periraw_data' = H_data[1:11,14],
                                          'd2H2O_data' = H_data[1:11,3], 'd2HT_sd' = H_data[1:11,15],

                                          # 'd2HTr' = H_data[1:11,13], 'd13CTr' = C_data[1:11,5],
                                          # 'd15NTr' = N_data[1:11,5]),
                                          'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                          'd15N_terr' = N_data[1:11,5]),
                              # 'd13C_enrich'=rep(.39, 31)),
                              inits = list(inits1s, inits2s, inits3s),
                              n.chains = 3,
                              n.adapt = 3000)

#burnin
update(mod_source_cor, n.iter=5e3)

#run
source_cor <- coda.samples(mod_source_cor, c('mean_d13CA_lake','mean_d13Cperi_lake',
                                             'mean_d15NA_lake','mean_d15Nperi_lake',
                                             'mean_d2HA_lake', 'mean_d2Hperi_lake', 'mean_eH_lake',
                                             'mean_cnA_lake', 'sd_d13Ca', 'sd_d15Na',
                                             'sd_d2H2O', 'sd_d13Cperi', 'sd_d15Nperi', 'sd_d2Hperi',
                                             'cnT', 'd13CT', 'd15NT', 'd2HT'),
                                             # 'sd_d13Ct', 'sd_d15Nt', 'sd_d2Ht'),
                              thin=100, n.iter=1e5)

#convergence
# gelman.plot(source_cor) #shouldn't vary from 1 by more than .05
g <- matrix(NA, nrow=nvar(source_cor), ncol=2)
for (v in 1:nvar(source_cor)) {
    g[v,] <- gelman.diag(source_cor[,v])$psrf
}; min(g); max(g)
which(g[,1] > 1.1 | g[,2] > 1.1)
g[g[,1] > 1.1 | g[,2] > 1.1,]
# windows(record=T)
# plot(cons_cor)

#save
save(source_cor, file='cor_sources.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], 'cor_sources.csv')

# consumer correction ####

#initialize the model
#for hierarchical (obsolete)
# inits1c <- list('C_enrich_spec' = rep(-1, nspecies), 'N_enrich_spec' = rep(-1, nspecies),
#                 'env_water_spec' = rep(0.01, nspecies), 'C_true_spec' = rep(-100, nspecies),
#                 'N_true_spec' = rep(-100, nspecies), 'H_true_spec' = rep(-500, nspecies))
#                 # 'trophic_lvl' = rep(1, nsamp))
# inits2c <- list('C_enrich_spec' = rep(1, nspecies), 'N_enrich_spec' = rep(1, nspecies),
#                 'env_water_spec' = rep(0.5, nspecies), 'C_true_spec' = rep(0, nspecies),
#                 'N_true_spec' = rep(0, nspecies), 'H_true_spec' = rep(0, nspecies))
#                 # 'trophic_lvl' = rep(3, nsamp))
# inits3c <- list('C_enrich_spec' = rep(5, nspecies), 'N_enrich_spec' = rep(5, nspecies),
#                 'env_water_spec' = rep(.99, nspecies), 'C_true_spec' = rep(50, nspecies),
#                 'N_true_spec' = rep(50, nspecies), 'H_true_spec' = rep(100, nspecies))
#                 # 'trophic_lvl' = rep(5, nsamp))

inits1c <- list('C_enrich_ind' = rep(-1, nsamp), 'N_enrich_ind' = rep(-1, nsamp),
                'env_water_ind' = rep(0.01, nsamp), 'C_true_ind' = rep(-100, nsamp),
                'N_true_ind' = rep(-100, nsamp), 'H_true_ind' = rep(-500, nsamp),
                # 'sd_C' = rgamma(nsamp,1,1), 'sd_N' = rgamma(nsamp,1,1),
                # 'sd_H' = rgamma(nsamp,1,1))
                # 'sd_C' = rgamma(nsamp,shape=.1,scale=1/.1), 'sd_N' = rgamma(nsamp,shape=.1,scale=1/.1),
                # 'sd_H' = rgamma(nsamp,shape=.1,scale=1/.1))
                'xi.1' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.1' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                'xi.2' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.2' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                'xi.3' = rnorm(nsamp, 0, 1/30^2), 'tau.eta.3' = rgamma(nsamp, shape=0.5, scale=1/0.5))
                # 'trophic_lvl' = rep(1, nsamp))
inits2c <- list('C_enrich_ind' = rep(1, nsamp), 'N_enrich_ind' = rep(1, nsamp),
                'env_water_ind' = rep(0.5, nsamp), 'C_true_ind' = rep(0, nsamp),
                'N_true_ind' = rep(0, nsamp), 'H_true_ind' = rep(0, nsamp),
                # 'sd_C' = rgamma(nsamp,1,1), 'sd_N' = rgamma(nsamp,1,1),
                # 'sd_H' = rgamma(nsamp,1,1))
                # 'sd_C' = rgamma(nsamp,shape=.1,scale=1/.1), 'sd_N' = rgamma(nsamp,shape=.1,scale=1/.1),
                # 'sd_H' = rgamma(nsamp,shape=.1,scale=1/.1))
                'xi.1' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.1' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                'xi.2' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.2' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                'xi.3' = rnorm(nsamp, 0, 1/30^2), 'tau.eta.3' = rgamma(nsamp, shape=0.5, scale=1/0.5))
                # 'trophic_lvl' = rep(3, nsamp))
inits3c <- list('C_enrich_ind' = rep(5, nsamp), 'N_enrich_ind' = rep(5, nsamp),
                'env_water_ind' = rep(.99, nsamp), 'C_true_ind' = rep(50, nsamp),
                'N_true_ind' = rep(50, nsamp), 'H_true_ind' = rep(100, nsamp),
                # 'sd_C' = rgamma(nsamp,1,1), 'sd_N' = rgamma(nsamp,1,1),
                # 'sd_H' = rgamma(nsamp,1,1))
                # 'sd_C' = rgamma(nsamp,shape=.1,scale=1/.1), 'sd_N' = rgamma(nsamp,shape=.1,scale=1/.1),
                # 'sd_H' = rgamma(nsamp,shape=.1,scale=1/.1))
                'xi.1' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.1' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                'xi.2' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.2' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                'xi.3' = rnorm(nsamp, 0, 1/30^2), 'tau.eta.3' = rgamma(nsamp, shape=0.5, scale=1/0.5))
                # 'trophic_lvl' = rep(5, nsamp))


mod_cons_cor <- jags.model('mod_correct_consumers.txt',
                              data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                          'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                          'H_data' = H_data[1:31,7],
                                          'trophic_lvl' = H_data[1:31,8]),
                              inits = list(inits1c, inits2c, inits3c),
                              n.chains = 3,
                              n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                              thin=100, n.iter=1e5)

#convergence
# gelman.plot(cons_cor) #shouldn't vary from 1 by more than .05
g <- matrix(NA, nrow=nvar(cons_cor), ncol=2)
for (v in 1:nvar(cons_cor)) {
    g[v,] <- gelman.diag(cons_cor[,v])$psrf
}; min(g); max(g)
which(g[,1] > 1.1 | g[,2] > 1.1)
g[g[,1] > 1.1 | g[,2] > 1.1,]
# windows(record=T)
# plot(cons_cor)

#save
save(cons_cor, file='cor_consumers.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], 'cor_consumers.csv')
