#Michael Vlah
#University of Washington
#vlahm13@gmail.com

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

library(rjags)
# library(RColorBrewer)
# library(car)
setwd("/home/mike/git/seven_lakes/final/SI/mixing_models")
H_data <- read.csv("raw_H.csv")
C_data <- read.csv("raw_C.csv")
N_data <- read.csv("raw_N.csv")
# source_names <- as.character(H_data[1:11,1])
# cons_names <- as.character(H_data[,5])

inits <- list('xi.1' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.1' = rgamma(nsamp, shape=0.5, scale=1/0.5),
              'xi.2' = rnorm(nsamp, 0, 1/4^2), 'tau.eta.2' = rgamma(nsamp, shape=0.5, scale=1/0.5),
              'xi.3' = rnorm(nsamp, 0, 1/30^2), 'tau.eta.3' = rgamma(nsamp, shape=0.5, scale=1/0.5))

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

print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}

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

# low phytosynthetic discrimination (eH) ####
#initialize model
mod_source_cor <- jags.model('sensitivity/sources_fixedeH.txt',
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

                                         'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                         'd15N_terr' = N_data[1:11,5],
                                         'fixed_eH'= rep(-204, nlakes)),
                             #'fixed_eH'= rep(-150, nlakes)),
                             #'fixed_eH'= rep(-96, nlakes)),
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
                           thin=100, n.iter=1e5)

#save
save(source_cor, file='sensitivity/source_loweH.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], file="sensitivity/source_loweH.csv")

# ave phytosynthetic discrimination (eH) ####
#initialize model
mod_source_cor <- jags.model('sensitivity/sources_fixedeH.txt',
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

                                         'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                         'd15N_terr' = N_data[1:11,5],
                                         # 'fixed_eH'= rep(-204, nlakes)),
                                         'fixed_eH'= rep(-150, nlakes)),
                             #'fixed_eH'= rep(-96, nlakes)),
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
                           thin=100, n.iter=1e5)

#save
save(source_cor, file='sensitivity/source_aveeH.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], file="sensitivity/source_aveeH.csv")

# hi phytosynthetic discrimination (eH) ####
#initialize model
mod_source_cor <- jags.model('sensitivity/sources_fixedeH.txt',
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

                                         'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                         'd15N_terr' = N_data[1:11,5],
                                         # 'fixed_eH'= rep(-204, nlakes)),
                                         #'fixed_eH'= rep(-150, nlakes)),
                                         'fixed_eH'= rep(-96, nlakes)),
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
                           thin=100, n.iter=1e5)

#save
save(source_cor, file='sensitivity/source_hieH.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], file="sensitivity/source_hieH.csv")

# low C:N of phytoplankton (CNa) ####
#initialize model
mod_source_cor <- jags.model('sensitivity/sources_fixedCNa.txt',
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

                                         'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                         'd15N_terr' = N_data[1:11,5],
                                         'fixed_CNa'=rep(3.7, nlakes)),
                             # 'fixed_CNa'=rep(6.7, nlakes)),
                             # 'fixed_CNa'=rep(9.7, nlakes)),
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
                           thin=100, n.iter=1e5)

#save
save(source_cor, file='sensitivity/source_lowCNa.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], file="sensitivity/source_lowCNa.csv")

# ave C:N of phytoplankton (CNa) ####
#initialize model
mod_source_cor <- jags.model('sensitivity/sources_fixedCNa.txt',
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

                                         'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                         'd15N_terr' = N_data[1:11,5],
                                         # 'fixed_CNa'=rep(3.7, nlakes)),
                                         'fixed_CNa'=rep(6.7, nlakes)),
                             # 'fixed_CNa'=rep(9.7, nlakes)),
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
                           thin=100, n.iter=1e5)

#save
save(source_cor, file='sensitivity/source_aveCNa.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], file="sensitivity/source_aveCNa.csv")

# hi C:N of phytoplankton (CNa) ####
#initialize model
mod_source_cor <- jags.model('sensitivity/sources_fixedCNa.txt',
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

                                         'd2H_terr' = H_data[1:11,13], 'd13C_terr' = C_data[1:11,5],
                                         'd15N_terr' = N_data[1:11,5],
                                         # 'fixed_CNa'=rep(3.7, nlakes)),
                                         # 'fixed_CNa'=rep(6.7, nlakes)),
                                         'fixed_CNa'=rep(9.7, nlakes)),
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
                           thin=100, n.iter=1e5)

#save
save(source_cor, file='sensitivity/source_hiCNa.rda')
source_ext <- mcmc_extractor(source_cor)
write.csv(source_ext[[1]], file="sensitivity/source_hiCNa.csv")

# low trophic fractionation of carbon (DC) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedDC.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       'fixed_DC' = rep(-1.89, nsamp)),
                           # 'fixed_DC' = rep(.39, nsamp)),
                           # 'fixed_DC' = rep(2.67, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_lowDC.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_lowDC.csv")

# ave trophic fractionation of carbon (DC) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedDC.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       # 'fixed_DC' = rep(-1.89, nsamp)),
                                       'fixed_DC' = rep(.39, nsamp)),
                           # 'fixed_DC' = rep(2.67, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_aveDC.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_aveDC.csv")

# hi trophic fractionation of carbon (DC) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedDC.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       # 'fixed_DC' = rep(-1.89, nsamp)),
                                       # 'fixed_DC' = rep(.39, nsamp)),
                                       'fixed_DC' = rep(2.67, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_hiDC.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_hiDC.csv")

# low trophic fractionation of nitrogen (DN) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedDN.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       'fixed_DN' = rep(1.42, nsamp)),
                           # 'fixed_DN' = rep(3.4, nsamp)),
                           # 'fixed_DN' = rep(5.38, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_lowDN.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_lowDN.csv")

# ave trophic fractionation of nitrogen (DN) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedDN.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       # 'fixed_DN' = rep(1.42, nsamp)),
                                       'fixed_DN' = rep(3.4, nsamp)),
                           # 'fixed_DN' = rep(5.38, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_aveDN.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_aveDN.csv")

# hi trophic fractionation of nitrogen (DN) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedDN.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       # 'fixed_DN' = rep(1.42, nsamp)),
                                       # 'fixed_DN' = rep(3.4, nsamp)),
                                       'fixed_DN' = rep(5.38, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_hiDN.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_hiDN.csv")

# low environmental water (w) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedw.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       'fixed_w' = rep(.04, nsamp)),
                           # 'fixed_w' = rep(.28, nsamp)),
                           # 'fixed_w' = rep(.52, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_loww.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_loww.csv")

# ave environmental water (w) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedw.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       # 'fixed_w' = rep(.04, nsamp)),
                                       'fixed_w' = rep(.28, nsamp)),
                           # 'fixed_w' = rep(.52, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_avew.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_avew.csv")

# hi environmental water (w) ####
#initialize model
mod_cons_cor <- jags.model('sensitivity/consumers_fixedw.txt',
                           data = list('nsamp' = nsamp, 'd2H2O' = H_data[1:31,9],
                                       'C_data' = C_data[1:31,6], 'N_data' = N_data[1:31,6],
                                       'H_data' = H_data[1:31,7],
                                       'trophic_lvl' = H_data[1:31,8],
                                       # 'fixed_w' = rep(.04, nsamp)),
                                       # 'fixed_w' = rep(.28, nsamp)),
                                       'fixed_w' = rep(.52, nsamp)),
                           inits = inits,
                           n.chains = 3,
                           n.adapt = 3000)

#burnin
update(mod_cons_cor, n.iter=5e3)

#run
cons_cor <- coda.samples(mod_cons_cor, c('C_true_ind', 'N_true_ind', 'H_true_ind', 'env_water_tot',
                                         'trophic_lvl', 'C_enrich_ind', 'N_enrich_ind', 'env_water_ind',
                                         'sd_C', 'sd_N', 'sd_H'),
                         thin=100, n.iter=1e5)

#save
save(cons_cor, file='sensitivity/cons_hiw.rda')
cons_ext <- mcmc_extractor(cons_cor)
write.csv(cons_ext[[1]], file="sensitivity/cons_hiw.csv")

# write combined scripts####

for(i in c('low','ave','hi')){
    for(j in c('eH','CNa')){

        source <- read.csv(paste0('sensitivity/source_', i, j, '.csv'))
        cons <- read.csv('cor_consumers.csv')

        consumers <- as.data.frame(cbind(as.character(H_data$lake_name), as.character(H_data$sample_type),
                                         cons$N_true_ind, cons$C_true_ind,
                                         cons$H_true_ind, cons$sd_N, cons$sd_C, cons$sd_H))
        colnames(consumers) <- c('lake_name','sample_type','cons_N_u','cons_C_u','cons_H_u','cons_N_sd','cons_C_sd','cons_H_sd')

        sources <- as.data.frame(cbind(as.character(H_data$Lake[1:11]), source$mean_d15NA_lake, source$mean_d15Nperi_lake,
                                       N_data$terr_d15N[1:11], source$mean_d13CA_lake, source$mean_d13Cperi_lake,
                                       C_data$terr_d13C[1:11], source$mean_d2HA_lake, source$mean_d2Hperi_lake,
                                       H_data$dD_terr_mean[1:11], source$sd_d15Na, source$sd_d15Nperi, N_data$d15NT_sd[1:11],
                                       source$sd_d13Ca, source$sd_d13Cperi, C_data$d13CT_sd[1:11], source$sd_d2H2O,
                                       source$sd_d2Hperi, H_data$d2HT_sd[1:11]))
        colnames(sources) <- c('lake_name','phyto_N_u','peri_N_u','terr_N_u','phyto_C_u','peri_C_u','terr_C_u',
                               'phyto_H_u','peri_H_u','terr_H_u','phyto_N_sd','peri_N_sd','terr_N_sd',
                               'phyto_C_sd','peri_C_sd','terr_C_sd','phyto_H_sd','peri_H_sd','terr_H_sd')

        merged <- merge.with.order(consumers, sources, all=TRUE, by='lake_name', keep_order=TRUE)
        sorted <- merged[,c(1:3,9:11,4,12:14,5,15:17,6,18:20,7,21:23,8,24:26)]
        rownames(sorted) <- 1:31
        sorted[,-(1:2)] <- apply(sorted[,-(1:2)], 2, as.numeric)

        write.csv(sorted, row.names=FALSE,
                  file=paste0('C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/mixing_models2/sensitivity/cor_',
                              i, j, '.csv'))
    }
}

for(i in c('low','ave','hi')){
    for(j in c('DC','DN', 'w')){

        source <- read.csv('cor_sources.csv')
        cons <- read.csv(paste0('sensitivity/cons_', i, j, '.csv'))

        consumers <- as.data.frame(cbind(as.character(H_data$lake_name), as.character(H_data$sample_type),
                                         cons$N_true_ind, cons$C_true_ind,
                                         cons$H_true_ind, cons$sd_N, cons$sd_C, cons$sd_H))
        colnames(consumers) <- c('lake_name','sample_type','cons_N_u','cons_C_u','cons_H_u','cons_N_sd','cons_C_sd','cons_H_sd')

        sources <- as.data.frame(cbind(as.character(H_data$Lake[1:11]), source$mean_d15NA_lake, source$mean_d15Nperi_lake,
                                       N_data$terr_d15N[1:11], source$mean_d13CA_lake, source$mean_d13Cperi_lake,
                                       C_data$terr_d13C[1:11], source$mean_d2HA_lake, source$mean_d2Hperi_lake,
                                       H_data$dD_terr_mean[1:11], source$sd_d15Na, source$sd_d15Nperi, N_data$d15NT_sd[1:11],
                                       source$sd_d13Ca, source$sd_d13Cperi, C_data$d13CT_sd[1:11], source$sd_d2H2O,
                                       source$sd_d2Hperi, H_data$d2HT_sd[1:11]))
        colnames(sources) <- c('lake_name','phyto_N_u','peri_N_u','terr_N_u','phyto_C_u','peri_C_u','terr_C_u',
                               'phyto_H_u','peri_H_u','terr_H_u','phyto_N_sd','peri_N_sd','terr_N_sd',
                               'phyto_C_sd','peri_C_sd','terr_C_sd','phyto_H_sd','peri_H_sd','terr_H_sd')

        merged <- merge.with.order(consumers, sources, all=TRUE, by='lake_name', keep_order=TRUE)
        sorted <- merged[,c(1:3,9:11,4,12:14,5,15:17,6,18:20,7,21:23,8,24:26)]
        rownames(sorted) <- 1:31
        sorted[,-(1:2)] <- apply(sorted[,-(1:2)], 2, as.numeric)

        write.csv(sorted, row.names=FALSE,
                  file=paste0('C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/mixing_models2/sensitivity/cor_',
                              i, j, '.csv'))
    }
}

# mixing ####
for(i in c('low','ave','hi')){
    for(j in c('eH','CNa','w','DC','DN')){
        #setup
        SI <- read.csv(paste0('sensitivity/cor_', i, j, '.csv'), header=TRUE)
        SI[,1:2] <- apply(SI[,1:2], MARGIN=2, FUN=as.character)

        isos <- c('C', 'N', 'H')
        which_prey <- c('phyto', 'peri', 'terr')
        npreytypes <- 3
        npreds <- nrow(SI)
        niso <- length(isos)
        alpha1 <- rep(1,npreytypes)

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

        #run model
        jags_data = list('niso'=niso, 'npreytypes'=npreytypes, 'npreds'=npreds,
                         'prey_mean_data1'=prey_mean_data, 'prey_mean_data2'=prey_mean_data,
                         'prey_var'=prey_var_data, 'preds_mean_data'=preds_mean_data,
                         'alpha'=alpha1)
        model = 'mod_mix.txt'
        jags_params = c('p', 'preds_mean', 'preds_mean_data','prey_mean_data2', 'prey_mean')

        mod <- jags.model(model, data=jags_data, inits=inits, n.chains=3, n.adapt=3000)
        update(mod, n.iter=5e3)
        j1 <- coda.samples(model=mod, variable.names=jags_params, n.iter=1e5, thin=100)

        #extract
        save(j1, file=paste0('sensitivity/allochth_', i, j, '.rda'))
        out <- reassembler_MCMClist(j1, 'p')
        MLEs <- out[[1]]
        dens_chains <- out[[2]]

        byspeclake_out <- byspeclake()
        save(byspeclake_out, file=paste0('sensitivity/byspeclake_', i, j, '.rda'))
        bylake_out <- bylake()
        save(bylake_out, file=paste0('sensitivity/bylake_', i, j, '.rda'))
    }
}

# plot ####

#just run setup and this. fetcher will err. no worries.
setwd('C:\\Users\\Mike\\git\\seven_lakes\\Analysis\\mixing_models\\mixing_models2\\sensitivity')
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

mat <- 1:4
counter <- 1
pars <- c('DC', 'DN', 'w', 'eH', 'CNa')
parconds <- c('low', 'ave', 'hi')
ps <- c('phyto', 'peri', 'terr')
lakes <- c('sm','lg','hi')

for(par in 1:5){
    for(parcond in 1:3){
        load(file=paste0('bylake_', parconds[parcond], pars[par], '.rda'))

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


lakenames <- c('Montane-Sm', 'Montane-Lg', 'Alpine')
# consnames <- c('Calanoida', 'Trichoptera', 'Cladocera')
parnames <- list(expression(bold(paste(Delta^13,'C',sep=''))),
                 expression(bold(paste(Delta^15,'N',sep=''))),
                 expression(bold(paste(omega,sep=''))),
                 expression(bold(paste('C:N'['a']))),
                 expression(bold(paste(epsilon,'H',sep=''))))

# tiff(type='cairo', res=300, width=7, height=6, units='in',
#      filename='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/official_figures/figure_4.tif')
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
mtext('Proportion t-OM', side=2, line=4.5, font=2, outer=TRUE, cex=1.15)
legend(x=3, y=1.58, xpd=NA, bty='n', lty=2:3, col=c('red', 'blue'), legend=c('Mean', 'Mode'), lwd=2, cex=1.3)
# dev.off()

# summary stats ####
w <- mat[grep('terr_\\w+_w',rownames(mat)),]
apply(w[1:3,], 2, mean) #low w
apply(w[7:9,], 2, mean) #hi w
cna <- mat[grep('terr_\\w+_CNa',rownames(mat)),] #mislabeled somehow. CNa = eH and vice-versa
apply(cna[1:3,], 2, mean)
apply(cna[7:9,], 2, mean)

cna <- mat[grep('terr_\\w+_eH',rownames(mat)),]
apply(cna[1:3,], 2, mean)
apply(cna[7:9,], 2, mean)
