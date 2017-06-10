#Michael Vlah
#University of Washington
#vlahm13@gmail.com
#last edit: 2/26/16

rm(list=ls(all=TRUE))

#setup####
nlakes <- 11
nsamp <- 31
group_vector_1 <- c(1,1,2,2,2,2,3,3,3,1,3)
group_vector_2 <- c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1)
species_vector <- c(rep(1,10),2,rep(3,6),4,4,rep(5,11),2)
nspecies <- length(levels(factor(species_vector)))
ngroups <- length(levels(factor(group_vector_1)))

# if (is.null(dev.list()) == FALSE){dev.off()}
# windows(record=T) #open plot window (on Windows only) --use quartz() or x11() for Mac
library(rjags)
library(RColorBrewer)
library(car)
setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/H_Isotopes/")
#setwd("~/Desktop/grad/Projects/Thesis/Seven Lakes Project 2014/Data/H_Isotopes")
H_data <- read.csv("H_data.csv")
C_data <- read.csv("C_data.csv")
N_data <- read.csv("N_data.csv")
#covariates <- read.csv("C:/Users/Mike/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Analysis/7lmeans.csv")
#iso <- read.csv("C:/Users/Mike/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Analysis/Isomeans.csv")
#CN_uncertainty <- read.csv("C:/Users/Mike/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Analysis/CN_uncertainty.csv")

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

max_dens <- function(mcmc.list){

    #extract all relevant info from mcmc.list
    nchains <- length(mcmc.list)
    chainlength <- length(mcmc.list[[1]][,1])
    paramlist <- gsub(" *\\[.*?\\] *", '', colnames(mcmc.list[[1]]))
    params <- levels(factor(paramlist))
    nparams <- length(params)
    param_lengths <- vector('numeric', length=nparams)
    for(i in 1:nparams){
        param_lengths[i] <- length(which(paramlist==paste(params[i])))
    }
    longest_param <- max(param_lengths)

    #create matrix of "MLE" values across all chains
    temp <- matrix(data=NA, nrow=chainlength, ncol=nchains)
    out <- matrix(data=NA, nrow=longest_param, ncol=nparams)
    for(j in 1:nparams){
        indices <- which(paramlist==params[j])
        for(i in 1:length(indices)){
            for(k in 1:nchains){
                temp[,k] <- mcmc.list[[k]][,indices[i]][1:length(mcmc.list[[k]][,1])]
            }
            temp2 <- as.vector(temp)
            out[i,j] <- density(temp2, n=10000)$x[which.max(density(temp2, n=10000)$y)]
        }
    }
    colnames(out) <- params

    return(out)
}

#estimate phyto dD, correct for dietary water####

#initialize the model
Hmod <- "C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\sensitivity_analysis\\hard_code_H.txt"
# H_estim_and_cor <- jags.model('H_phytoByLake_consBySpecLake.txt',
H_estim_and_cor <- jags.model(Hmod,
                            data = list('nlakes' = nlakes, 'd2H2O_pt1' = H_data[1:11,3], 'nsamp' = nsamp,
                                        'd2H2O_pt2' = H_data[,9], 'cons_meas' = H_data[,7],
                                        'trophic_lvls' = H_data[,8], 'grps1' = group_vector_1,
                                        'ngroups' = ngroups, 'nspecies' = nspecies,
                                        'specs' = species_vector, 'grps2' = group_vector_2,
                                        'mean_w_sample'=rep(.28, 31)),  #hard-code here
                            inits = function ()
                            {
                              list('mean_phyto_global' = rnorm(1, -206.78, 22.07),
                                   #'mean_eH_global' = rnorm(1, -150, 28),
                                   'xi' = rnorm(1, 0, 1/25^2), #25 is scale parameter on half-Cauchy
                                   'tau.eta' = rgamma(1, shape=0.5, scale=1/0.5),
                                   #'mean_w_global' = rnorm(1, 0.27, 0.14),
                                   'mean_samp_global' = runif(1, -500, 1),
                                   'xi.2' = rnorm(1, 0, 1/33^2),
                                   'tau.eta.2' = rgamma(1, shape=0.5, scale=1/0.5),
                                   'xi.3' = rnorm(ngroups, 0, 1/24^2),
                                   'tau.eta.3' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                   'xi.4' = rnorm(ngroups, 0, 1/33^2),
                                   'tau.eta.4' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                   'xi.6' = rnorm(1, 0, 1/26^2),
                                   'tau.eta.6' = rgamma(1, shape=0.5, scale=1/0.5),
                                   'xi.7' = rnorm(ngroups, 0, 1/.30^2),
                                   'tau.eta.7' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                   'sd_w_species' = runif(1, 0, 1),
                                   'sd_w_sample' = runif(nspecies, 0, 1),
                                   'xi.5' = rnorm(nsamp, 0, 1/30^2),
                                   'tau.eta.5' = rgamma(nsamp, shape=0.5, scale=1/0.5))
                            },
                            n.chains = 3,
                            n.adapt = 2000)

#burnin
update(H_estim_and_cor, n.iter=5e3)

#run
H_samples.out <- coda.samples(H_estim_and_cor, c('mean_phyto_lake', 'sd_phyto_lake',
                                                 'mean_eH_lake', 'd2H2O_pred', 'samp_correct',
                                                 'sd_samp', 'mean_w_sample',
                                                 'sd_w_sample'), n.iter=1e5, thin=100)

#diagnostics (more in Git branch "testing")
#head(H_samples.out.keep[[1]])
#colnames(H_samples.out.keep[[1]])[-c(65:95)]
#plot(H_samples.out.keep)
#gelman.plot(H_samples.out.keep[,1]) #shouldn't vary from 1 by more than .05

#extract values with highest likelihoods for all outputs
#verify convergence and unimodality first
H_dens <- max_dens(H_samples.out)

#write to CSV
names1 <- as.character(H_data[,1])
names2 <- as.character(H_data[,5])
write.csv(Reduce(cbind, list(names1, names2, H_dens)),
          file="C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\sensitivity_analysis\\ave_eH.csv")
          # file="C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\H_phytoByLake_consBySpecLake.csv")

#estimate phyto d13C, correct for trophic enrichment and compounding####

#initialize the model
Cmod <- "C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\sensitivity_analysis\\hard_code_C.txt"
C_estim_and_cor <- jags.model(Cmod,
                              data = list('nlakes' = nlakes, 'cn_terr' = C_data[1:11,3],
                                          'd13C_terr' = C_data[1:11,5], 'cnPOM' = C_data[1:11,2],
                                          'd13C_POM' = C_data[1:11,4], 'cons_meas' = C_data[,6],
                                          'trophic_lvls' = C_data[,7], 'grps1' = group_vector_1,
                                          'nsamp' = nsamp, 'ngroups' = ngroups,
                                          'd13C_enrich'=rep(.39, 31)), #hard-code here
                              inits = function ()
                              {
                                  list(#'mean_cnA_global' = rnorm(1, 4.74, 1.2),
                                       'cnT' = rnorm(nlakes, 49.72, 10),
                                       'mean_d13CA_global' = rnorm(1, -28.62, 2.56),
                                       'd13CT' = rnorm(nlakes, -27.12, 5),
                                       'xi' = rnorm(1, 0, 1/1.7^2), #1.7 is scale parameter on half-Cauchy
                                       'tau.eta' = rgamma(1, shape=0.5, scale=1/0.5),
                                       'xi.2' = rnorm(1, 0, 1/2^2),
                                       'tau.eta.2' = rgamma(1, shape=0.5, scale=1/0.5),
                                       'xi.3' = rnorm(ngroups, 0, 1/2^2),
                                       'tau.eta.3' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                       'xi.4' = rnorm(ngroups, 0, 1/2.8^2),
                                       'tau.eta.4' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                       'xi.5' = rnorm(nsamp, 0, 1/8^2),
                                       'tau.eta.5' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                                      # 'd13C_enrich' = rnorm(nsamp, 0.39, 1.14),
                                       'mean_cons_true' = runif(nsamp, -50, 0))
                              },
                              n.chains = 3,
                              n.adapt = 2000)

#burnin
update(C_estim_and_cor, n.iter=5e3)

#run
C_samples.out <- coda.samples(C_estim_and_cor, c('mean_cnA_lake', 'cnT', 'mean_d13CA_lake',
                                                 'd13C_POM_pred', 'd13CT',
                                                 'fracT', 'mean_cons_true', 'd13C_enrich',
                                                 'sd_cons_meas', 'sd_d13CA_lake'),
                              thin=100, n.iter=1e5)

#diagnostics (more in Git branch "testing")
# colnames(C_samples.out.keep[[1]])
# plot(C_samples.out.keep)
# gelman.plot(C_samples.out.keep) #shouldn't vary from 1 by more than .05

#extract values with highest likelihoods for all outputs
#verify convergence and unimodality first
C_dens <- max_dens(C_samples.out)

#write to CSV
names1 <- as.character(H_data[,1])
names2 <- as.character(H_data[,5])
write.csv(Reduce(cbind, list(names1, names2, C_dens)),
          file="C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\sensitivity_analysis\\tfracC_ave.csv")

#estimate phyto d15N, correct for trophic enrichment and compounding####

#initialize the model
Nmod <- "C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\sensitivity_analysis\\hard_code_N.txt"
# N_estim_and_cor <- jags.model('N_phytoByLake_indivEnrich.txt',
N_estim_and_cor <- jags.model(Nmod,
                              data = list('nlakes' = nlakes, 'cn_terr' = N_data[1:11,3],
                                          'd15N_terr' = N_data[1:11,5], 'cnPOM' = N_data[1:11,2],
                                          'd15N_POM' = N_data[1:11,4], 'cons_meas' = N_data[,6],
                                          'trophic_lvls' = N_data[,7], 'grps1' = group_vector_1,
                                          'nsamp' = nsamp, 'ngroups' = ngroups,
                                          'd15N_enrich'=rep(3.4, 31)), #hard-code here
                              inits = function ()
                              {
                                  list(#'mean_cnA_global' = rnorm(1, 4.74, 1.2),
                                       'cnT' = rnorm(nlakes, 49.72, 10),
                                       'mean_d15NA_global' = rnorm(1, 0.68, 1.23),
                                       'd15NT' = rnorm(nlakes, -27.12, 5),
                                       'xi' = rnorm(1, 0, 1/1.7^2), #1.7 is scale parameter on half-Cauchy
                                       'tau.eta' = rgamma(1, shape=0.5, scale=1/0.5),
                                       'xi.2' = rnorm(1, 0, 1/1^2),
                                       'tau.eta.2' = rgamma(1, shape=0.5, scale=1/0.5),
                                       'xi.3' = rnorm(ngroups, 0, 1/2^2),
                                       'tau.eta.3' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                       'xi.4' = rnorm(ngroups, 0, 1/1.5^2),
                                       'tau.eta.4' = rgamma(ngroups, shape=0.5, scale=1/0.5),
                                       'xi.5' = rnorm(nsamp, 0, 1/2.8^2),
                                       'tau.eta.5' = rgamma(nsamp, shape=0.5, scale=1/0.5),
                                       #'d15N_enrich' = rnorm(nsamp, 3.4, 0.99),
                                       'mean_cons_true' = runif(nsamp, -15, 10))
                              },
                              n.chains = 3,
                              n.adapt = 2000)

#burnin
update(N_estim_and_cor, n.iter=5e3)

#run
N_samples.out <- coda.samples(N_estim_and_cor, c('mean_cnA_lake', 'cnT', 'mean_d15NA_lake',
                                                 'd15N_POM_pred', 'd15NT',
                                                 'fracT', 'mean_cons_true', 'd15N_enrich',
                                                 'sd_cons_meas', 'sd_d15NA_lake'),
                              thin=100, n.iter=1e5)

#diagnostics (more in Git branch "testing")
# colnames(N_samples.out.keep[[1]])
# plot(N_samples.out.keep)
# gelman.plot(N_samples.out.keep[,1]) #shouldn't vary from 1 by more than .05

#extract values with highest likelihoods for all outputs
#verify convergence and unimodality first

N_dens <- max_dens(N_samples.out)

#write to CSV
names1 <- as.character(H_data[,1])
names2 <- as.character(H_data[,5])
write.csv(Reduce(cbind, list(names1, names2, N_dens)),
          file="C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\Seven Lakes Project 2014\\Analysis\\mixing_models\\sensitivity_analysis\\tfracN_ave.csv")

#########obsolete after this point##################

#extract posterior means and sds from jags.model output
H_s <- summary(H_samples.out.keep)
H_m <- H_s$statistics[,'Mean']  #this extracts a vector of named numbers representing the means of the posteriors
H_sd <- H_s$statistics[,'SD']  #same thing for SDs

#rearrange into a new data.frame
H_cons_cor <- unname(H_m[1:31])
H_d2H2O_pred <- unname(H_m[32:42])
H_eH <- unname(H_m[43:53])
H_mean_phyto <- unname(H_m[54:64])
H_prec_cons <- unname(H_m[65:95])
H_sd_cons <- unname(H_m[96:126])
H_sd_phyto <- unname(H_m[127:137])
H_mean_terr <- H_data[,11]
H_sd_terr <- H_data[,12]

temp1 <- as.data.frame(cbind(as.character(H_data$lake_name), as.character(H_data$sample_type), H_cons_cor,
                             H_prec_cons, H_mean_terr, H_sd_terr))
temp2 <- as.data.frame(cbind(as.character(H_data$Lake[1:11]), H_mean_phyto, H_sd_phyto))
H_data_part2 <- merge.with.order(temp1, temp2, by='V1', all=TRUE, sort=FALSE, keep_order=1)
colnames(H_data_part2)[1:2] <- c('lake', 'cons_type')
rownames(H_data_part2) <- 1:31
H_data_part2[,-(1:2)] <- apply(H_data_part2[,-(1:2)], MARGIN=2, FUN=as.numeric)

#covariates added here:
# H_data_part2 <- merge.with.order(x=H_data_part2, y=covariates, by='lake', all.x=TRUE, sort=FALSE, keep_order=1)
# rownames(H_data_part2) <- 1:31
# for (i in 10:52){
#   plot(H_data_part2$cons_mean, H_data_part2[,i], ylab=colnames(H_data_part2)[i])
# }
# mod1 <- lm(H_data_part2$cons_mean[-c(3,14)] ~ H_data_part2$color_chlA[-c(3,14)])
# summary(mod1)
# anova(mod1)
# plot(H_data_part2$cons_mean[-c(3,14)] ~ H_data_part2$color_chlA[-c(3,14)])
# abline(mod1)
#
# mod1 <- lm(H_data_part2$cons_mean ~ H_data_part2$DOC)
# summary(mod1)
# anova(mod1)
#
# mod1 <- lm(H_data_part2$cons_mean ~ H_data_part2$chlA)
# summary(mod1)
# anova(mod1)
#
# mod1 <- lm(H_data_part2$cons_mean ~ H_data_part2$dD_H2O)
# summary(mod1)
# anova(mod1)
#
# mod1 <- lm(H_data_part2$cons_mean ~ H_data_part2$d18O_H2O)
# summary(mod1)
# anova(mod1)

#estimate allochthony for zoop, caddis, fish consumers####

allochth <- jags.model('allochthony.bug',
                       data = list('mean_phyto' = H_data_part2[,7], 'sd_phyto' = H_data_part2[,8],
                                   'mean_terr' = H_data_part2[,5], 'sd_terr' = H_data_part2[,6],
                                   'cons_cor' = H_data_part2[,3], 'prec_cons' = H_data_part2[,4],
                                   'nsamp' = nsamp),
                       inits = function ()
                       {
                         list('d2H_phyto' = rnorm(nsamp, mean(H_mean_phyto), 0.1),#mean(H_sd_phyto)),
                              'd2H_terr' = rnorm(nsamp, mean(H_mean_terr), 0.1),#mean(H_sd_terr)),
                              'phi_phyto' = runif(nsamp, -5, 5))
                       },
                       n.chains = 4,
                       n.adapt = 2000)

#run, thin, burn, etc.
samples.out.2 <- coda.samples(allochth, c('phi_phyto', 'd2H_phyto',
                                          'd2H_terr', 'cons_pred'), n.iter=1e5)
samples.out.keep.2 <- window(samples.out.2, start=5e3, end=1e5, thin=100)

#plot diagnostics
plot(samples.out.keep.2)

#extract posterior means and sds from jags.model output
s.2 <- summary(samples.out.keep.2)
m.2 <- s.2$statistics[,'Mean']  #this extracts a vector of named numbers representing the means of the posteriors
sd.2 <- s.2$statistics[,'SD']  #same thing for SDs

#rearrange into a final data.frame
cons_mean_final <- unname(m.2[1:31])
phyto_mean_final <- unname(m.2[32:62])
terr_mean_final <- unname(m.2[63:93])
phi_phyto_mean <- unname(m.2[94:124])
phi_terr_mean <- unname(m.2[125:155])

cons_sd_final <- unname(sd.2[1:31])
phyto_sd_final <- unname(sd.2[32:62])
terr_sd_final <- unname(sd.2[63:93])
phi_phyto_sd <- unname(sd.2[94:124])
phi_terr_sd <- unname(sd.2[125:155])

output.allochth <- as.data.frame(cbind(as.character(H_data_part2$lake), cons_mean_final, phyto_mean_final, terr_mean_final, phi_phyto_mean,
                                       phi_terr_mean, cons_sd_final, phyto_sd_final, terr_sd_final, phi_phyto_sd, phi_terr_sd))
colnames(output.allochth)[1] <- 'lake'
output.allochth[,-1] <- apply(output.allochth[,-1], MARGIN=2, FUN=as.numeric)

#plot output, return table####

#add other isotopes to frame
temp <- CN_uncertainty[c(1:17,54:66,75),c(2,5,8,10,13)]
output.allochth <- cbind(output.allochth, temp)
rownames(output.allochth) <- 1:31

write.csv(x = output.allochth, row.names=FALSE,
          file = "C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/H_Isotopes/allochth_out.csv")

#plot
# cols <- brewer.pal(11, 'Paired')
# col_temp <- cbind(levels(output.allochth$lake), cols)
# col_temp <- merge.with.order(output.allochth[,1:2], col_temp, by.x='lake', by.y='V1', all.x=TRUE, sort=FALSE, keep_order=1)
# col_vec <- as.character(col_temp[,3])

cols <- brewer.pal(5, 'Set1')
col_temp <- cbind(levels(factor(as.character(output.allochth$type))), cols)
col_temp2 <- merge.with.order(output.allochth[,c(1,12)], col_temp, by.x='type', by.y='V1', all.x=TRUE, sort=FALSE, keep_order=1)
col_vec <- as.character(col_temp2[,3])

par(mfrow=c(4,3))
for (i in levels(output.allochth$lake)){
  plot(output.allochth$d15N_cor[output.allochth$lake==i],
       output.allochth$cons_mean_final[output.allochth$lake==i],
       col='black', pch=21, cex=3, bg=col_vec[output.allochth$lake==i], main=paste(i),
       xlab='d15N', ylab='d2H', ylim=c(min(output.allochth[,2:4]), max(output.allochth[,2:4])))
  text(output.allochth$d15N_cor[output.allochth$lake==i],
       output.allochth$cons_mean_final[output.allochth$lake==i],
       labels=rownames(output.allochth)[output.allochth$lake==i])
  points(output.allochth$d15N_cor[output.allochth$lake==i],
         output.allochth$terr_mean_final[output.allochth$lake==i],
         col='black', pch=22, cex=2, bg=col_vec[output.allochth$lake==i])
  points(output.allochth$d15N_cor[output.allochth$lake==i],
         output.allochth$phyto_mean_final[output.allochth$lake==i],
         col='black', pch=24, cex=2, bg=col_vec[output.allochth$lake==i])
  #   error.bars(x=output.allochth$d15N_cor[output.allochth$lake==i],
  #              y=output.allochth$cons_mean_final[output.allochth$lake==i],
  #              upper=output.allochth$cons_sd_final[output.allochth$lake==i],
  #              cap.length=0.05, horiz=F)
  # #   error.bars(x=output.allochth$d15N_cor[output.allochth$lake==i],
  # #              y=output.allochth$cons_mean_final[output.allochth$lake==i],
  # #              upper=output.allochth$d15N_sd_tot[output.allochth$lake==i],
  # #              cap.length=0.05, horiz=T)
  #   error.bars(x=output.allochth$d15N_cor[output.allochth$lake==i],
  #              y=output.allochth$phyto_mean_final[output.allochth$lake==i],
  #              upper=output.allochth$phyto_sd_final[output.allochth$lake==i],
  #              cap.length=0.05, horiz=F)
  # #   error.bars(x=output.allochth$d15N_cor[output.allochth$lake==i],
  # #              y=output.allochth$phyto_mean_final[output.allochth$lake==i],
  # #              upper=output.allochth$d15N_sd_tot[output.allochth$lake==i],
  # #              cap.length=0.05, horiz=T)
  #   error.bars(x=output.allochth$d15N_cor[output.allochth$lake==i],
  #              y=output.allochth$terr_mean_final[output.allochth$lake==i],
  #              upper=output.allochth$terr_sd_final[output.allochth$lake==i],
  #              cap.length=0.05, horiz=F)
  # #   error.bars(x=output.allochth$d15N_cor[output.allochth$lake==i],
  # #              y=output.allochth$terr_mean_final[output.allochth$lake==i],
  # #              upper=output.allochth$d15N_sd_tot[output.allochth$lake==i],
  # #              cap.length=0.05, horiz=T)
}
plot(1,1, type='n', xaxt='n', yaxt='n', xlab='', ylab='')
legend(x='center', legend=c('terr', 'cons', 'phyto', 'cal', 'chaob', 'clad', 'fish', 'caddis'),
       pch=c(22,21,24,15,15,15,15,15), bg=c('white', 'white', 'white', NULL,NULL,NULL,NULL,NULL),
       col=c('black', 'black', 'black', col_temp[,2]))
par(mfrow=c(1,1))

# }
run()
