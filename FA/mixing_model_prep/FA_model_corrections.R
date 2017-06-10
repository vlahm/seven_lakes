
#setup####
setwd("C:/Users/Mike/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/mixing_model_prep")
# setwd("~/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/mixing_model_prep")

rm(list=ls())
# library(graph)
# library(Rgraphviz)
# library(fastinR)
library(knitr)
library(stringr)
library(rjags)
# library(rgr)

gmean <- function (x) {
    if (is.null(dim(x))) {
        exp(mean(log(x)))
    } else {
        t(apply(x, 1, function(y) {
            exp(mean(log(x)))
        }))
    }
}

#FA conversion coefficients####

#take two
nsamp <- 3
data <- read.csv('fa.csv', header=T, row.names=1)
prey_data <- matrix(nrow=nrow(data), ncol=nsamp*10)
pred_data <- matrix(nrow=nrow(data), ncol=nsamp*10)

index <- rep(1:3, 10)
for(j in 1:ncol(prey_data)){
    for(i in 1:nrow(prey_data)){
        prey_data[i,j] <- abs(rnorm(1, data[i,index[j]], data[i,(index[j]+6)]))
        pred_data[i,j] <- abs(rnorm(1, data[i,(index[j]+3)], data[i,(index[j]+9)]))
    }
}

index2 <- seq(1,30,3)
prey_data <- t(Reduce(cbind, list(prey_data[,index2],
                                 prey_data[,(index2+1)], prey_data[,(index2+2)])))
pred_data <- t(Reduce(cbind, list(pred_data[,index2],
                                 pred_data[,(index2+1)], pred_data[,(index2+2)])))
rownames(prey_data) <- rep(c('peri_prey', 'phyto_prey', 'terr_prey'), each=10)
rownames(pred_data) <- rep(c('peri_pred', 'phyto_pred', 'terr_pred'), each=10)
colnames(prey_data) <- rownames(data)
colnames(pred_data) <- rownames(data)

typevec <- rep(1:3, each=10)
ntypes <- 3
nsamp <- 30
nfas <- nrow(data)

#data load and setup
# pred_data <- t(read.csv("allochthony_FA_consumer.csv", header = T, row.names = 1))
# prey_data <- t(read.csv("allochthony_FA_producer.csv", header = T,
#                    stringsAsFactors = F, row.names = 1))
#
# prey_types <- factor(str_match(rownames(prey_data), '(\\w*)')[,2])
# pred_types <- factor(str_match(rownames(pred_data), '\\_(\\w+)')[,2])
# pred_lakes <- factor(str_match(rownames(pred_data), '(\\w+)\\_')[,2])
# prey_lakes <- rep(pred_lakes, each=3)
# FAs <- colnames(pred_data)
# nprey <- length(levels(prey_types))
# nfas <- ncol(pred_data)
# npred <- nrow(pred_data)

# replace NAs in compositions with 0.00001
for (i in 1:nfas) {
    prey_data[is.na(prey_data[,i]) == TRUE, i] <- 0.00001
}
for (i in 1:nfas) {
    pred_data[is.na(pred_data[,i]) == TRUE, i] <- 0.00001
}

# combine FA profiles by prey type and apply closure (sum-to-1) operation

agg_profiles <- aggregate(prey_data, by=list(typevec), gmean)
agg_samp_types <- as.character(agg_profiles[,1])
agg_profiles <- agg_profiles[,-1]

closed_profiles <- matrix(nrow=nrow(agg_profiles), ncol=ncol(agg_profiles))
for(i in 1:nrow(agg_profiles)){
    for(j in 1:ncol(agg_profiles)){
        closed_profiles[i,j] <- (1/sum(agg_profiles[i,1:ncol(agg_profiles)])) * agg_profiles[i,j]
    }
}

# for Wishart hyperprior on conv_coef
# R <- diag(0.01, nfas)

#assemble and run model
data <- list('nsamp'=nsamp, 'nfas'=nfas, 'typevec'=typevec,
             'prey_gmean'=closed_profiles, 'prey_data'=prey_data, 'pred_data'=pred_data)

mod <- jags.model("FA_conversion_coeffs.txt", n.chains = 3, data=data)
update(mod, 10000)
out <- coda.samples(model=mod, variable.names=c('conv_coef'),#, 'prey_est', 'pred_est'),
                         n.iter = 1e5, thin = 100)

######continue editing from this point
############################

plot(out)
crosscorr.plot(out)

# display posterior summary
summary(samps.FA)

# get estimated discrimination from all chains:
r.samps <- do.call("rbind", samps.FA)
dim(r.samps)
fish.cc.samples <- r.samps[, seq(1, 2 * n.fats, 2)]

shrimp.cc.samples <- r.samps[, seq(2, 2 * n.fats, 2)]
fish.cc <- colMeans(fish.cc.samples)
shrimp.cc <- colMeans(shrimp.cc.samples)
# combine and write to file, repeat for 3 fish
# species in final analysis
ccs <- rbind(fish.cc, fish.cc, fish.cc, shrimp.cc)
rownames(ccs) <- unique(prey.ix)[-4] #ix should be prey_type.  he had separate vars for numeric and character
colnames(ccs) <- colnames(pred.table.FA)
write.csv(ccs, file = "cc_FA.csv")
# using independent ccs seems warranted here given
# cross-corr plot
fish.cc.var <- apply(fish.cc.samples, 2, var)
shrimp.cc.var <- apply(shrimp.cc.samples, 2, var)
# combine
ccs.var <- rbind(fish.cc.var, fish.cc.var, fish.cc.var,
                 shrimp.cc.var)
rownames(ccs.var) <- unique(prey.ix)[-4]
colnames(ccs.var) <- colnames(pred.table.FA)
write.csv(ccs.var, file = "cc_FA_var.csv")
