#Michael Vlah (vlahm13@gmail.com)
#University of Washington
#data formatter for Seven Lakes allochthony models
#last edit: 7/15/16

#change input csv names and output csv name

rm(list=ls())

# setup ####
setwd("/home/mike/git/seven_lakes/final/SI/mixing_models")
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

    # tmp <- function(x) x==1; 1    # why we must check what to do if it is missing or not...
    # tmp()

    if(!missing(keep_order))
    {
        if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
        if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
        # if you didn't get "return" by now - issue a warning.
        warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
    } else {return(merge(x=x,y=y,..., sort = sort))}
}
source <- read.csv('cor_sources.csv')
cons <- read.csv('cor_consumers.csv')
raw_H <- read.csv("raw_H.csv")
raw_C <- read.csv("raw_C.csv")
raw_N <- read.csv("raw_N.csv")

# arrange ####
#consumers
consumers <- as.data.frame(cbind(as.character(raw_H$lake_name), as.character(raw_H$sample_type),
                                 cons$N_true_ind, cons$C_true_ind,
                                 cons$H_true_ind, cons$sd_N, cons$sd_C, cons$sd_H))
colnames(consumers) <- c('lake_name','sample_type','cons_N_u','cons_C_u','cons_H_u','cons_N_sd','cons_C_sd','cons_H_sd')

#sources

sources <- as.data.frame(cbind(as.character(raw_H$Lake[1:11]), source$mean_d15NA_lake, source$mean_d15Nperi_lake,
                               raw_N$terr_d15N[1:11], source$mean_d13CA_lake, source$mean_d13Cperi_lake,
                               raw_C$terr_d13C[1:11], source$mean_d2HA_lake, source$mean_d2Hperi_lake,
                               raw_H$dD_terr_mean[1:11], source$sd_d15Na, source$sd_d15Nperi, raw_N$d15NT_sd[1:11],
                               source$sd_d13Ca, source$sd_d13Cperi, raw_C$d13CT_sd[1:11], source$sd_d2H2O,
                               source$sd_d2Hperi, raw_H$d2HT_sd[1:11]))
colnames(sources) <- c('lake_name','phyto_N_u','peri_N_u','terr_N_u','phyto_C_u','peri_C_u','terr_C_u',
                       'phyto_H_u','peri_H_u','terr_H_u','phyto_N_sd','peri_N_sd','terr_N_sd',
                       'phyto_C_sd','peri_C_sd','terr_C_sd','phyto_H_sd','peri_H_sd','terr_H_sd')

# merge and order ####
merged <- merge.with.order(consumers, sources, all=TRUE, by='lake_name', keep_order=TRUE)
sorted <- merged[,c(1:3,9:11,4,12:14,5,15:17,6,18:20,7,21:23,8,24:26)]
rownames(sorted) <- 1:31
sorted[,-(1:2)] <- apply(sorted[,-(1:2)], 2, as.numeric)

# get mean and SD by lake type and consumer type ####
# ltypes <- factor(c(2,1,2,3,3,1,2,3,3,2,3,2,1,2,1,2,2,3,3,2,2,1,3,3,1,1,1,3,3,2,1))
# levels(ltypes) <- c('lg','sm','hi')
# means <- aggregate(sorted[,3:14], by=list(ltypes, sorted$sample_type), FUN=mean)
# sds <- aggregate(sorted[,3:14], by=list(ltypes, sorted$sample_type), FUN=sd)
# sds2 <- aggregate(sorted[,c(18,22,26)], by=list(ltypes, sorted$sample_type), FUN=mean)
#
# means <- means[-c(4,5,8),]
# sds <- sds[-c(4,5,8),]
# sds2 <- sds2[-c(4,5,8),]
#
# rownames(means) <- 1:8
# inds <- c(4:6,8:10,12:14)
# inds2 <- c(16:18,20:22,24:26)
# temp <- matrix(NA,nrow=3,ncol=length(inds))
# temp2 <- matrix(NA,nrow=3,ncol=length(inds))
# names <- vector(length=length(inds))
# names2 <- vector(length=length(inds))
# for(i in 1:length(inds)){
#     temp[1,i] <- mean(sorted[c(1,3,7,10),inds[i]])
#     temp[2,i] <- mean(sorted[c(2,6,27),inds[i]])
#     temp[3,i] <- mean(sorted[c(4,5,8,9),inds[i]])
#     temp2[1,i] <- sd(sorted[c(1,3,7,10),inds2[i]])
#     temp2[2,i] <- sd(sorted[c(2,6,27),inds2[i]])
#     temp2[3,i] <- sd(sorted[c(4,5,8,9),inds2[i]])
#     names[i] <- colnames(sorted)[inds[i]]
#     names2[i] <- colnames(sorted)[inds2[i]]
# }
# colnames(temp) <- names
# rownames(temp) <- c('sm','lg','hi')
# colnames(temp2) <- names2
# rownames(temp2) <- c('sm','lg','hi')
# means[1,c(4:6,8:10,12:14)] <- temp[2,1:9]
# means[4,c(4:6,8:10,12:14)] <- temp[2,1:9]
# means[6,c(4:6,8:10,12:14)] <- temp[2,1:9]
# means[2,c(4:6,8:10,12:14)] <- temp[1,1:9]
# means[5,c(4:6,8:10,12:14)] <- temp[1,1:9]
# means[7,c(4:6,8:10,12:14)] <- temp[1,1:9]
# means[3,c(4:6,8:10,12:14)] <- temp[3,1:9]
# means[8,c(4:6,8:10,12:14)] <- temp[3,1:9]
# sds[1,c(4:6,8:10,12:14)] <- temp2[2,1:9]
# sds[4,c(4:6,8:10,12:14)] <- temp2[2,1:9]
# sds[6,c(4:6,8:10,12:14)] <- temp2[2,1:9]
# sds[2,c(4:6,8:10,12:14)] <- temp2[1,1:9]
# sds[5,c(4:6,8:10,12:14)] <- temp2[1,1:9]
# sds[7,c(4:6,8:10,12:14)] <- temp2[1,1:9]
# sds[3,c(4:6,8:10,12:14)] <- temp2[3,1:9]
# sds[8,c(4:6,8:10,12:14)] <- temp2[3,1:9]
#
# agg <- cbind(means, sds[3:ncol(sds)])
# agg[,18] <- sds2[,3]; agg[,22] <- sds2[,4]; agg[,26] <- sds2[,5]
# colnames(agg)[1:2] <- c('lake_type', 'sample_type')
# colnames(agg)[15:ncol(agg)] <- colnames(sorted)[15:ncol(agg)]

# write output files ####
# write.csv(agg, row.names=FALSE,
#           file="C:/Users/Mike/git/seven_lakes/Analysis/mixing_models/mixing_models2/cor_agg.csv")
write.csv(sorted, row.names=FALSE, file='cor_all.csv')
