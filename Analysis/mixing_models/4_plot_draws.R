dim(j1$samples[[1]])
dimnames(j1$samples[[1]])

#start witht he part after the divider

phytox <- j1$samples[[1]][,224][1]
phytoy <- j1$samples[[1]][,227][1]
perix <- j1$samples[[1]][,225][1]
periy <- j1$samples[[1]][,228][1]
terrx <- j1$samples[[1]][,226][1]
terry <- j1$samples[[1]][,229][1]
consx <- j1$samples[[1]][,157][1]
consy <- j1$samples[[1]][,188][1]
xmin <- min(phytox, perix, terrx, consx, j1$samples[[1]][,95][1:1000])
xmax <- max(phytox, perix, terrx, consx, j1$samples[[1]][,95][1:1000])
ymin <- min(phytoy, periy, terry, consy, j1$samples[[1]][,125][1:1000])
ymax <- max(phytoy, periy, terry, consy, j1$samples[[1]][,125][1:1000])

plot(y=j1$samples[[1]][,95][1:1000], x=j1$samples[[1]][,125][1:1000],
     pch=1, xlim=c(ymin, ymax), ylim=c(xmin,xmax), xlab='d2H', ylab='d15N')
points(y=j1$samples[[2]][,95][1:1000], x=j1$samples[[2]][,125][1:1000])
points(y=j1$samples[[3]][,95][1:1000], x=j1$samples[[3]][,125][1:1000])
points(y=phytox, x=phytoy,
       pch=20, cex=3, col='steelblue1')
points(y=perix, x=periy,
       pch=20, cex=3, col='darkgreen')
points(y=terrx, x=terry,
       pch=20, cex=3, col='sienna4')
points(y=consx, x=consy,
       pch=20, cex=3, col='yellow')


# ####
# num=12
thingy <- function(num){
phytox <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data[1,1,',num,']'))][1]
phytoy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data[1,2,',num,']'))][1]
perix <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data[2,1,',num,']'))][1]
periy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data[2,2,',num,']'))][1]
terrx <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data[3,1,',num,']'))][1]
terry <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data[3,2,',num,']'))][1]
consx <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',1]'))][1]
consy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',2]'))][1]
drawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000],
           j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000],
           j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000])
drawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000],
           j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000],
           j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000])
xmin <- min(phytox, perix, terrx, consx, drawx)
xmax <- max(phytox, perix, terrx, consx, drawx)
ymin <- min(phytoy, periy, terry, consy, drawy)
ymax <- max(phytoy, periy, terry, consy, drawy)

plot(y=drawx, x=drawy,
     pch=1, xlim=c(ymin, ymax), ylim=c(xmin,xmax), xlab='d2H', ylab='d15N', main=paste(SI[num,1],SI[num,2]))
points(y=phytox, x=phytoy,
       pch=20, cex=3, col='steelblue1')
points(y=perix, x=periy,
       pch=20, cex=3, col='darkgreen')
points(y=terrx, x=terry,
       pch=20, cex=3, col='sienna4')
points(y=consx, x=consy,
       pch=20, cex=3, col='yellow')
}
pdf(file="C:/Users/Mike/Desktop/drawsNH_scaled.pdf", onefile=T)
for(i in 1:length(predrows)){
    thingy(i)
}
dev.off()

#for fully bayesian + resid####
thingy2 <- function(num){
    phytox <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[1,1,',num,']'))][1]
    phytoy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[1,2,',num,']'))][1]
    perix <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[2,1,',num,']'))][1]
    periy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[2,2,',num,']'))][1]
    terrx <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[3,1,',num,']'))][1]
    terry <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[3,2,',num,']'))][1]
    consx <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',1]'))][1]
    consy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',2]'))][1]
    drawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000],
               j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000],
               j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000])
    drawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000],
               j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000],
               j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000])
    phytodrawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,1,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,1,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,1,',num,']'))][1:1000])
    peridrawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,1,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,1,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,1,',num,']'))][1:1000])
    terrdrawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,1,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,1,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,1,',num,']'))][1:1000])
    phytodrawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000])
    peridrawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000])
    terrdrawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000])
    xmin <- min(phytox, perix, terrx, consx, drawx, peridrawx, phytodrawx, terrdrawx, na.rm=T)
    xmax <- max(phytox, perix, terrx, consx, drawx, peridrawx, phytodrawx, terrdrawx, na.rm=T)
    ymin <- min(phytoy, periy, terry, consy, drawy, peridrawy, phytodrawy, terrdrawy, na.rm=T)
    ymax <- max(phytoy, periy, terry, consy, drawy, peridrawy, phytodrawy, terrdrawy, na.rm=T)

    plot(y=drawx, x=drawy,
         pch=1, xlim=c(ymin, ymax), ylim=c(xmin,xmax), xlab='d2H', ylab='d15N', main=paste(SI[num,1],SI[num,2]))
    points(y=phytodrawx, x=phytodrawy, col='steelblue1')
    points(y=peridrawx, x=peridrawy, col='darkgreen')
    points(y=terrdrawx, x=terrdrawy, col='sienna4')
    points(y=phytox, x=phytoy,
           pch=21, cex=3, bg='steelblue1', col='white')
    points(y=perix, x=periy,
           pch=21, cex=3, bg='darkgreen', col='white')
    points(y=terrx, x=terry,
           pch=21, cex=3, bg='sienna4', col='white')
    points(y=consx, x=consy,
           pch=21, cex=3, bg='yellow', col='white')
}
pdf(file="C:/Users/Mike/Desktop/allindiv_recon.pdf", onefile=T)
for(i in 1:npreds){
    thingy2(i)
}
dev.off()

#for fully bayesian + resid (all three isotopes) has issue ####
thingy3 <- function(num){
    phytox <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[1,1,',num,']'))][1]
    phytoy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[1,2,',num,']'))][1]
    phytoz <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[1,3,',num,']'))][1]
    perix <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[2,1,',num,']'))][1]
    periy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[2,2,',num,']'))][1]
    periz <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[2,3,',num,']'))][1]
    terrx <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[3,1,',num,']'))][1]
    terry <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[3,2,',num,']'))][1]
    terrz <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean_data2[3,3,',num,']'))][1]
    consx <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',1]'))][1]
    consy <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',2]'))][1]
    consz <- j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean_data[',num,',3]'))][1]
    drawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000],
               j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000],
               j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',1]'))][1:1000])
    drawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000],
               j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000],
               j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',2]'))][1:1000])
    drawz <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',3]'))][1:1000],
               j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',3]'))][1:1000],
               j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('preds_mean[',num,',3]'))][1:1000])
    phytodrawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,1,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,1,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,1,',num,']'))][1:1000])
    peridrawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,1,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,1,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,1,',num,']'))][1:1000])
    terrdrawx <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,1,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,1,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,1,',num,']'))][1:1000])
    phytodrawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000])
    peridrawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000])
    terrdrawy <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000])
    phytodrawz <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,3,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[1,2,',num,']'))][1:1000])
    peridrawz <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,3,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[2,2,',num,']'))][1:1000])
    terrdrawz <- c(j1$samples[[1]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,3,',num,']'))][1:1000])
               # j1$samples[[2]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000],
               # j1$samples[[3]][,which(dimnames(j1$samples[[1]])[[2]]==paste0('prey_mean[3,2,',num,']'))][1:1000])
    xmin <- min(phytox, perix, terrx, consx, drawx, peridrawx, phytodrawx, terrdrawx, na.rm=T)
    xmax <- max(phytox, perix, terrx, consx, drawx, peridrawx, phytodrawx, terrdrawx, na.rm=T)
    ymin <- min(phytoy, periy, terry, consy, drawy, peridrawy, phytodrawy, terrdrawy, na.rm=T)
    ymax <- max(phytoy, periy, terry, consy, drawy, peridrawy, phytodrawy, terrdrawy, na.rm=T)
    zmin <- min(phytoz, periz, terrz, consz, drawz, peridrawz, phytodrawz, terrdrawz, na.rm=T)
    zmax <- max(phytoz, periz, terrz, consz, drawz, peridrawz, phytodrawz, terrdrawz, na.rm=T)

    plot(y=drawx, x=drawz,
         pch=1, xlim=c(zmin, zmax), ylim=c(xmin,xmax), xlab='d2H', ylab='d13C', main=paste(count, SI[num,1],SI[num,2],'HxC'))
    points(y=phytodrawx, x=phytodrawz, col='steelblue1')
    points(y=peridrawx, x=peridrawz, col='darkgreen')
    points(y=terrdrawx, x=terrdrawz, col='sienna4')
    points(y=phytox, x=phytoz,
           pch=21, cex=3, bg='steelblue1', col='white')
    points(y=perix, x=periz,
           pch=21, cex=3, bg='darkgreen', col='white')
    points(y=terrx, x=terrz,
           pch=21, cex=3, bg='sienna4', col='white')
    points(y=consx, x=consz,
           pch=21, cex=3, bg='yellow', col='white')

    plot(y=drawy, x=drawz,
         pch=1, xlim=c(zmin, zmax), ylim=c(ymin,ymax), xlab='d2H', ylab='d15N', main=paste(count, SI[num,1],SI[num,2],'HxN'))
    points(y=phytodrawy, x=phytodrawz, col='steelblue1')
    points(y=peridrawy, x=peridrawz, col='darkgreen')
    points(y=terrdrawy, x=terrdrawz, col='sienna4')
    points(y=phytoy, x=phytoz,
           pch=21, cex=3, bg='steelblue1', col='white')
    points(y=periy, x=periz,
           pch=21, cex=3, bg='darkgreen', col='white')
    points(y=terry, x=terrz,
           pch=21, cex=3, bg='sienna4', col='white')
    points(y=consy, x=consz,
           pch=21, cex=3, bg='yellow', col='white')

    plot(x=drawx, y=drawy,
         pch=1, xlim=c(xmin, xmax), ylim=c(ymin,ymax), xlab='d13C', ylab='d15N', main=paste(count, SI[num,1],SI[num,2],'CxN'))
    points(x=phytodrawx, y=phytodrawy, col='steelblue1')
    points(x=peridrawx, y=peridrawy, col='darkgreen')
    points(x=terrdrawx, y=terrdrawy, col='sienna4')
    points(x=phytox, y=phytoy,
           pch=21, cex=3, bg='steelblue1', col='white')
    points(x=perix, y=periy,
           pch=21, cex=3, bg='darkgreen', col='white')
    points(x=terrx, y=terry,
           pch=21, cex=3, bg='sienna4', col='white')
    points(x=consx, y=consy,
           pch=21, cex=3, bg='yellow', col='white')
}
pdf(file="C:/Users/Mike/Desktop/fixterr.pdf", onefile=T)
count = 1
for(i in 1:npreds){
    thingy3(i)
    count <- count + 1
}
dev.off()
