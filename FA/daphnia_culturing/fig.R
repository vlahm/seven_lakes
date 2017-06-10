mort = read.csv("C:/Users/Mike/git/seven_lakes/FA/daphnia_culturing/mortality_data.csv", skip=7)[,4:6]

# pdf(width=7, height=6, file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/fig6.pdf',
#     compress=FALSE)
tiff(type='cairo', res=300, width=6, height=6, units='in',
     filename='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/official_figures/figure_6.tif')
barplot(t(as.matrix(mort[2:3])), beside=TRUE, col=c('white','gray50'), width=.5,
        names.arg=c('C. ozolinii','Mont.','Mont. pr.','Alpine','Alpine pr.'),
        space=c(0,.5), ylim=c(0,80), ylab='', xlab='', yaxt='n')
mtext('Treatment', 1, font=2, line=2.5)
mtext('Number of individuals', 2, font=2, line=2.5)
axis(2, seq(0,80,10), las=2)
# lines(c(0,6),c(0,0))
legend(x=0.2, y=80, legend=c('Total mortality','Total harvest'), bty='n', fill=c('white','gray50'),
       x.intersp=.5)
dev.off()
# shell('D:/Dropbox/Grad/Projects/Thesis/Seven^ Lakes^ Project^ 2014/Manuscript/figures/fig6.pdf')
