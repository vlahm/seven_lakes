rm(list=ls())
#first attempt with binomial data - not enough df ####
data <- read.csv("C:/Users/Mike/git/seven_lakes/FA/daphnia_culturing/for_analysis.csv")
                 row.names=1)
trt <- factor(rownames(data))

mod1 <- glm(prop.mort~trt, weights=total.ind, family=binomial, data=data)
anodev1<-anova(mod1); anodev1
fstat1<-(anodev1$Dev[2]/anodev1$Df[2])/(anodev1[2,"Resid. Dev"]/anodev1[2,"Resid. Df"]); fstat1
1-pf(fstat1,4,0)

###not enough df for this.  why do i only have 4 df?

harv <- data$Total.Harvest+10
mort <- data$Total.Mortality+10
fac <- factor(c(1,2,3,4,4))

mod1 <- glm(cbind(mort,harv)~fac, family=binomial(link="logit"))
anodev1<-anova(mod1); anodev1
fstat1<-(anodev1$Dev[2]/anodev1$Df[2])/(anodev1[2,"Resid. Dev"]/anodev1[2,"Resid. Df"]); fstat1
1-pf(fstat1,3,1)

#second attempt with bernoulli data ####
#i artificially added one survivor to each of the raw treatments just so the model would work
data_bern <- read.csv("C:/Users/Mike/git/seven_lakes/FA/daphnia_culturing/for_analysis.csv/for_analysis_bernoulli.csv")

mod1 <- glm(cbind(mort,harv)~trt, family=binomial(link="cloglog"), data=data_bern)
anodev1<-anova(mod1); anodev1
fstat1<-(anodev1$Dev[2]/anodev1$Df[2])/(anodev1[2,"Resid. Dev"]/anodev1[2,"Resid. Df"]); fstat1
1-pf(fstat1,anodev1$Df[2],anodev1$`Resid. Df`[2])

if(!require(multcomp)){install.packages(multcomp)}; library(multcomp)

summary(glht(mod1, mcp(trt='Tukey')))


#testing - looks like full mortality for a treatment breaks the model ####

z <- rep(c(0,1), times=c(70,30))
a <- rep(c(0,1), each=50)
b <- rep(c(0,1), times=c(20,80))
c <- rep(c(0,1), times=c(0,100))
test <- as.data.frame(rbind(cbind('a',a), cbind('b',b), cbind('c',c), cbind('z',z)))
# test <- as.data.frame(rbind(cbind('c',c), cbind('b',b), cbind('a',a)))
test[,2] <- as.integer(as.vector(test[,2]))
test[,3] <- 1-test[,2]
colnames(test) <- c('trt', 'mort', 'harv')

mod <- glm(cbind(mort,harv)~trt, family=binomial, data=test)
anodev<-anova(mod)
fstat<-(anodev$Dev[2]/anodev$Df[2])/(anodev[2,"Resid. Dev"]/anodev[2,"Resid. Df"])
1-pf(fstat,anodev$Df[2],anodev$`Resid. Df`[2])
summary(glht(mod, mcp(trt='Tukey')))
