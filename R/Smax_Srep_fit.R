#========================================================
#Fitting Ricker models with fixed Smax and fixed Srep
#Author: Catarina wor
#Date February 2022
#========================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)
library(cowplot)
library(tmbstan)
library(reshape)

model_dir <- paste(here::here(), "/TMB", sep="")

#data
dat <- read.csv("../data/Harrison_simples_survchi.csv")



#simple model for parameter guesses
srm <- lm(log(dat$R/dat$S_adj)~ dat$S_adj)
a_srm<-srm$coefficients[1]
b_srm<--srm$coefficients[2]


#========================================
#Model with fixed Smax

compile("../TMB/Ricker_tva_Smax.cpp")
dyn.load(dynlib("../TMB/Ricker_tva_Smax"))

Smaxdata<-list(obs_logR=log(dat$R),obs_S=dat$S_adj, prbeta1=3,prbeta2=3)


parameters<- list(
  alphao=a_srm,
  logSmax = log(1/b_srm),
  rho=.2,
  logvarphi= 0.1,
  alpha=rep(0.9,length(dat$R))
  )




obj <- MakeADFun(Smaxdata,parameters,DLL="Ricker_tva_Smax",random="alpha")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()

opt <- nlminb(obj$par,obj$fn,obj$gr)
#,lower=c(0.1,log(min(dat$R)*.5),0.0,-3.0,rep(0.1,length(dat$R))),        
#              upper=c(4.0,log(1000000),1.0,5.0,rep(5.0,length(dat$R))),)
opt
rep <- obj$report()
rep

SmaxAIC<-2*3-2*-opt$objective


#look at posteriors
fitmcmc1 <- tmbstan(obj, chains=3,
              iter=10000, init="random",
              lower=c(0.1,log(min(dat$R)*.5),0.0,-3.0,rep(0.1,length(dat$R))),        
              upper=c(4.0,log(5000000),1.0,5.0,rep(4.0,length(dat$R))),
               control = list(adapt_delta = 0.98))

mc <- extract(fitmcmc1, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)
    
fit_summary <- summary(fitmcmc1)

posterior <- as.array(fitmcmc1)

mainrun <- melt(posterior)



pm <- ggplot(mainrun)
pm <- pm + geom_density(aes(x=value, color=chains), size=1.2)
pm <- pm + facet_wrap(~parameters, scales="free")
pm <- pm + theme_bw(16)+labs(colour = "Prior")
pm <- pm + scale_color_viridis_d(end = 0.8,option = "A")
pm

predSmax<-matrix(NA,nrow=length(dat$S_adj),ncol=length(rep$alpha))

s_pred<-seq(0,max(dat$S_adj)*1.05,length=length(dat$S_adj))

for(i in 1:length(rep$alpha)){
  predSmax[,i]<- s_pred*exp(rep$alpha[i]-rep$beta*s_pred)
}



df1<-data.frame(S=rep(dat$S_adj,length(rep$alpha)),
    predS=rep(s_pred,length(rep$alpha)),
    predR=c(predSmax),
    R=rep(dat$R,length(rep$alpha)),
    a=rep(rep$alpha,each=length(rep$alpha)),
    umsy=rep(rep$umsy,each=length(rep$umsy)),
    ayr=as.factor(rep(dat$BroodYear,each=length(rep$alpha))),
    model="Smax fixed",
    BroodYear=rep(dat$BroodYear,length(rep$alpha)))




df1<-df1[sort(order(df1$S)),]


p <- ggplot(df1)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_line(aes(x=predS,y=predR, color=ayr), size=2, alpha=0.6)
p <- p + geom_text(aes(x=S,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + theme_bw(16) + scale_color_viridis_d(end = 0.9, option="D")
p <- p + coord_cartesian(xlim=c(0,max(df$S)/1000*1.05))
p <- p + labs(title = "Recursive Bayes Ricker model", x = "Spawners (1,000s)", y = "Recruits (1,000s)", color = "Year\n") 
p



dfa<-data.frame(broodyear=dat$BroodYear,alpha=(rep$alpha),model="Smax")


pa <- ggplot(dfa)
pa <- pa + geom_line(aes(x=broodyear,y=alpha), size=2)
pa <- pa + geom_point(aes(x=broodyear,y=alpha),size=4)
pa <- pa + theme_bw(16)
pa <- pa + labs(title = "Kalman Filter model - alpha time series", y = expression(alpha), x = "Brood year") 
pa



#====================================================
#Model with fixed Srep


compile("../TMB/Ricker_tva_Srep.cpp")
dyn.load(dynlib("../TMB/Ricker_tva_Srep"))


Srepdata<-Smaxdata


parametersSrep<- list(
  alphao=a_srm,
  logSrep = log(a_srm/b_srm),
  rho=.2,
  logvarphi= 0.1,
  alpha=rep(a_srm,length(dat$R))
  )


objSrep <- MakeADFun(Srepdata,parametersSrep,DLL="Ricker_tva_Srep",random="alpha")

newtonOption(objSrep, smartsearch=FALSE)
objSrep$fn()
objSrep$gr()

optSrep <- nlminb(objSrep$par,objSrep$fn,objSrep$gr)
repSrep <- objSrep$report()
repSrep

SrepAIC<-2*3-2*-optSrep$objective


predSrep<-matrix(NA,nrow=length(dat$S_adj),ncol=length(repSrep$alpha))

for(i in 1:length(repSrep$alpha)){
  predSrep[,i]<- s_pred*exp(repSrep$alpha[i]-repSrep$beta[i]*s_pred)
}



df2<-data.frame(S=rep(dat$S_adj,length(repSrep$alpha)),
    predS=rep(s_pred,length(repSrep$alpha)),
    predR=c(predSrep),
    R=rep(dat$R,length(repSrep$alpha)),
    a=rep(repSrep$alpha,each=length(repSrep$alpha)),
    umsy=rep(repSrep$umsy,each=length(repSrep$umsy)),
    ayr=as.factor(rep(dat$BroodYear,each=length(repSrep$alpha))),
    model="Srep fixed",
    BroodYear=rep(dat$BroodYear,length(repSrep$alpha)))


df <- rbind(df1,df2)


p <- ggplot(df)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_line(aes(x=predS/1000,y=predR/1000, color=ayr), size=2, alpha=0.6)
p <- p + geom_text(aes(x=S/1000,y=R/1000,label=BroodYear ),hjust=0, vjust=0)
p <- p + geom_abline(slope=1)
p <- p + theme_bw(16) + scale_color_viridis_d(end = 0.9, option="D")
p <- p + facet_wrap(~model)
p <- p + coord_cartesian(xlim=c(0,max(df$S)/1000*1.05))
p <- p + labs(title = "Recursive Bayes Ricker model", x = "Spawners (1,000s)", y = "Recruits (1,000s)", color = "Year\n") 
p





dfa2<-data.frame(broodyear=dat$BroodYear,alpha=(repSrep$alpha), model="Srep")

da<-rbind(dfa,dfa2)

pa <- ggplot(da)
pa <- pa + geom_line(aes(x=broodyear,y=alpha), size=2)
pa <- pa + geom_point(aes(x=broodyear,y=alpha),size=4)
pa <- pa + theme_bw(16)
pa <- pa + facet_wrap(~model)
pa <- pa + labs(title = "recursive Bayes model - alpha time series", y = expression(alpha), x = "Brood year") 
pa




