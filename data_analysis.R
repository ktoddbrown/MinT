#################################################MinT#################################################
##libraries
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
#library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
library(openxlsx)
library(DEoptim)
library(gridExtra)

##ggplot theme
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=14, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=14, colour="black"),
                 axis.title=element_text(size=14, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=16, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=16, face="bold", hjust=-0.05))

#######################################################################################################

##data 
#respiration rates in umol/ml/h
resp<-read.xlsx(xlsxFile = c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/data_checked_respiration.xlsx"),
                sheet = 1)

#arrange the data
resp.ordered<-resp[order(resp$Sample, resp$Day), ]

#biological data
biology<-read.xlsx(xlsxFile = c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/data_checked_biology.xlsx"),
                   sheet = 1)

#arrange the data in the same way
biology.ordered<-biology[order(biology$Sample, biology$Day), -c(2,3)]
biology.ordered$Sample<-paste0("CSub3_", biology.ordered$Sample)

#merging both data sets
m0<-merge(resp.ordered, biology.ordered, by.x=c("Sample", "Day"), by.y = c("Sample", "Day"), all=T)

####################################proxy based Biomass calculations##################################
#There are plenty of possible conversions
#I will create a new data frame to store them all at one place for later use

conversions<-data.frame(Sample=m0[,1])

#1. protein based conversions
#Henriksen et al., 1996
conversions$Cmic<-m0$Prot.in/12.01/4
conversions$Proxy<-c("Cellular protein")
conversions$Reference<-c("Henriksen et al., \n 1996")
conversions$Organism<-c("Penicillium chrysogenum")

#Hanegraaf and Muller 2001 presented different data
#Paracoccus denitrificans
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.57*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Hanegraaf and \n Muller, 2001"),
                                           Organism=rep("Paracoccus denitrificans")))

#E. Coli
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.82*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Hanegraaf and \n Muller, 2001"),
                                           Organism=rep("Escherichia coli")))
#van Duuren et al., 2013
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.553*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("van Duuren et al., \n 2013"),
                                           Organism=rep("Pseudomonas putida")))
#Baart et al., 2008
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.688*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Baart et al., 2008"),
                                           Organism=rep("Neisseria meningitidis")))
#Beck at al., 2018
#Alicyclobacillus acidocaldarius
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.385*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Beck et al., 2018"),
                                           Organism=rep("Alicyclobacillus acidocaldarius")))
#Synechococcus 7002
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.272*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Beck et al., 2018"),
                                           Organism=rep("Synechococcus 7002")))


#Make a graph
ggplot(conversions, aes(Reference, Cmic))+geom_boxplot(cex=0.8, aes(fill=Organism), show.legend = T)+
  theme_min+theme(axis.title.x = element_blank(), axis.text.x=element_text(size=14, colour="black"), legend.position = c(0.8,0.8))+
  ylab(expression(paste("Microbial biomass carbon (", mu, "mol ", ml^{-1},")")))

###################################Data for modelling#########################################

d<-subset(m0, Substrate=="Glucose" | Substrate=="Cellobiose")
summary(d)

#removing outliers
#respiration rate
ggplot(d, aes(Time, r))+geom_point(cex=6)+facet_wrap(Structure~Substrate, scales="free")

d[(d$Substrate=="Glucose" & d$Structure=="Glass wool" & d$Time>20 & d$Time<30 & d$r<0.4 & !is.na(d$r)), "r"]<-NA

#Proteins
ggplot(d, aes(Time, Prot.in))+geom_point(cex=6)+facet_wrap(Structure~Substrate, scales="free")

d[(d$Substrate=="Cellobiose" & d$Structure=="Broth" &  d$Time>75), "Prot.in"]<-NA
d[(d$Substrate=="Glucose" & d$Structure=="Broth" & d$Time>75), "Prot.in"]<-NA
d[(d$Substrate=="Cellobiose" & d$Structure=="Mixed glass" & d$Time>75), "Prot.in"]<-NA
d[(d$Substrate=="Cellobiose" & d$Structure=="Glass wool" & d$Time>75 & d$Prot.in>110 & !is.na(d$Prot.in)), "Prot.in"]<-NA

###########################Use literature data to calculate biomass###########################
#Calculating Cmic from intercellular protein assuming 0.548% of biomass
d$Cmic<-d$Prot.in/0.548*0.45/12.01/4

#DEB require knowing only protein C concentration
d$Protc<-d$Prot.in*0.46/12.01/4

d$DNAc<-d$DNA*0.51/12.01/4
#############################################################################################
#First, all models are calibrated against r only

#############################################################################################


###############################################################################################
#####################################Monod growth##############################################
###############################################################################################
source("../monod_r.R")

#across all treatments
monod_r1<-monod_r(data=d, FACT = 4)
monod_r1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod_r2<-monod_r(data=d, FACT = 2)
monod_r2$goodness

#for different substrates
monod_r3<-monod_r(data=d, FACT = 1)
monod_r3$goodness

#for each separately 
monod_r4<-monod_r(data=d, FACT = 3)
monod_r4$goodness


stopImplicitCluster()

monod_r1$goodness
monod_r2$goodness
monod_r3$goodness
monod_r4$goodness

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(as.numeric(monod_r1$goodness[1,4])-as.numeric(monod_r2$goodness[1,4])), df=10)

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
source("../mend_r.R")

#across all treatments
mend_r1<-mend_r(data=d, FACT = 4)
mend_r1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_r2<-mend_r(data=d, FACT = 2)
mend_r2$goodness

#for different substrates
mend_r3<-mend_r(data=d, FACT = 1)
mend_r3$goodness

#for each separately 
mend_r4<-mend_r(data=d, FACT = 3)
mend_r4$goodness


stopImplicitCluster()

mend_r1$goodness
mend_r2$goodness
mend_r3$goodness
mend_r4$goodness
###############################################################################################
###########################################DEB model###########################################
###############################################################################################
source("../deb_r.R")

#across all treatments
deb_r1<-deb_r(data=d, FACT = 4)
deb_r1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
deb_r2<-deb_r(data=d, FACT = 2)
deb_r2$goodness

#for different substrates
deb_r3<-deb_r(data=d, FACT = 1)
deb_r3$goodness

#for each separately 
deb_r4<-deb_r(data=d, FACT = 3)
deb_r4$goodness


stopImplicitCluster()

deb_r1$goodness
deb_r2$goodness
deb_r3$goodness
deb_r4$goodness
####################################################################################################
####################################################################################################
####################################################################################################
#All models are further calibrated against the respiration rate and protein concentration in cell

###############################################################################################
#####################################Monod growth##############################################
###############################################################################################
source("../monod.R")

#across all treatments
monod1<-monod(data=d, FACT = 4)
monod1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod2<-monod(data=d, FACT = 2)
monod2$goodness

#for different substrates
monod3<-monod(data=d, FACT = 1)
monod3$goodness

#for each separately 
monod4<-monod(data=d, FACT = 3)
monod4$goodness


stopImplicitCluster()

monod1$goodness
monod2$goodness
monod3$goodness
monod4$goodness

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-ll_model1--ll_model2), df=number of parameters difference)

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
source("../mend.R")

#across all treatments
mend1<-mend(data=d, FACT = 4)
mend1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend2<-mend(data=d, FACT = 2)
mend2$goodness

#for different substrates
mend3<-mend(data=d, FACT = 1)
mend3$goodness

#for each separately 
mend4<-mend(data=d, FACT = 3)
mend4$goodness


stopImplicitCluster()

mend1$goodness
mend2$goodness
mend3$goodness
mend4$goodness
###############################################################################################
###########################################DEB model###########################################
###############################################################################################
source("../deb.R")

#across all treatments
deb1<-deb(data=d, FACT = 4)
deb1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
deb2<-deb(data=d, FACT = 2)
deb2$goodness

#for different substrates
deb3<-deb(data=d, FACT = 1)
deb3$goodness

#for each separately 
deb4<-deb(data=d, FACT = 3)
deb4$goodness


stopImplicitCluster()

deb1$goodness
deb2$goodness
deb3$goodness
deb4$goodness

#`````````````````````````````````````````````````````````````````````````````````````````#
###################Abundance of proteins in cell is estimated as one parameter#############

#`````````````````````````````````````````````````````````````````````````````````````````#
###############################################################################################
#####################################Monod growth##############################################
###############################################################################################
source("../monod_i.R")
source("../monod_i_alternative.R")

#across all treatments
monod_i1<-monod_i(data=d, FACT = 4)
monod_i1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod_i2<-monod_i(data=d, FACT = 2)
monod_i2$goodness

#for different substrates
monod_i3<-monod_i(data=d, FACT = 1)
monod_i3$goodness

#for each separately 
monod_i4<-monod_i(data=d, FACT = 3)
monod_i4$goodness

monod_i4_alternative<-monod_i_alternative(data=d, FACT = 3)
monod_i4_alternative$goodness


stopImplicitCluster()

monod_i1$goodness
monod_i2$goodness
monod_i3$goodness
monod_i4$goodness

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-ll_model1--ll_model2), df=number of parameters difference)

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
source("../mend_i.R")
source("../mend_i_alternative.R")


#across all treatments
mend_i1<-mend_i(data=d, FACT = 4)
mend_i1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_i2<-mend_i(data=d, FACT = 2)
mend_i2$goodness

#for different substrates
mend_i3<-mend_i(data=d, FACT = 1)
mend_i3$goodness

#for each separately 
mend_i4<-mend_i(data=d, FACT = 3)
mend_i4$goodness

mend_i4_alternative<-mend_i_alternative(data=d, FACT = 3)
mend_i4_alternative$goodness


stopImplicitCluster()

mend_i1$goodness
mend_i2$goodness
mend_i3$goodness
mend_i4$goodness
###############################################################################################
###########################################DEB model###########################################
###############################################################################################
source("../deb_i.R")
source("../deb_i_alternative.R")

#across all treatments
deb_i1<-deb_i(data=d, FACT = 4)
deb_i1$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
deb_i2<-deb_i(data=d, FACT = 2)
deb_i2$goodness

#for different substrates
deb_i3<-deb_i(data=d, FACT = 1)
deb_i3$goodness

#for each separately 
deb_i4<-deb_i(data=d, FACT = 3)
deb_i4$goodness

deb_i4_alternative<-deb_i_alternative(data=d, FACT = 3)
deb_i4_alternative$goodness

stopImplicitCluster()

deb_i1$goodness
deb_i2$goodness
deb_i3$goodness
deb_i4$goodness

###########################################################################################
###########################################################################################
deb_pars<-as.data.frame(rbind(deb_i4[[1]]$pars, deb_i4[[2]]$pars, deb_i4[[3]]$pars,
                              deb_i4[[4]]$pars, deb_i4[[5]]$pars, deb_i4[[6]]$pars))
deb_pars$Substrate<-c(rep("Cellobiose", times=3),
                       rep("Glucose", times=3))

deb_pars$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)
Deb_pars<-melt(deb_pars, id.vars=c("Substrate", "Structure"))

deb_pars_sd<-as.data.frame(rbind(summary(deb_i4[[1]]$par_prof)[2,], summary(deb_i4[[2]]$par_prof)[2,], 
                                  summary(deb_i4[[3]]$par_prof)[2,], summary(deb_i4[[4]]$par_prof)[2,], 
                                  summary(deb_i4[[5]]$par_prof)[2,], summary(deb_i4[[6]]$par_prof)[2,]))
deb_pars_sd$Substrate<-c(rep("Cellobiose", times=3),
                          rep("Glucose", times=3))

deb_pars_sd$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)

Deb_pars$sd<-melt(deb_pars_sd, id.vars=c("Substrate", "Structure"))[,4]


ggplot(Deb_pars, aes(Substrate, value))+geom_point(cex=6, aes(colour=Structure))+
  facet_wrap(~variable, scales="free")+geom_errorbar(aes(ymax=value+sd, ymin=value-sd, colour=Structure))

mean(deb_pars$R_0)
mean(deb_pars$S_0)

source("../deb_i_fix.R")

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

deb_i4_fix<-deb_i_fix(data=d,FACT = 3)
deb_i4_fix$goodness

stopImplicitCluster()

#Show the parameters

deb_fpars<-as.data.frame(rbind(deb_i4_fix[[1]]$pars, deb_i4_fix[[2]]$pars, deb_i4_fix[[3]]$pars,
                               deb_i4_fix[[4]]$pars, deb_i4_fix[[5]]$pars, deb_i4_fix[[6]]$pars))
deb_fpars$Substrate<-c(rep("Cellobiose", times=3),
                      rep("Glucose", times=3))

deb_fpars$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)
Deb_fpars<-melt(deb_fpars, id.vars=c("Substrate", "Structure"))

deb_fpars_sd<-as.data.frame(rbind(summary(deb_i4_fix[[1]]$par_prof)[2,], summary(deb_i4_fix[[2]]$par_prof)[2,], 
                                 summary(deb_i4_fix[[3]]$par_prof)[2,], summary(deb_i4_fix[[4]]$par_prof)[2,], 
                                 summary(deb_i4_fix[[5]]$par_prof)[2,], summary(deb_i4_fix[[6]]$par_prof)[2,]))
deb_fpars_sd$Substrate<-c(rep("Cellobiose", times=3),
                         rep("Glucose", times=3))

deb_fpars_sd$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)

Deb_fpars$sd<-melt(deb_fpars_sd, id.vars=c("Substrate", "Structure"))[,4]


ggplot(Deb_fpars, aes(Substrate, value))+geom_point(cex=6, aes(colour=Structure))+
  facet_wrap(~variable, scales="free")+geom_errorbar(aes(ymax=value+sd, ymin=value-sd, colour=Structure))


#keep m0 fixed
mean(deb_fpars$m0)

source("../deb_i_fix2.R")

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

deb_i4_fix2<-deb_i_fix2(data=d,FACT = 3)
deb_i4_fix2$goodness

stopImplicitCluster()

#Show the parameters

deb_fpars2<-as.data.frame(rbind(deb_i4_fix2[[1]]$pars, deb_i4_fix2[[2]]$pars, deb_i4_fix2[[3]]$pars,
                               deb_i4_fix2[[4]]$pars, deb_i4_fix2[[5]]$pars, deb_i4_fix2[[6]]$pars))
deb_fpars2$Substrate<-c(rep("Cellobiose", times=3),
                       rep("Glucose", times=3))

deb_fpars2$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)
Deb_fpars2<-melt(deb_fpars2, id.vars=c("Substrate", "Structure"))

deb_fpars2_sd<-as.data.frame(rbind(summary(deb_i4_fix2[[1]]$par_prof)[2,], summary(deb_i4_fix2[[2]]$par_prof)[2,], 
                                  summary(deb_i4_fix2[[3]]$par_prof)[2,], summary(deb_i4_fix2[[4]]$par_prof)[2,], 
                                  summary(deb_i4_fix2[[5]]$par_prof)[2,], summary(deb_i4_fix2[[6]]$par_prof)[2,]))
deb_fpars2_sd$Substrate<-c(rep("Cellobiose", times=3),
                          rep("Glucose", times=3))

deb_fpars2_sd$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)

Deb_fpars2$sd<-melt(deb_fpars2_sd, id.vars=c("Substrate", "Structure"))[,4]


ggplot(Deb_fpars2, aes(Substrate, value))+
  geom_point(cex=6, aes(colour=Structure), position = position_dodge(width = 1))+
  facet_wrap(~variable, scales="free")+theme_min+
  geom_errorbar(aes(ymax=value+sd, ymin=value-sd, colour=Structure), 
                position = position_dodge(width = 1))


################################################################################################
#######################################Figures##################################################
################################################################################################
deb_fpars2

#deb model
deb_hard<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    #Carbon uptake
    Cu=Vmax*C*S/(Km+C)
    
    #maintnance
    m=S*0.003121921
    
    #reserve mobilization rate for growth and enzyme production
    an=f*R-m
    
    #calibrated variables
    #protein abundance in reserves and structures is mean value reported 
    #by Hanegraaf and Muller, 2001
    r=m+pmax(an*(1-Yu), 0)-pmin(0, an)
    
    #Protc=0.7095*S+0.6085*R
    Protc=fpr*R+fps*S
    
    dR<-Cu-f*R
    dS<-pmax(an*Yu,0)+pmin(0, an)
    dC<--Cu
    
    
    return(list(c(dR, dS, dC), r=r, Protein=Protc))
    
  })
}

out<-as.data.frame(ode(y=c(R=0.02668775, S=0.4782267, C=25), 
                       parms=deb_fpars2[6, c(1:6)], 
         deb_hard, times=seq(0,170)))


out2<-out
ggplot(out, aes(R, r))+geom_point()

out2$Mic<-out2$R+out2$S
Out2<-melt(out2[, c(1,2,3,6,7)], id.vars = c("time"))

ggplot(Out2, aes(time, value))+geom_line(lwd=1.2, aes(colour=variable))+theme_min+
  xlab("Time")+ylab(expression(paste("Microbial Biomass Constituents (", mu, "mol ", ml^{-1}, ")")))+
  theme(legend.position = c(0.8,0.4),
        legend.title = element_blank())+
  ggtitle("a)")
  

out3<-data.frame(time=out2$time, Conv=out2$Protein/out2$Mic)


ggplot(out3, aes(time, Conv))+geom_line(lwd=1.2)+theme_min+
  ggtitle("b)")+xlab("Time")+ylab("Conversion Factor")+ylim(0,0.4)

grid.arrange(ggplot(Out2, aes(time, value))+geom_line(aes(colour=variable))+theme_min+
               xlab("Time")+ylab("Microbial Biomass Pools")+
               theme(legend.position = c(0.8,0.3),
                     legend.title = element_blank())+
               ggtitle("a)"),
             ggplot(Out3, aes(time, Yield))+geom_line(aes(colour=Method))+theme_min+
               ggtitle("b)")+xlab("Time")+theme(legend.position = c(0.8,0.8))+ylim(0,0.7),
             ncol=2)

#####all conversion factors
#first one
out<-as.data.frame(ode(y=c(R=0.02668775, S=0.4782267, C=25), 
                       parms=deb_fpars2[1, c(1:6)], 
                       deb_hard, times=seq(0,125)))
out$Mic<-out$R+out$S
out1<-data.frame(time=out$time, Conv=out$Protein/out$Mic, 
                 Substrate=rep(deb_fpars2[1, 7], times=nrow(out)),
                 Structure=rep(deb_fpars2[1, 8], times=nrow(out)))

#the rest
for(i in 2:nrow(deb_fpars2)){
  out<-as.data.frame(ode(y=c(R=0.02668775, S=0.4782267, C=25), 
                         parms=deb_fpars2[i, c(1:6)], 
                         deb_hard, times=seq(0,125)))
  out$Mic<-out$R+out$S
  out.else<-data.frame(time=out$time, Conv=out$Protein/out$Mic, 
                   Substrate=rep(deb_fpars2[i, 7], times=nrow(out)),
                   Structure=rep(deb_fpars2[i, 8], times=nrow(out)))
  
  out1<-rbind(out1, out.else)
}

ggplot(subset(out1, time<60), aes(time, Conv))+geom_line(lwd=1.2, aes(colour=Structure, linetype=Substrate))+theme_min+
  xlab("Time")+ylab("Conversion Factor")

#all all
#first one
alal<-as.data.frame(ode(y=c(R=0.02668775, S=0.4782267, C=25), 
                       parms=deb_fpars2[1, c(1:6)], 
                       deb_hard, times=seq(0,125)))
alal$Mic<-alal$R+alal$S
alal$Substrate<-rep(deb_fpars2[1, 7], times=nrow(alal))
alal$Structure<-rep(deb_fpars2[1, 8], times=nrow(alal))


#the rest
for(i in 2:nrow(deb_fpars2)){
  alal2<-as.data.frame(ode(y=c(R=0.02668775, S=0.4782267, C=25), 
                          parms=deb_fpars2[i, c(1:6)], 
                          deb_hard, times=seq(0,125)))
  alal2$Mic<-alal2$R+alal2$S
  alal2$Substrate<-rep(deb_fpars2[i, 7], times=nrow(alal2))
  alal2$Structure<-rep(deb_fpars2[i, 8], times=nrow(alal2))
  
  alal<-rbind(alal, alal2)
}

ggplot(alal, aes(Mic, r))+geom_point(cex=4, aes(colour=Structure, shape=Substrate))+theme_min+
  ylab(expression(paste("Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Microbial Biomass (", mu, "mol ", ml^{-1}, ")")))

ggplot(alal, aes(S, r))+geom_point(cex=4, aes(colour=Structure, shape=Substrate))+theme_min+
  ylab(expression(paste("Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Structures (", mu, "mol ", ml^{-1}, ")")))

ggplot(alal, aes(R, r))+geom_point(cex=4, aes(colour=Structure, shape=Substrate))+theme_min+
  ylab(expression(paste("Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Reserves (", mu, "mol ", ml^{-1}, ")")))

#############################################OvP MONOD MEND R################################
#Monod (4.5 to 6.5 inches - landscape)
ggplot(subset(monod_r4$OvP, variable=="r"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  ylim(0,1)+xlim(0,1)+geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Measured Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))

#Mend (4.5 to 6.5 inches - landscape)
ggplot(subset(mend_r4$OvP, variable=="r"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  ylim(0,1)+xlim(0,1)+geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Measured Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))

#############################################OvP MONOD MEND Cmic################################
#Monod (4.5 to 6.5 inches - landscape)
ggplot(subset(monod_r4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#Mend (4.5 to 6.5 inches - landscape)
ggplot(subset(mend_r4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#############################################OvP MONOD MEND Cmic2################################
#Monod (4.5 to 6.5 inches - landscape)
ggplot(subset(monod4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#Mend (4.5 to 6.5 inches - landscape)
ggplot(subset(mend4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#############################################OvP MONOD MEND Cmic3################################
#Monod (4.5 to 6.5 inches - landscape)
ggplot(subset(monod_i4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#Mend (4.5 to 6.5 inches - landscape)
ggplot(subset(mend_i4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#############################################OvP DEB r Cmic################################
#CO2 (4.5 to 6.5 inches - landscape)
ggplot(subset(deb_i4$OvP, variable=="r"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+ylim(0,1)+xlim(0,1)+
  ylab(expression(paste("Predicted Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Measured Respiration Rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))

#Protc (4.5 to 6.5 inches - landscape)
ggplot(subset(deb_i4$OvP, variable=="Protc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+ylim(0.4,1.4)+xlim(0.4,1.4)+
  ylab(expression(paste("Predicted Cellular Protein (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured Cellular Protein (", mu, "mol ", ml^{-1}, ")")))

#DNAc (4.5 to 6.5 inches - landscape)
ggplot(subset(deb_i4_alternative$OvP, variable=="DNAc"), aes(obs, value))+
  geom_point(cex=4, aes(colour=Structure, shape=Substrate), alpha=0.8)+theme_min+
  geom_abline(intercept = 0, slope=1, lwd=1.2, colour="grey", lty=2)+
  ylab(expression(paste("Predicted DNA content (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Measured DNA content (", mu, "mol ", ml^{-1}, ")")))+
  ylim(0, 0.1)+xlim(0,0.1)


source("../deb_full_test.R")

d_test<-subset(m0, Substrate=="Celluloze")
#DEB require knowing only protein C concentration
d_test$Protc<-d_test$Prot.in*0.46/12.01/4
d_test$E<-d_test$Prot.out*0.46/12.01/4

d_test$DNAc<-d_test$DNA*0.51/12.01/4

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for each separately 
deb_full<-deb_full_test(data=d_test, FACT = 2)
deb_full$goodness

deb_celluloze<-deb_i_fix(data=d_test, FACT = 2)
deb_celluloze$goodness

stopImplicitCluster()

d_test2<-subset(m0, Substrate=="Mix")
#DEB require knowing only protein C concentration
d_test2$Protc<-d_test2$Prot.in*0.46/12.01/4
d_test2$E<-d_test2$Prot.out*0.46/12.01/4

d_test2$DNAc<-d_test2$DNA*0.51/12.01/4

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for each separately 
deb_full<-deb_full_test(data=d_test, FACT = 2)
deb_full$goodness

deb_mix<-deb_i_fix(data=d_test, FACT = 2)
deb_mix$goodness

stopImplicitCluster()
