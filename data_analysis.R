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
                 plot.title=element_text(size=28, face="bold", hjust=-0.05))

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
conversions$Reference<-c("Henriksen et al., 1996")
conversions$Organism<-c("Penicillium chrysogenum")

#Hanegraaf and Muller 2001 presented different data
#Paracoccus denitrificans
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.57*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Hanegraaf and Muller, 2001"),
                                           Organism=rep("Paracoccus denitrificans")))

#E. Coli
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.82*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Hanegraaf and Muller, 2001"),
                                           Organism=rep("Escherichia coli")))
#van Duuren et al., 2013
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.553*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("van Duuren et al., 2013"),
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
                                           Reference=rep("Beck at al., 2018"),
                                           Organism=rep("Alicyclobacillus acidocaldarius")))
#Synechococcus 7002
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$Prot.in/0.272*0.45/12.01/4,
                                           Proxy=rep("Cellular protein"),
                                           Reference=rep("Beck at al., 2018"),
                                           Organism=rep("Synechococcus 7002")))


#Make a graph
ggplot(conversions, aes(Reference, Cmic))+geom_boxplot(cex=0.8, aes(colour=Organism), show.legend = F)+
  facet_grid(.~Proxy)+coord_flip()+theme_min+theme(axis.title.y = element_blank())+
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
#1-pchisq(-2*(-ll_model1--ll_model2), df=number of parameters difference)

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


stopImplicitCluster()

mend_i1$goodness
mend_i2$goodness
mend_i3$goodness
mend_i4$goodness
###############################################################################################
###########################################DEB model###########################################
###############################################################################################
source("../deb_i.R")

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

stopImplicitCluster()

deb_i1$goodness
deb_i2$goodness
deb_i3$goodness
deb_i4$goodness

###########################################################################################
###########################################################################################
deb_pars<-as.data.frame(rbind(deb_i4.2[[1]]$pars, deb_i4.2[[2]]$pars, deb_i4.2[[3]]$pars,
                              deb_i4.2[[4]]$pars, deb_i4.2[[5]]$pars, deb_i4.2[[6]]$pars))
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


