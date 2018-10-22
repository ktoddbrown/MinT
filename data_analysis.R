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
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
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
                 plot.title=element_text(size=18, face="bold", hjust=-0.05))

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

d2<-subset(m0, Substrate=="Glucose" | Substrate=="Cellobiose")

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
d2$Protc<-d2$Prot.in*0.46/12.01/4

d$DNAc<-d$DNA*0.51/12.01/4
d2$DNAc<-d2$DNA*0.51/12.01/4
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

#Proteins
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

#DNA
source("../monod_i_DNA.R")
source("../monod_i_DNA_fixed.R")
#across all treatments
monod_i1_DNA<-monod_i(data=d, FACT = 4)
monod_i1_DNA$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod_i2_DNA<-monod_i_DNA(data=d, FACT = 2)
monod_i2_DNA$goodness

monod_i2_DNA_fixed<-monod_i_DNA_fixed(data=d, FACT = 2)
monod_i2_DNA_fixed$goodness

#for different substrates
monod_i3_DNA<-monod_i_DNA(data=d, FACT = 1)
monod_i3_DNA$goodness

#for each separately 
monod_i4_DNA<-monod_i_DNA(data=d, FACT = 3)
monod_i4_DNA$goodness

stopImplicitCluster()

monod_i1_DNA$goodness
monod_i2_DNA$goodness
monod_i3_DNA$goodness
monod_i4_DNA$goodness

monod_pars<-as.data.frame(rbind(monod_i4_DNA[[1]]$pars, monod_i4_DNA[[2]]$pars, monod_i4_DNA[[3]]$pars,
                              monod_i4_DNA[[4]]$pars, monod_i4_DNA[[5]]$pars, monod_i4_DNA[[6]]$pars))

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-ll_model1--ll_model2), df=number of parameters difference)

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
source("../mend_i.R")

#Proteins
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

#DNA
source("../mend_i_DNA.R")

#across all treatments
mend_i1_DNA<-mend_i_DNA(data=d, FACT = 4)
mend_i1_DNA$goodness

no_cors<-detectCores()
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_i2_DNA<-mend_i_DNA(data=d, FACT = 2)
mend_i2_DNA$goodness

#for different substrates
mend_i3_DNA<-mend_i_DNA(data=d, FACT = 1)
mend_i3_DNA$goodness

#for each separately 
mend_i4_DNA<-mend_i_DNA(data=d, FACT = 3)
mend_i4_DNA$goodness

stopImplicitCluster()

mend_i1_DNA$goodness
mend_i2_DNA$goodness
mend_i3_DNA$goodness
mend_i4_DNA$goodness

###############################################################################################
###########################################DEB model###########################################
###############################################################################################
source("../deb_i_all.R")
source("../deb_i_all_fix.R")

#across all treatments
deb_i1_all<-deb_i_all(data=d, FACT = 4)
deb_i1_all$goodness

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
deb_i2_all<-deb_i_all(data=d, FACT = 2)
deb_i2_all$goodness

deb_i2_all_fix<-deb_i_all_fix(data=d, FACT = 2)
deb_i2_all_fix$goodness

#fixing Yu
source("../deb_i_all_fix_Yu.R")
deb_i2_all_fix_Yu<-deb_i_all_fix_Yu(data=d, FACT = 2)
deb_i2_all_fix_Yu$goodness

plot(c(deb_i2_all_fix_Yu[[1]]$par_prof$pars[,"fpr"],
       deb_i2_all_fix_Yu[[2]]$par_prof$pars[,"fpr"],
       deb_i2_all_fix_Yu[[3]]$par_prof$pars[,"fpr"])~c(deb_i2_all_fix_Yu[[1]]$par_prof$pars[,"fps"],
                                                        deb_i2_all_fix_Yu[[2]]$par_prof$pars[,"fps"],
                                                        deb_i2_all_fix_Yu[[3]]$par_prof$pars[,"fps"]))

anova(lm(c(deb_i2_all_fix_Yu[[1]]$par_prof$pars[,"Vmax"],
       deb_i2_all_fix_Yu[[2]]$par_prof$pars[,"Vmax"],
       deb_i2_all_fix_Yu[[3]]$par_prof$pars[,"Vmax"])~c(deb_i2_all_fix_Yu[[1]]$par_prof$pars[,"Km"],
                                                       deb_i2_all_fix_Yu[[2]]$par_prof$pars[,"Km"],
                                                       deb_i2_all_fix_Yu[[3]]$par_prof$pars[,"Km"])))

################################################################################################
#for different substrates
deb_i3_all<-deb_i_all(data=d, FACT = 1)
deb_i3_all$goodness

#for each separately 
deb_i4_all<-deb_i_all(data=d, FACT = 3)
deb_i4_all$goodness

deb_i4_all_fix<-deb_i_all_fix(data=d, FACT = 3)
deb_i4_all_fix$goodness



stopImplicitCluster()

deb_i1_all$goodness
deb_i2_all$goodness
deb_i3_all$goodness
deb_i4_all$goodness



# source("../deb_i.R")
# 
# #across all treatments
# deb_i1<-deb_i(data=d, FACT = 4)
# deb_i1$goodness
# 
# no_cors<-detectCores()-1
# cl<-makeCluster(no_cors)
# registerDoParallel(cl)
# 
# #for different structures
# deb_i2<-deb_i(data=d, FACT = 2)
# deb_i2$goodness
# 
# #for different substrates
# deb_i3<-deb_i(data=d, FACT = 1)
# deb_i3$goodness
# 
# #for each separately 
# deb_i4<-deb_i(data=d, FACT = 3)
# deb_i4$goodness
# 
# stopImplicitCluster()
# 
# deb_i1$goodness
# deb_i2$goodness
# deb_i3$goodness
# deb_i4$goodness

###########################################################################################
###########################################################################################




deb_pars<-as.data.frame(rbind(deb_i4_all[[1]]$pars, deb_i4_all[[2]]$pars, deb_i4_all[[3]]$pars,
                              deb_i4_all[[4]]$pars, deb_i4_all[[5]]$pars, deb_i4_all[[6]]$pars))
deb_pars$Substrate<-c(rep("Cellobiose", times=3),
                       rep("Glucose", times=3))

deb_pars$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)
Deb_pars<-melt(deb_pars, id.vars=c("Substrate", "Structure"))

deb_pars_sd<-as.data.frame(rbind(summary(deb_i4_all[[1]]$par_prof)[2,], summary(deb_i4_all[[2]]$par_prof)[2,], 
                                  summary(deb_i4_all[[3]]$par_prof)[2,], summary(deb_i4_all[[4]]$par_prof)[2,], 
                                  summary(deb_i4_all[[5]]$par_prof)[2,], summary(deb_i4_all[[6]]$par_prof)[2,]))
deb_pars_sd$Substrate<-c(rep("Cellobiose", times=3),
                          rep("Glucose", times=3))

deb_pars_sd$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"), times=2)

Deb_pars$sd<-melt(deb_pars_sd, id.vars=c("Substrate", "Structure"))[,4]


ggplot(Deb_pars, aes(Substrate, value))+geom_point(cex=6, aes(colour=Structure))+
  facet_wrap(~variable, scales="free")+geom_errorbar(aes(ymax=value+sd, ymin=value-sd, colour=Structure))

mean(deb_pars$R_0)
mean(deb_pars$S_0)


##############################################################################################
deb_p<-as.data.frame(rbind(deb_i2_all_fix[[1]]$pars, deb_i2_all_fix[[2]]$pars, deb_i2_all_fix[[3]]$pars))

deb_p$Structure<-c("Broth", "Glass wool", "Mixed glass")
Deb_p<-melt(deb_p, id.vars=c("Structure"))

deb_p_sd<-as.data.frame(rbind(summary(deb_i2_all_fix[[1]]$par_prof)[2,], summary(deb_i2_all_fix[[2]]$par_prof)[2,], 
                                 summary(deb_i2_all_fix[[3]]$par_prof)[2,]))
deb_p_sd$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"))

Deb_p$sd<-melt(deb_p_sd, id.vars=c("Structure"))[,3]


ggplot(Deb_p, aes(Structure, value))+geom_point(cex=6)+
  facet_wrap(~variable, scales="free")+geom_errorbar(aes(ymax=value+sd, ymin=value-sd, colour=Structure))

mean(c(deb_i2_all_fix[[1]]$par_prof$pars[,"Yu"],
       deb_i2_all_fix[[2]]$par_prof$pars[,"Yu"],
       deb_i2_all_fix[[3]]$par_prof$pars[,"Yu"]))
sd(c(deb_i2_all_fix[[1]]$par_prof$pars[,"Yu"],
       deb_i2_all_fix[[2]]$par_prof$pars[,"Yu"],
       deb_i2_all_fix[[3]]$par_prof$pars[,"Yu"]))


###

deb_p2<-as.data.frame(rbind(deb_i2_all_fix_Yu[[1]]$pars, deb_i2_all_fix_Yu[[2]]$pars, deb_i2_all_fix_Yu[[3]]$pars))

deb_p2$Structure<-c("Broth", "Glass wool", "Mixed glass")
Deb_p2<-melt(deb_p2, id.vars=c("Structure"))

deb_p_sd2<-as.data.frame(rbind(summary(deb_i2_all_fix_Yu[[1]]$par_prof)[2,], 
                               summary(deb_i2_all_fix_Yu[[2]]$par_prof)[2,], 
                               summary(deb_i2_all_fix_Yu[[3]]$par_prof)[2,]))
deb_p_sd2$Structure<-rep(c("Broth", "Glass wool", "Mixed glass"))

Deb_p2$sd<-melt(deb_p_sd2, id.vars=c("Structure"))[,3]


ggplot(Deb_p2, aes(Structure, value))+geom_point(cex=6)+
  facet_wrap(~variable, scales="free")+geom_errorbar(aes(ymax=value+sd, ymin=value-sd, colour=Structure))

Deb_p2[Deb_p2$variable=="Vmax", c("Structure", "value", "sd")]
Deb_p2[Deb_p2$variable=="Km", c("Structure", "value", "sd")]
Deb_p2[Deb_p2$variable=="f", c("Structure", "value", "sd")]
Deb_p2[Deb_p2$variable=="m0", c("Structure", "value", "sd")]
Deb_p2[Deb_p2$variable=="fpr", c("Structure", "value", "sd")]
Deb_p2[Deb_p2$variable=="fps", c("Structure", "value", "sd")]
Deb_p2[Deb_p2$variable=="fds", c("Structure", "value", "sd")]


################################################################################################
#######################################Figures##################################################
################################################################################################
###Fig. 2 - time course of Proteins and DNA with DEB simulation
f2<-d2[,c(2, 3, 4, 14, 15)]
f2.2<-f2 %>% group_by(Day, Substrate, Structure) %>% summarize(prot.mean=mean(Protc, na.rm = T),
                                                               prot.sd=sd(Protc, na.rm = T), 
                                                               dna.mean=mean(DNAc, na.rm = T),
                                                               dna.sd=sd(DNAc, na.rm = T))
f2.2$Structure<-recode(f2.2$Structure, "Broth" = "BROTH",
                       "Glass wool" = "WOOL",
                       "Mixed glass" = "GLASS")

#deb_model
deb_model<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    Cu=Vmax*C*S/(Km+C)
    m=S*m0
    an=f*R-m
    r=m+pmax(an*(1-0.6944324), 0)
    Protc=fpr*R+fps*S
    DNAc=fds*S
    dR<-Cu-f*R
    dS<-pmax(an*0.6944324,0)+pmin(0, an)
    dC<--Cu
    return(list(c(dR, dS, dC), r=r, Protc=Protc, DNAc=DNAc))
  })
}


#DEB simulation
deb_broth<-as.data.frame(ode(y=c(R=0.02988069, S=0.4958934, C=25), 
                                 parms=c(deb_p2[1,-8]), func=deb_model, times=seq(1,125)))
deb_broth$Structure<-rep("BROTH", times=nrow(deb_broth))

deb_wool<-as.data.frame(ode(y=c(R=0.02988069, S=0.4958934, C=25), 
                             parms=c(deb_p2[2,-8]), func=deb_model, times=seq(1,125)))
deb_wool$Structure<-rep("WOOL", times=nrow(deb_broth))

deb_glass<-as.data.frame(ode(y=c(R=0.02988069, S=0.4958934, C=25), 
                            parms=c(deb_p2[3,-8]), func=deb_model, times=seq(1,125)))
deb_glass$Structure<-rep("GLASS", times=nrow(deb_broth))

deb_f2<-rbind(deb_broth, deb_wool, deb_glass)

####################################################

#Monod simulation
monod_model<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    Protc=fp*Cmic
    dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
    dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
    return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C), Protc=Protc))
    
  })
}

monod_model_DNA<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    DNAc=fd*Cmic
    dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
    dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
    return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C), DNAc=DNAc))
    
  })
}

monod_broth<-as.data.frame(ode(y=c(Cmic=monod_i2[[1]]$pars[["Cmic_0"]], C=25), 
                             parms=monod_i2[[1]]$pars, func=monod_model, times=seq(1,125)))
monod_broth$Structure<-rep("BROTH", times=nrow(monod_broth))

monod_broth$DNAc<-as.data.frame(ode(y=c(Cmic=monod_i2_DNA[[1]]$pars[["Cmic_0"]], C=25), 
                  parms=monod_i2_DNA[[1]]$pars, func=monod_model_DNA, times=seq(1,125)))[,4]

monod_wool<-as.data.frame(ode(y=c(Cmic=monod_i2[[2]]$pars[["Cmic_0"]], C=25), 
                               parms=monod_i2[[2]]$pars, func=monod_model, times=seq(1,125)))
monod_wool$Structure<-rep("WOOL", times=nrow(monod_wool))

monod_wool$DNAc<-as.data.frame(ode(y=c(Cmic=monod_i2_DNA[[2]]$pars[["Cmic_0"]], C=25), 
                                    parms=monod_i2_DNA[[2]]$pars, func=monod_model_DNA, times=seq(1,125)))[,4]

monod_glass<-as.data.frame(ode(y=c(Cmic=monod_i2[[3]]$pars[["Cmic_0"]], C=25), 
                              parms=monod_i2[[3]]$pars, func=monod_model, times=seq(1,125)))
monod_glass$Structure<-rep("GLASS", times=nrow(monod_glass))

monod_glass$DNAc<-as.data.frame(ode(y=c(Cmic=monod_i2_DNA[[3]]$pars[["Cmic_0"]], C=25), 
                                   parms=monod_i2_DNA[[3]]$pars, func=monod_model_DNA, times=seq(1,125)))[,4]

monod_f2<-rbind(monod_broth, monod_wool, monod_glass)
############################################################

#MEND simulation
mend_model<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    F1=1/CUE*(Vmax+mr)*Cmic*C/(Km+C)
    F4=(1/CUE-1)*Vmax*Cmic*C/(Km+C)
    F5=(1/CUE-1)*mr*Cmic*C/(Km+C)
    F8=mr*Cmic
    Protc=fp*Cmic
    dCmic<-F1-(F4+F5)-F8
    dC<--F1+F8
    return(list(c(dCmic, dC), r=F4+F5, Protc=Protc))
    })
}

mend_model_DNA<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    F1=1/CUE*(Vmax+mr)*Cmic*C/(Km+C)
    F4=(1/CUE-1)*Vmax*Cmic*C/(Km+C)
    F5=(1/CUE-1)*mr*Cmic*C/(Km+C)
    F8=mr*Cmic
    DNAc=fd*Cmic
    dCmic<-F1-(F4+F5)-F8
    dC<--F1+F8
    return(list(c(dCmic, dC), r=F4+F5, DNAc=DNAc))
    })
}

mend_broth<-as.data.frame(ode(y=c(Cmic=mend_i2[[1]]$pars[["Cmic_0"]], C=25), 
                               parms=mend_i2[[1]]$pars, func=mend_model, times=seq(1,125)))
mend_broth$Structure<-rep("BROTH", times=nrow(monod_broth))

mend_broth$DNAc<-as.data.frame(ode(y=c(Cmic=mend_i2_DNA[[1]]$pars[["Cmic_0"]], C=25), 
                                    parms=mend_i2_DNA[[1]]$pars, func=mend_model_DNA, times=seq(1,125)))[,4]

mend_wool<-as.data.frame(ode(y=c(Cmic=mend_i2[[2]]$pars[["Cmic_0"]], C=25), 
                              parms=mend_i2[[2]]$pars, func=mend_model, times=seq(1,125)))
mend_wool$Structure<-rep("WOOL", times=nrow(monod_wool))

mend_wool$DNAc<-as.data.frame(ode(y=c(Cmic=mend_i2_DNA[[2]]$pars[["Cmic_0"]], C=25), 
                                   parms=mend_i2_DNA[[2]]$pars, func=mend_model_DNA, times=seq(1,125)))[,4]

mend_glass<-as.data.frame(ode(y=c(Cmic=mend_i2[[3]]$pars[["Cmic_0"]], C=25), 
                               parms=mend_i2[[3]]$pars, func=mend_model, times=seq(1,125)))
mend_glass$Structure<-rep("GLASS", times=nrow(monod_glass))

mend_glass$DNAc<-as.data.frame(ode(y=c(Cmic=mend_i2_DNA[[3]]$pars[["Cmic_0"]], C=25), 
                                    parms=mend_i2_DNA[[3]]$pars, func=mend_model_DNA, times=seq(1,125)))[,4]

mend_f2<-rbind(mend_broth, mend_wool, mend_glass)
############################################################

#7 to 7 inches frame when exported to pdf

#DEB simulation
ggplot(f2.2, aes(x=Day*24))+
  geom_point(cex=6, aes(y=prot.mean, colour="Cellular Proteins", shape=Substrate))+
  geom_errorbar(width=0.1, aes(ymax=prot.mean+prot.sd, ymin=prot.mean-prot.sd, colour="Cellular Proteins"))+
  geom_point(cex=6, aes(y=dna.mean*30, colour="DNA", shape=Substrate))+
  facet_wrap(~Structure, dir="v", strip.position="top")+
  geom_errorbar(width=0.1, aes(ymax=dna.mean*30+dna.sd*30, ymin=dna.mean*30-dna.sd*30, colour="DNA"))+
  geom_line(lwd=0.8, data=deb_f2, aes(x=time, y=Protc, colour="Cellular Proteins"))+
  geom_line(lwd=0.8, data=deb_f2, aes(x=time, y=DNAc*30, colour="DNA"))+
  scale_y_continuous(sec.axis = sec_axis(~./30, 
                                         name = expression(paste("DNA (", mu, "mol " (C[DNA])~ml^{-1}, ")"))))+
  scale_colour_manual(values = c("black", "grey"))+
  labs(y = expression(paste("Cellular Proteins (", mu, "mol " (C[Proteins])~ml^{-1}, ")")),
       x = "Time (h)",
       colour = "Legend")+theme_min+theme(strip.background = element_blank(),
                                          legend.title = element_blank(),
                                          legend.key.height = unit(0.01, "in"),
                                          legend.margin = margin(0,0,0,0, unit="in"),
                                          legend.direction = c("horizontal"),
                                          legend.position = c(0.70,0.58),
                                          axis.title.x = element_text(size=18),
                                          axis.title.y = element_text(size=18),
                                          axis.text.x = element_text(size=18),
                                          axis.text.y = element_text(size=18))+
  scale_x_continuous(expand = c(0,0))

#Monod simulation
ggplot(f2.2, aes(x=Day*24))+
  geom_point(cex=6, aes(y=prot.mean, colour="Cellular Proteins", shape=Substrate))+
  geom_errorbar(width=0.1, aes(ymax=prot.mean+prot.sd, ymin=prot.mean-prot.sd, colour="Cellular Proteins"))+
  geom_point(cex=6, aes(y=dna.mean*30, colour="DNA", shape=Substrate))+
  facet_wrap(~Structure, dir="v", strip.position="top")+
  geom_errorbar(width=0.1, aes(ymax=dna.mean*30+dna.sd*30, ymin=dna.mean*30-dna.sd*30, colour="DNA"))+
  geom_line(lwd=0.8, data=monod_f2, aes(x=time, y=Protc, colour="Cellular Proteins"))+
  geom_line(lwd=0.8, data=monod_f2, aes(x=time, y=DNAc*30, colour="DNA"))+
  scale_y_continuous(sec.axis = sec_axis(~./30, 
                                         name = expression(paste("DNA (", mu, "mol " (C[DNA])~ml^{-1}, ")"))))+
  scale_colour_manual(values = c("black", "grey"))+
  labs(y = expression(paste("Cellular Proteins (", mu, "mol " (C[Proteins])~ml^{-1}, ")")),
       x = "Time (h)",
       colour = "Legend")+theme_min+theme(strip.background = element_blank(),
                                          legend.title = element_blank(),
                                          legend.key.height = unit(0.01, "in"),
                                          legend.margin = margin(0,0,0,0, unit="in"),
                                          legend.direction = c("horizontal"),
                                          legend.position = c(0.70,0.58),
                                          axis.title.x = element_text(size=18),
                                          axis.title.y = element_text(size=18),
                                          axis.text.x = element_text(size=18),
                                          axis.text.y = element_text(size=18))+
  scale_x_continuous(expand = c(0,0))

#MEND simulation
ggplot(f2.2, aes(x=Day*24))+
  geom_point(cex=6, aes(y=prot.mean, colour="Cellular Proteins", shape=Substrate))+
  geom_errorbar(width=0.1, aes(ymax=prot.mean+prot.sd, ymin=prot.mean-prot.sd, colour="Cellular Proteins"))+
  geom_point(cex=6, aes(y=dna.mean*30, colour="DNA", shape=Substrate))+
  facet_wrap(~Structure, dir="v", strip.position="top")+
  geom_errorbar(width=0.1, aes(ymax=dna.mean*30+dna.sd*30, ymin=dna.mean*30-dna.sd*30, colour="DNA"))+
  geom_line(lwd=0.8, data=mend_f2, aes(x=time, y=Protc, colour="Cellular Proteins"))+
  geom_line(lwd=0.8, data=mend_f2, aes(x=time, y=DNAc*30, colour="DNA"))+
  scale_y_continuous(sec.axis = sec_axis(~./30, 
                                         name = expression(paste("DNA (", mu, "mol " (C[DNA])~ml^{-1}, ")"))))+
  scale_colour_manual(values = c("black", "grey"))+
  labs(y = expression(paste("Cellular Proteins (", mu, "mol " (C[Proteins])~ml^{-1}, ")")),
       x = "Time (h)",
       colour = "Legend")+theme_min+theme(strip.background = element_blank(),
                                          legend.title = element_blank(),
                                          legend.key.height = unit(0.01, "in"),
                                          legend.margin = margin(0,0,0,0, unit="in"),
                                          legend.direction = c("horizontal"),
                                          legend.position = c(0.70,0.58),
                                          axis.title.x = element_text(size=18),
                                          axis.title.y = element_text(size=18),
                                          axis.text.x = element_text(size=18),
                                          axis.text.y = element_text(size=18))+
  scale_x_continuous(expand = c(0,0))



############Fig. 3
#monod and mend model parameters as affected by the variables being used for calibration
#monod parameters are monod_i2, monod_i2_DNA and monod_r2
#mend parameters are mend_i2, mend_i2_DNA and mend_r2

monod_r_par<-as.data.frame(rbind(monod_r2[[1]]$pars[-5], 
                                 monod_r2[[2]]$pars[-5], 
                                 monod_r2[[3]]$pars[-5]))
monod_r_par$Structure<-c("BROTH", "WOOL", "GLASS")
monod_r_par$Legend<-c("Respiration")

monod_p_par<-as.data.frame(rbind(monod_i2[[1]]$pars[-c(5,6)], 
                                 monod_i2[[2]]$pars[-c(5,6)], 
                                 monod_i2[[3]]$pars[-c(5,6)]))
monod_p_par$Structure<-c("BROTH", "WOOL", "GLASS")
monod_p_par$Legend<-c("Respiration & Cellular Proteins")

monod_d_par<-as.data.frame(rbind(monod_i2_DNA[[1]]$pars[-c(5,6)], 
                                 monod_i2_DNA[[2]]$pars[-c(5,6)], 
                                 monod_i2_DNA[[3]]$pars[-c(5,6)]))
monod_d_par$Structure<-c("BROTH", "WOOL", "GLASS")
monod_d_par$Legend<-c("Respiration & DNA")


monod_par<-rbind(monod_r_par, monod_p_par, monod_d_par)
Monod_par<-melt(monod_par, id.vars = c("Structure", "Legend"))
Monod_par$Model<-c("Monod")

#parameters differencess
#Vmax
anova(lm(value~Legend+Structure, subset(Monod_par, variable=="Vmax"),
           weights = sd^2))
#Km
anova(lm(value~Structure+Legend, subset(Monod_par, variable=="Km"),
         weights = sd^2))
#CUE
anova(lm(value~Legend+Structure, subset(Monod_par, variable=="CUE"),
         weights = sd^2))
#k
anova(lm(value~Legend+Structure, subset(Monod_par, variable=="k"),
         weights = sd^2))

monod_r_sd<-as.data.frame(rbind(summary(monod_r2[[1]]$par_prof)[2,c(1:4)], 
                                 summary(monod_r2[[2]]$par_prof)[2,c(1:4)], 
                                 summary(monod_r2[[3]]$par_prof)[2,c(1:4)]))
monod_r_sd$Structure<-c("BROTH", "WOOL", "GLASS")
monod_r_sd$Legend<-c("Respiration")

monod_p_sd<-as.data.frame(rbind(summary(monod_i2[[1]]$par_prof)[2,c(1:4)], 
                                summary(monod_i2[[2]]$par_prof)[2,c(1:4)], 
                                summary(monod_i2[[3]]$par_prof)[2,c(1:4)]))
monod_p_sd$Structure<-c("BROTH", "WOOL", "GLASS")
monod_p_sd$Legend<-c("Respiration & Cellular Proteins")

monod_d_sd<-as.data.frame(rbind(summary(monod_i2_DNA[[1]]$par_prof)[2,c(1:4)], 
                                summary(monod_i2_DNA[[2]]$par_prof)[2,c(1:4)], 
                                summary(monod_i2_DNA[[3]]$par_prof)[2,c(1:4)]))
monod_d_sd$Structure<-c("BROTH", "WOOL", "GLASS")
monod_d_sd$Legend<-c("Respiration & DNA")


monod_sd<-rbind(monod_r_sd, monod_p_sd, monod_d_sd)
Monod_par$sd<-melt(monod_sd, id.vars = c("Structure", "Legend"))[,4]

mend_r_par<-as.data.frame(rbind(mend_r2[[1]]$pars[-5], 
                                mend_r2[[2]]$pars[-5], 
                                mend_r2[[3]]$pars[-5]))
mend_r_par$Structure<-c("BROTH", "WOOL", "GLASS")
mend_r_par$Legend<-c("Respiration")

mend_p_par<-as.data.frame(rbind(mend_i2[[1]]$pars[-c(5,6)], 
                                mend_i2[[2]]$pars[-c(5,6)], 
                                mend_i2[[3]]$pars[-c(5,6)]))
mend_p_par$Structure<-c("BROTH", "WOOL", "GLASS")
mend_p_par$Legend<-c("Respiration & Cellular Proteins")

mend_d_par<-as.data.frame(rbind(mend_i2_DNA[[1]]$pars[-c(5,6)], 
                                mend_i2_DNA[[2]]$pars[-c(5,6)], 
                                mend_i2_DNA[[3]]$pars[-c(5,6)]))
mend_d_par$Structure<-c("BROTH", "WOOL", "GLASS")
mend_d_par$Legend<-c("Respiration & DNA")


mend_par<-rbind(mend_r_par, mend_p_par, mend_d_par)
Mend_par<-melt(mend_par, id.vars = c("Structure", "Legend"))
Mend_par$Model<-c("MEND")

mend_r_sd<-as.data.frame(rbind(summary(mend_r2[[1]]$par_prof)[2,c(1:4)], 
                                summary(mend_r2[[2]]$par_prof)[2,c(1:4)], 
                                summary(mend_r2[[3]]$par_prof)[2,c(1:4)]))
mend_r_sd$Structure<-c("BROTH", "WOOL", "GLASS")
mend_r_sd$Legend<-c("Respiration")

mend_p_sd<-as.data.frame(rbind(summary(mend_i2[[1]]$par_prof)[2,c(1:4)], 
                                summary(mend_i2[[2]]$par_prof)[2,c(1:4)], 
                                summary(mend_i2[[3]]$par_prof)[2,c(1:4)]))
mend_p_sd$Structure<-c("BROTH", "WOOL", "GLASS")
mend_p_sd$Legend<-c("Respiration & Cellular Proteins")

mend_d_sd<-as.data.frame(rbind(summary(mend_i2_DNA[[1]]$par_prof)[2,c(1:4)], 
                                summary(mend_i2_DNA[[2]]$par_prof)[2,c(1:4)], 
                                summary(mend_i2_DNA[[3]]$par_prof)[2,c(1:4)]))
mend_d_sd$Structure<-c("BROTH", "WOOL", "GLASS")
mend_d_sd$Legend<-c("Respiration & DNA")


mend_sd<-rbind(mend_r_sd, mend_p_sd, mend_d_sd)
Mend_par$sd<-melt(mend_sd, id.vars = c("Structure", "Legend"))[,4]

Monod_par$Variable<-Monod_par$variable
levels(Monod_par$Variable)<-c("V[MAX]", "K[M]", "CUE", "k[MB]")

f3a<-ggplot(Monod_par, aes(x=Structure, y=value))+
  geom_point(cex=6, aes(shape=Legend, colour=Legend), position = position_dodge(width = 1), show.legend = T)+
  geom_errorbar(width=0.1, aes(ymin=value-sd, ymax=value+sd, colour=Legend),
                position = position_dodge(width = 1), show.legend = T)+
  facet_wrap(~Variable, scales="free", nrow=1, strip.position = "top", labeller = label_parsed)+theme_min+
  scale_shape_manual(values=c(16, 16, 1))+
  scale_colour_manual(values = c("black", "grey", "black"))+
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18, face=c("bold.italic")),
        legend.position = c("bottom"),
        legend.title = element_blank())+
  ylab("Parameter Value")+ggtitle("(a)")

Mend_par$Variable<-Mend_par$variable
levels(Mend_par$Variable)<-c("V[MAX]", "K[M]", "CUE", "m[R]")

f3b<-ggplot(Mend_par, aes(x=Structure, y=value))+
  geom_point(cex=6, aes(shape=Legend, colour=Legend), position = position_dodge(width = 1), show.legend = T)+
  geom_errorbar(width=0.1, aes(ymin=value-sd, ymax=value+sd, colour=Legend),
                position = position_dodge(width = 1), show.legend = T)+
  facet_wrap(~Variable, scales="free", nrow=1, strip.position = "top", labeller = label_parsed)+theme_min+
  scale_shape_manual(values=c(16, 16, 1))+
  scale_colour_manual(values = c("black", "grey", "black"))+
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c("bottom"),
        legend.title = element_blank())+
  ylab("Parameter Value")+ggtitle("(b)")

#8 per 14 inches
grid.arrange(f3a, f3b, nrow=2)

###Fig 4 
#respiration rate as a function of Cmic and R
f4a<-ggplot(deb_f2, aes(R+S, r))+geom_point(cex=6, aes(colour=Structure, shape=Structure), show.legend = F)+
  theme_min+scale_shape_manual(values=c(16, 16, 1))+
  scale_colour_manual(values = c("black", "grey", "black"))+
  xlab(expression(paste(C[MB], " (", mu, "mol ", ml^{-1}, ")")))+
  ylab(expression(paste(R[H], " (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  ggtitle("(a)")
  
f4b<-ggplot(deb_f2, aes(R, r))+geom_point(cex=6, aes(colour=Structure, shape=Structure), 
                                     show.legend = T)+
  theme_min+scale_shape_manual(values=c(16, 16, 1))+
  scale_colour_manual(values = c("black", "grey", "black"))+
  xlab(expression(paste("R (", mu, "mol ", ml^{-1}, ")")))+
  ylab(expression(paste(R[H], " (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  ggtitle("(b)")+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2))

#5 per 10 inches
grid.arrange(f4a, f4b, nrow=1)

###Fig 5
#conversion factors
ggplot(deb_f2, aes(x=time))+
  facet_wrap(~Structure, dir="v", strip.position="top")+
  geom_line(lwd=0.8, aes(y=Protc/(R+S)*100, colour="Cellular Proteins"))+
  geom_line(lwd=0.8, aes(x=time, y=DNAc/(R+S)*15*100, colour="DNA"))+
  scale_y_continuous(sec.axis = sec_axis(~./15, 
                                         name = expression(paste(frac(C[DNA], C[MB])%*%100))))+
  scale_colour_manual(values = c("black", "grey"))+
  labs(y = expression(paste(frac(C[Proteins], C[MB])%*%100)),
       x = "Time (h)",
       colour = "Legend")+theme_min+theme(strip.background = element_blank(),
                                          legend.title = element_blank(),
                                          legend.key.height = unit(0.01, "in"),
                                          legend.margin = margin(0,0,0,0, unit="in"),
                                          legend.direction = c("horizontal"),
                                          legend.position = c(0.50,0.58),
                                          axis.title.x = element_text(size=18),
                                          axis.title.y = element_text(size=18),
                                          axis.text.x = element_text(size=18),
                                          axis.text.y = element_text(size=18),
                                          plot.margin = unit(c(0.05, 0.2, 0.05, 0.05), "in"))+
  scale_x_continuous(expand = c(0,0))


#################################################################################################
#################################################################################################
#################################################################################################
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

######################################community data##########################################
#bacteria
bac<-t(read.xlsx(xlsxFile=c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/CSUB3_NOR16S_whid_ncs (002).xlsx"),
              1, colNames = F))
na<-as.vector(bac[1,])
bac<-bac[-1,]

colnames(bac)<-na
library(vegan)

bac.dat<-bac[, -c(1:4)]
bac.dat<-apply(bac.dat, 2, as.numeric)
bac.env<-bac[, c(1:4)]
bac.env<-as.data.frame(bac.env)
bac.env<-apply(bac.env, 2, as.factor)


bac.dat2<-decostand(bac.dat, method=c("normalize"))

bac.dist<-vegdist(bac.dat2)

bac.aov<-adonis(bac.dat2~Day, bac.env)

bac.aov$aov

#pca analysis and scores extraction

bac.pca<-rda(bac.dat2)

summary(bac.pca)

ordixyplot(bac.pca, cex=3, groups=bac.env$Structure, pch=16)

bac.sc<-scores(bac.pca, choices=c(1:3), display="species")
summary(bac)

#fungi
fu<-t(read.xlsx(xlsxFile=c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/CSUB3_NormITS_whid_NCS.xlsx"),
                 1, colNames = F))
naf<-as.vector(fu[1,])
fu<-fu[-1,]

colnames(fu)<-naf

fu.dat<-fu[, -c(1:4)]
fu.dat<-apply(fu.dat, 2, as.numeric)
fu.env<-fu[, c(1:4)]
fu.env<-as.data.frame(fu.env)
fu.env<-apply(bac.env, 2, as.factor)


fu.dat2<-decostand(fu.dat, method=c("normalize"))

fu.dist<-vegdist(fu.dat2)

fu.aov<-adonis(fu.dat2~Substrate, fu.env)

fu.aov$aov

#pca analysis and scores extraction

fu.pca<-rda(fu.dat)

ordixyplot(fu.pca, cex=3, groups=fu.env$Structure, pch=16)

bac.sc<-scores(bac.pca, choices=c(1:3), display="species")
summary(bac)


##########################################marstrop data#####################################
mar<-read.xlsx(xlsxFile = c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/marstrop.xlsx"),
               2)
source("../deb_mar.R")

deb_mar_res<-deb_mar(data=mar)
deb_mar_res$fit$Gfit
deb_mar_res$pars

ggplot(deb_mar_res$fit$Yhat, aes(time, obs))+geom_point(cex=6)+
  facet_wrap(~variable, scales = "free")+
  geom_line(aes(x=time, y=value), lwd=1.5)
