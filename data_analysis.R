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
#######################Conversion factors are estimated as one parameter###################

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



###########################################################################################
###########################################################################################
