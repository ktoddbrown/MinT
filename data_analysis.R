#################################################MinT#################################################
##libraries
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
library(data.table)
library(openxlsx)
library(plyr)

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


#2. DNA based conversions
#Christensen et al., 1993 and 1995
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.0264*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Christensen et al., 1993 and 1995"),
                                           Organism=rep("Soil community")))


#Yokoyama et al., 2017
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.098/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Yokoyama et al., 2017"),
                                           Organism=rep("Soil community")))



#Marstrop et al., 2000
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA*2.1458/12.01/4/0.45,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Marstrop et al., 2000"),
                                           Organism=rep("Soil community")))

#Makino et al., 2003
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.0214*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Makino et al., 2003"),
                                           Organism=rep("Escherichia coli")))
#Van Putten et al., 1995
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=(m0$DNA*0.2503+15)*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Van Putten et al., 1995"),
                                           Organism=rep("Bacillus licheniformis")))


#Henriksen et al., 1996
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.0096*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Henriksen et al., 1996"),
                                           Organism=rep("Penicillium chrysogenum")))

#Hanegraaf and Muller, 2001 present data from several bacterial strains
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.045*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Hanegraaf and Muller, 2001"),
                                           Organism=rep("Three different species")))

#van Duuren et al., 2013
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.027*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("van Duuren et al., 2013"),
                                           Organism=rep("Pseudomonas putida")))

#Baart et al., 2008
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.012*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Baart et al., 2008"),
                                           Organism=rep("Neisseria meningitidis")))


#Beck at al., 2018
conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.01*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Beck at al., 2018"),
                                           Organism=rep("Escherichia coli")))

conversions<-rbind(conversions, data.frame(Sample=m0[,1],
                                           Cmic=m0$DNA/0.004*0.45/12.01/4,
                                           Proxy=rep("DNA"),
                                           Reference=rep("Beck at al., 2018"),
                                           Organism=rep("Synechococcus 7002")))

#Make a graph
ggplot(conversions, aes(Reference, Cmic))+geom_boxplot(cex=0.8, aes(colour=Organism), show.legend = F)+
  facet_grid(.~Proxy)+coord_flip()+theme_min+theme(axis.title.y = element_blank())+
  ylab(expression(paste("Microbial biomass carbon (", mu, "mol ", ml^{-1},")")))

######################################################################################################
#carefull here, there is slightly different recalculation
m0$Cmic.dna.init<-m0$DNA.init/0.0494*0.44/12.01/3*0.25
######################################################################################################


#extracellular protein to extracellular protein carbon (46% of carbon in protein - Vrede et al., 2004)
m0$E<-m0$Prot.out*0.46/12.01/4


#extracting only the columns of interest
mint<-m0[,c("Structure", "Substrate", "r", "Time", "DOCinit", "Cmic.dna", "Cmic.dna.init", "E", "Cmic.prot.filled")]
mint<-mint[!is.na(mint$Substrate), ]


dat<-subset(mint, Substrate!="Free" & Substrate!="Celluloze")
summary(dat)

#######################################################################################
##############################First order decay########################################
#######################################################################################
source("../first_order_decay_function.R")

#across all
decay_all_results<-first_order_function(data = dat, SUB = TRUE, FACT = 4, Niter = 10000, Vars = c("Substrate", "r", "Time"))

#no_cors<-detectCores()
#cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
decay_structures_results<-first_order_function(data = dat, SUB = TRUE, FACT = 2, Niter = 10000, Vars = c("Substrate", "r", "Time"))
#for different substrates
decay_substrates_results<-first_order_function(data = dat, SUB = FALSE, FACT = 1, Niter = 10000, Vars = c("Time", "r"))
#for each separately
decay_unique_results<-first_order_function(data = dat, SUB = FALSE, FACT = 3, Niter = 10000, Vars = c("Time", "r"))


stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

decay_all_results$likelihood
decay_structures_results$likelihood
decay_substrates_results$likelihood
decay_unique_results$likelihood

###############################################################################################
#####################################Monod growth##############################################
###############################################################################################
source("../monod_growth_function.R")


#across all treatments
monod_all_results<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Time"),
                      ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                      VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#no_cors<-detectCores()-1
#cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod_structures_results<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Vars = c("Substrate", "r", "Time"),
                                         ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                         VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#for different substrates
monod_substrates_results<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Vars = c("Time", "r"),
                                         ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                         VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

#for each separately 
monod_unique_results<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Vars = c("Time", "r"),
                                         ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                         VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

decay_all_results$likelihood
decay_structures_results$likelihood
decay_substrates_results$likelihood
decay_unique_results$likelihood

monod_all_results$likelihood
monod_structures_results$likelihood
monod_substrates_results$likelihood
monod_unique_results$likelihood

###############################################################################################
#####################################Microbial - enzyme model##################################
###############################################################################################
source("../mem_function.R")


#across all treatments
mem_all_results<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Time"),
                                         ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                         VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#no_cors<-detectCores()-1
#cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mem_structures_results<-mem_function(data=dat, SUB = TRUE, FACT = 2, Vars = c("Substrate", "r", "Time"),
                                                ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                                VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#for different substrates
mem_substrates_results<-mem_function(data=dat, SUB = FALSE, FACT = 1, Vars = c("Time", "r"),
                                                ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                                VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

#for each separately 
mem_unique_results<-mem_function(data=dat, SUB = FALSE, FACT = 3, Vars = c("Time", "r"),
                                            ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                            VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

decay_all_results$likelihood
decay_structures_results$likelihood
decay_substrates_results$likelihood
decay_unique_results$likelihood

monod_all_results$ll_r
monod_structures_results$ll_r
monod_substrates_results$ll_r
monod_unique_results$ll_r

mem_all_results$ll_r
mem_structures_results$ll_r
mem_substrates_results$ll_r
mem_unique_results$ll_r

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
source("../mend_function.R")

#across all treatments
mend_all_results<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Time"),
                              ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                              VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#no_cors<-detectCores()-1
#cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_structures_results<-mend_function(data=dat, SUB = TRUE, FACT = 2, Vars = c("Substrate", "r", "Time"),
                                            ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                            VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#for different substrates
mend_substrates_results<-mend_function(data=dat, SUB = FALSE, FACT = 1, Vars = c("Time", "r"),
                                            ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                            VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

#for each separately 
mend_unique_results<-mend_function(data=dat, SUB = FALSE, FACT = 3, Vars = c("Time", "r"),
                                        ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                        VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

decay_all_results$likelihood
decay_structures_results$likelihood
decay_substrates_results$likelihood
decay_unique_results$likelihood

monod_all_results$ll_r
monod_structures_results$ll_r
monod_substrates_results$ll_r
monod_unique_results$ll_r

mem_all_results$ll_r
mem_structures_results$ll_r
mem_substrates_results$ll_r
mem_unique_results$ll_r

mend_all_results$ll_r
mend_structures_results$ll_r
mend_substrates_results$ll_r
mend_unique_results$ll_r


###############################################################################################
###########################Metabolic Microbial - enzyme model##################################
###############################################################################################
source("../mmem_function.R")


#across all treatments
mmem_all_results<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Time"),
                                ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)


#no_cors<-detectCores()-1
#cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mmem_structures_results<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Vars = c("Substrate", "r", "Time"),
                                              ColM = c("time", "r", "r1", "r2"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                              VarsCmic = c("Substrate", "Cmic.dna", "Time"), Niter = 10000)

#for different substrates
mmem_substrates_results<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Vars = c("Time", "r"),
                                              ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                              VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

#for each separately 
mmem_unique_results<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Vars = c("Time", "r"),
                                          ColV = c("time", "r"),Cmic=mean(dat$Cmic.dna.init, na.rm = T),
                                          VarsCmic = c("Time", "Cmic.dna"), Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

decay_all_results$likelihood
decay_structures_results$likelihood
decay_substrates_results$likelihood
decay_unique_results$likelihood

monod_all_results$ll_r
monod_structures_results$ll_r
monod_substrates_results$ll_r
monod_unique_results$ll_r

mem_all_results$ll_r
mem_structures_results$ll_r
mem_substrates_results$ll_r
mem_unique_results$ll_r

mend_all_results$ll_r
mend_structures_results$ll_r
mend_substrates_results$ll_r
mend_unique_results$ll_r

mmem_all_results$ll_r
mmem_structures_results$ll_r
mmem_substrates_results$ll_r
mmem_unique_results$ll_r

mmem_all_results_m$ll_r

plot(obs_E~mod_E, mmem_unique_results$OvP_E)
####################################################################################################
####################################################################################################
####################################################################################################
#Monod, MEM, MEND and MMEM are further calibrated against respiration and microbial biomass as well
#since several conversion factors between DNA or protein content and micrbial biomass C exist, models
#are calibrated against 4 different datasets with 4 different conversion factors used:
#highest and lowest DNA based microbial biomass and highest and lowest protein based microbial biomass

#1. Highest DNA
#measured microbial biomass
dat$Cmic.dna_H<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]/0.004*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_H<-mean(m0$DNA.init/0.004*0.45/12.01/3*0.25, na.rm=T)

##################modeling
#monod growth
monod_all_results2<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H","Time"),
                                         ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_H,
                                         VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

#mem
mem_all_results2<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H","Time"),
                              ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_H,
                              VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

#mend
mend_all_results2<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H","Time"),
                                ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_H,
                                VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

#mmem
mmem_all_results2<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H","Time"),
                                ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_H,
                                VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

#results
monod_all_results2$ll_r
mem_all_results2$ll_r
mend_all_results2$ll_r
mmem_all_results2$ll_r

monod_all_results2$ll_Cmic
mem_all_results2$ll_Cmic
mend_all_results2$ll_Cmic
mmem_all_results2$ll_Cmic

#2. Lowest DNA
#measured microbial biomass
dat$Cmic.dna_L<-(m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]*0.2503+15)*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_L<-mean((m0$DNA.init*0.2503+15)*0.45/12.01/3*0.25, na.rm=T)

##################modeling
#monod growth
monod_all_results3<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L","Time"),
                                          ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_L,
                                          VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)

#mem
mem_all_results3<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L","Time"),
                               ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_L,
                               VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)

#mend
mend_all_results3<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_L,
                                 VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)

#mmem
mmem_all_results3<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_L,
                                 VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)

#results
monod_all_results2$ll_r
mem_all_results2$ll_r
mend_all_results2$ll_r
mmem_all_results2$ll_r

monod_all_results2$ll_Cmic
mem_all_results2$ll_Cmic
mend_all_results2$ll_Cmic
mmem_all_results2$ll_Cmic

monod_all_results3$ll_r
mem_all_results3$ll_r
mend_all_results3$ll_r
mmem_all_results3$ll_r

monod_all_results3$ll_Cmic
mem_all_results3$ll_Cmic
mend_all_results3$ll_Cmic
mmem_all_results3$ll_Cmic

#3. Median DNA - this is 2.39% of biomass - the best choice of a modeller
#measured microbial biomass
dat$Cmic.dna_M<-(m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]/0.0239)*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_M<-mean((m0$DNA.init/0.0239)*0.45/12.01/3*0.25, na.rm=T)

##################modeling
#monod growth
monod_all_results4<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M","Time"),
                                          ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_M,
                                          VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)

#mem
mem_all_results4<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M","Time"),
                               ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_M,
                               VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)

#mend
mend_all_results4<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_M,
                                 VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)

#mmem
mmem_all_results4<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.dna.init_M,
                                 VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)

#results
monod_all_results2$ll_r
mem_all_results2$ll_r
mend_all_results2$ll_r
mmem_all_results2$ll_r

monod_all_results2$ll_Cmic
mem_all_results2$ll_Cmic
mend_all_results2$ll_Cmic
mmem_all_results2$ll_Cmic

monod_all_results3$ll_r
mem_all_results3$ll_r
mend_all_results3$ll_r
mmem_all_results3$ll_r

monod_all_results3$ll_Cmic
mem_all_results3$ll_Cmic
mend_all_results3$ll_Cmic
mmem_all_results3$ll_Cmic

monod_all_results4$ll_r
mem_all_results4$ll_r
mend_all_results4$ll_r
mmem_all_results4$ll_r

monod_all_results4$ll_Cmic
mem_all_results4$ll_Cmic
mend_all_results4$ll_Cmic
mmem_all_results4$ll_Cmic

#5. Lowest protein 
#gap filling
ggplot(m0, aes(DNA, Prot.in))+geom_point()+geom_smooth(method=lm)
summary(lm(Prot.in~DNA, m0))
fill_coefs<-coef(lm(Prot.in~DNA, m0))

m0[is.na(m0$Prot.in), "Prot.in"]<-m0[is.na(m0$Prot.in), "DNA"]*fill_coefs[2]+fill_coefs[1]


#measured microbial biomass
dat$Cmic.prot_L<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "Prot.in"]/0.272*0.45/12.01/4
#measured initial microbial biomass
Cmic.prot.init_L<-mean((m0$DNA.init*fill_coefs[2]+fill_coefs[1])/0.272*0.45/12.01/3*0.25, na.rm=T)

##################modeling
#monod growth
monod_all_results5<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L","Time"),
                                          ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_L,
                                          VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

#mem
mem_all_results5<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L","Time"),
                               ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_L,
                               VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

#mend
mend_all_results5<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_L,
                                 VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

#mmem
mmem_all_results5<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_L,
                                 VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

#results
monod_all_results2$ll_r
mem_all_results2$ll_r
mend_all_results2$ll_r
mmem_all_results2$ll_r

monod_all_results2$ll_Cmic
mem_all_results2$ll_Cmic
mend_all_results2$ll_Cmic
mmem_all_results2$ll_Cmic

monod_all_results3$ll_r
mem_all_results3$ll_r
mend_all_results3$ll_r
mmem_all_results3$ll_r

monod_all_results3$ll_Cmic
mem_all_results3$ll_Cmic
mend_all_results3$ll_Cmic
mmem_all_results3$ll_Cmic

monod_all_results4$ll_r
mem_all_results4$ll_r
mend_all_results4$ll_r
mmem_all_results4$ll_r

monod_all_results4$ll_Cmic
mem_all_results4$ll_Cmic
mend_all_results4$ll_Cmic
mmem_all_results4$ll_Cmic

monod_all_results5$ll_r
mem_all_results5$ll_r
mend_all_results5$ll_r
mmem_all_results5$ll_r

monod_all_results5$ll_Cmic
mem_all_results5$ll_Cmic
mend_all_results5$ll_Cmic
mmem_all_results5$ll_Cmic

#6. Highest protein 
#gap filling
#measured microbial biomass
dat$Cmic.prot_H<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "Prot.in"]/0.82*0.45/12.01/4
#measured initial microbial biomass
Cmic.prot.init_H<-mean((m0$DNA.init*fill_coefs[2]+fill_coefs[1])/0.82*0.45/12.01/3*0.25, na.rm=T)

##################modeling
#monod growth
monod_all_results6<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H","Time"),
                                          ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_H,
                                          VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)

#mem
mem_all_results6<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H","Time"),
                               ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_H,
                               VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)

#mend
mend_all_results6<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_H,
                                 VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)

#mmem
mmem_all_results6<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H","Time"),
                                 ColM = c("time", "r", "Cmic","r1", "Cmic1","r2", "Cmic2"),Cmic=Cmic.prot.init_H,
                                 VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)

#results
monod_all_results2$ll_r
mem_all_results2$ll_r
mend_all_results2$ll_r
mmem_all_results2$ll_r

monod_all_results2$ll_Cmic
mem_all_results2$ll_Cmic
mend_all_results2$ll_Cmic
mmem_all_results2$ll_Cmic

monod_all_results3$ll_r
mem_all_results3$ll_r
mend_all_results3$ll_r
mmem_all_results3$ll_r

monod_all_results3$ll_Cmic
mem_all_results3$ll_Cmic
mend_all_results3$ll_Cmic
mmem_all_results3$ll_Cmic

monod_all_results4$ll_r
mem_all_results4$ll_r
mend_all_results4$ll_r
mmem_all_results4$ll_r

monod_all_results4$ll_Cmic
mem_all_results4$ll_Cmic
mend_all_results4$ll_Cmic
mmem_all_results4$ll_Cmic

monod_all_results5$ll_r
mem_all_results5$ll_r
mend_all_results5$ll_r
mmem_all_results5$ll_r

monod_all_results5$ll_Cmic
mem_all_results5$ll_Cmic
mend_all_results5$ll_Cmic
mmem_all_results5$ll_Cmic

monod_all_results6$ll_r
mem_all_results6$ll_r
mend_all_results6$ll_r
mmem_all_results6$ll_r

monod_all_results6$ll_Cmic
mem_all_results6$ll_Cmic
mend_all_results6$ll_Cmic
mmem_all_results6$ll_Cmic


######################################################################################################
############################################Enzymes###################################################
######################################################################################################

#MEM, MEND and MMEM predicts the enzyme concentration as well
#Model calibration is thus further constrained by the enzyme concentration
#I use 4 different microbial biomass data sets

#1. Highest DNA
#measured microbial biomass
dat$Cmic.dna_H<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]/0.004*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_H<-mean(m0$DNA.init/0.004*0.45/12.01/3*0.25, na.rm=T)

##################modeling
#mem
mem_all_results7<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H", "E", "Time"),
                               ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_H,
                               VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

#mend
mend_all_results7<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_H,
                                 VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

#mmem
mmem_all_results7<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_H", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_H,
                                 VarsCmic = c("Substrate", "Cmic.dna_H", "Time"), Niter = 100000)

ggplot(dat, aes(Time, Cmic.dna_H))+geom_point(aes(colour=Structure))

mem_all_results7$ll_r
mend_all_results7$ll_r
mmem_all_results7$ll_r


mem_all_results7$ll_Cmic
mend_all_results7$ll_Cmic
mmem_all_results7$ll_Cmic


mem_all_results7$ll_E
mend_all_results7$ll_E
mmem_all_results7$ll_E

#2. Median DNA
#measured microbial biomass
dat$Cmic.dna_M<-(m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]/0.0239)*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_M<-mean((m0$DNA.init/0.0239)*0.45/12.01/3*0.25, na.rm=T)


##################modeling
#mem
mem_all_results8<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M", "E", "Time"),
                               ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_M,
                               VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)

#mend
mend_all_results8<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_M,
                                 VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)

#mmem
mmem_all_results8<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_M", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_M,
                                 VarsCmic = c("Substrate", "Cmic.dna_M", "Time"), Niter = 100000)


mem_all_results7$ll_r
mend_all_results7$ll_r
mmem_all_results7$ll_r


mem_all_results7$ll_Cmic
mend_all_results7$ll_Cmic
mmem_all_results7$ll_Cmic


mem_all_results7$ll_E
mend_all_results7$ll_E
mmem_all_results7$ll_E


mem_all_results8$ll_r
mend_all_results8$ll_r
mmem_all_results8$ll_r


mem_all_results8$ll_Cmic
mend_all_results8$ll_Cmic
mmem_all_results8$ll_Cmic


mem_all_results8$ll_E
mend_all_results8$ll_E
mmem_all_results8$ll_E

#3. Lowest protein
##################modeling
#mem
mem_all_results9<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L", "E", "Time"),
                               ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.prot.init_L,
                               VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

#mend
mend_all_results9<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.prot.init_L,
                                 VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

#mmem
mmem_all_results9<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_L", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.prot.init_L,
                                 VarsCmic = c("Substrate", "Cmic.prot_L", "Time"), Niter = 100000)

mem_all_results7$ll_r
mend_all_results7$ll_r
mmem_all_results7$ll_r


mem_all_results7$ll_Cmic
mend_all_results7$ll_Cmic
mmem_all_results7$ll_Cmic


mem_all_results7$ll_E
mend_all_results7$ll_E
mmem_all_results7$ll_E


mem_all_results8$ll_r
mend_all_results8$ll_r
mmem_all_results8$ll_r


mem_all_results8$ll_Cmic
mend_all_results8$ll_Cmic
mmem_all_results8$ll_Cmic


mem_all_results8$ll_E
mend_all_results8$ll_E
mmem_all_results8$ll_E


mem_all_results9$ll_r
mend_all_results9$ll_r
mmem_all_results9$ll_r


mem_all_results9$ll_Cmic
mend_all_results9$ll_Cmic
mmem_all_results9$ll_Cmic


mem_all_results9$ll_E
mend_all_results9$ll_E
mmem_all_results9$ll_E

#4. Highest protein
##################modeling
#mem
mem_all_results10<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H", "E", "Time"),
                               ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.prot.init_H,
                               VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)

#mend
mend_all_results10<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.prot.init_H,
                                 VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)

#mmem
mmem_all_results10<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.prot_H", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.prot.init_H,
                                 VarsCmic = c("Substrate", "Cmic.prot_H", "Time"), Niter = 100000)


mem_all_results7$ll_r
mend_all_results7$ll_r
mmem_all_results7$ll_r


mem_all_results7$ll_Cmic
mend_all_results7$ll_Cmic
mmem_all_results7$ll_Cmic


mem_all_results7$ll_E
mend_all_results7$ll_E
mmem_all_results7$ll_E


mem_all_results8$ll_r
mend_all_results8$ll_r
mmem_all_results8$ll_r


mem_all_results8$ll_Cmic
mend_all_results8$ll_Cmic
mmem_all_results8$ll_Cmic


mem_all_results8$ll_E
mend_all_results8$ll_E
mmem_all_results8$ll_E


mem_all_results9$ll_r
mend_all_results9$ll_r
mmem_all_results9$ll_r


mem_all_results9$ll_Cmic
mend_all_results9$ll_Cmic
mmem_all_results9$ll_Cmic


mem_all_results9$ll_E
mend_all_results9$ll_E
mmem_all_results9$ll_E

mem_all_results10$ll_r
mend_all_results10$ll_r
mmem_all_results10$ll_r


mem_all_results10$ll_Cmic
mend_all_results10$ll_Cmic
mmem_all_results10$ll_Cmic


mem_all_results10$ll_E
mend_all_results10$ll_E
mmem_all_results10$ll_E

mem_all_results11$ll_r
mend_all_results11$ll_r
mmem_all_results11$ll_r


mem_all_results11$ll_Cmic
mend_all_results11$ll_Cmic
mmem_all_results11$ll_Cmic


mem_all_results11$ll_E
mend_all_results11$ll_E
mmem_all_results11$ll_E

#Lowest DNA

##################modeling
#mem
mem_all_results11<-mem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L", "E", "Time"),
                               ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_L,
                               VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)

#mend
mend_all_results11<-mend_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_L,
                                 VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)

#mmem
mmem_all_results11<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Vars = c("Substrate", "r", "Cmic.dna_L", "E", "Time"),
                                 ColM = c("time", "r", "Cmic","E","r1", "Cmic1","E1","r2", "Cmic2", "E2"),Cmic=Cmic.dna.init_L,
                                 VarsCmic = c("Substrate", "Cmic.dna_L", "Time"), Niter = 100000)


