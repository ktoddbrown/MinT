#################################################MinT#################################################
##libraries
library(deSolve)
#library(dplyr)
library(FME)
library(reshape)
#library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
#library(data.table)
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

#extracellular protein to extracellular protein carbon (46% of carbon in protein - Vrede et al., 2004)
m0$E<-m0$Prot.out*0.46/12.01/4


#extracting only the columns of interest
mint<-m0[,c("Structure", "Substrate", "r", "Time", "DOCinit", "DNA", "Prot.in", "E")]
mint<-mint[!is.na(mint$Substrate), ]


dat<-subset(mint, Substrate!="Free" & Substrate!="Celluloze")
summary(dat)

#Four datasets of microbial biomass are defined 

#1. Highest DNA
#measured microbial biomass
dat$Cmic.dna_H<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]/0.004*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_H<-mean(m0$DNA.init/0.004*0.45/12.01/3*0.25, na.rm=T)

#2. Median DNA - this is 2.39% of biomass - the best choice of a modeller
#measured microbial biomass
dat$Cmic.dna_M<-(m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "DNA"]/0.0239)*0.45/12.01/4
#measured initial microbial biomass
Cmic.dna.init_M<-mean((m0$DNA.init/0.0239)*0.45/12.01/3*0.25, na.rm=T)

#3. Lowest protein 
#gap filling
ggplot(m0, aes(DNA, Prot.in))+geom_point()+geom_smooth(method=lm)
summary(lm(Prot.in~DNA, m0))
fill_coefs<-coef(lm(Prot.in~DNA, m0))

m0[is.na(m0$Prot.in), "Prot.in"]<-m0[is.na(m0$Prot.in), "DNA"]*fill_coefs[2]+fill_coefs[1]


#measured microbial biomass
dat$Cmic.prot_L<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "Prot.in"]/0.272*0.45/12.01/4
#measured initial microbial biomass
Cmic.prot.init_L<-mean((m0$DNA.init*fill_coefs[2]+fill_coefs[1])/0.272*0.45/12.01/3*0.25, na.rm=T)



#4. Highest protein 
#gap filling
#measured microbial biomass
dat$Cmic.prot_H<-m0[(m0$Substrate!="Free" & m0$Substrate!="Celluloze" & !is.na(m0$Substrate)), "Prot.in"]/0.82*0.45/12.01/4
#measured initial microbial biomass
Cmic.prot.init_H<-mean((m0$DNA.init*fill_coefs[2]+fill_coefs[1])/0.82*0.45/12.01/3*0.25, na.rm=T)



#5. DEB require protein and DNA C concentration
dat$DNAc<-dat$DNA*0.51/12.01/4
dat$Protc<-dat$Prot.in*0.46/12.01/4

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
#a stands for highest DNA, b for median, c for highest protein and d for lowest protein
monod_all_a<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
monod_all_b<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean = TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
monod_all_c<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
monod_all_d<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod_structures_a<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
monod_structures_b<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
monod_structures_c<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
monod_structures_d<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

#for different substrates
monod_substrates_a<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                                Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                                ColM = c("time", "r"),
                                                Cmic=Cmic.dna.init_H,
                                                Niter = 10000)
monod_substrates_b<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                                 ColM = c("time", "r"),
                                                 Cmic=Cmic.dna.init_M,
                                                 Niter = 10000)
monod_substrates_c<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
monod_substrates_d<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

#for each separately 
monod_each_a<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
monod_each_b<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
monod_each_c<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
monod_each_d<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

decay_all_results$likelihood
decay_structures_results$likelihood
decay_substrates_results$likelihood
decay_unique_results$likelihood

monod_all_a$ll
monod_all_b$ll
monod_all_c$ll
monod_all_d$ll

monod_substrates_a$ll
monod_substrates_b$ll
monod_substrates_c$ll
monod_substrates_d$ll

monod_structures_a$ll
monod_structures_b$ll
monod_structures_c$ll
monod_structures_d$ll

monod_each_a$ll
monod_each_b$ll
monod_each_c$ll
monod_each_d$ll

###############################################################################################
#####################################Microbial - enzyme model##################################
###############################################################################################
source("../mem_function.R")


#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mem_all_a<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                   Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                   ColM = c("time", "r", "r1", "r2"),
                                   Cmic=Cmic.dna.init_H,
                                   Niter = 10000)
mem_all_b<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                   Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                   ColM = c("time", "r", "r1", "r2"),
                                   Cmic=Cmic.dna.init_M,
                                   Niter = 10000)
mem_all_c<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                   Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                   ColM = c("time", "r", "r1", "r2"),
                                   Cmic=Cmic.prot.init_H,
                                   Niter = 10000)
mem_all_d<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                                   Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                   ColM = c("time", "r", "r1", "r2"),
                                   Cmic=Cmic.prot.init_L,
                                   Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mem_structures_a<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
mem_structures_b<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
mem_structures_c<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
mem_structures_d<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "r", "r1", "r2"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

#for different substrates
mem_substrates_a<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
mem_substrates_b<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
mem_substrates_c<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
mem_substrates_d<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "r"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

#for each separately 
mem_each_a<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "r"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
mem_each_b<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "r"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
mem_each_c<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "r"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
mem_each_d<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "r"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mem_all_a$ll
mem_all_b$ll
mem_all_c$ll
mem_all_d$ll

mem_substrates_a$ll
mem_substrates_b$ll
mem_substrates_c$ll
mem_substrates_d$ll

mem_structures_a$ll
mem_structures_b$ll
mem_structures_c$ll
mem_structures_d$ll

mem_each_a$ll
mem_each_b$ll
mem_each_c$ll
mem_each_d$ll

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
source("../mend_function.R")

#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mend_all_a<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                        Vars = c("Time", "r", "E", "Cmic.dna_H"),
                        ColM = c("time", "r", "r1", "r2"),
                        Cmic=Cmic.dna.init_H,
                        Niter = 10000)
mend_all_b<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                        Vars = c("Time", "r", "E", "Cmic.dna_M"),
                        ColM = c("time", "r", "r1", "r2"),
                        Cmic=Cmic.dna.init_M,
                        Niter = 10000)
mend_all_c<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                        Vars = c("Time", "r", "E", "Cmic.prot_H"),
                        ColM = c("time", "r", "r1", "r2"),
                        Cmic=Cmic.prot.init_H,
                        Niter = 10000)
mend_all_d<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                        Vars = c("Time", "r", "E", "Cmic.prot_L"),
                        ColM = c("time", "r", "r1", "r2"),
                        Cmic=Cmic.prot.init_L,
                        Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_structures_a<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.dna_H"),
                               ColM = c("time", "r", "r1", "r2"),
                               Cmic=Cmic.dna.init_H,
                               Niter = 10000)
mend_structures_b<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.dna_M"),
                               ColM = c("time", "r", "r1", "r2"),
                               Cmic=Cmic.dna.init_M,
                               Niter = 10000)
mend_structures_c<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.prot_H"),
                               ColM = c("time", "r", "r1", "r2"),
                               Cmic=Cmic.prot.init_H,
                               Niter = 10000)
mend_structures_d<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.prot_L"),
                               ColM = c("time", "r", "r1", "r2"),
                               Cmic=Cmic.prot.init_L,
                               Niter = 10000)

#for different substrates
mend_substrates_a<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.dna_H"),
                               ColM = c("time", "r"),
                               Cmic=Cmic.dna.init_H,
                               Niter = 10000)
mend_substrates_b<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.dna_M"),
                               ColM = c("time", "r"),
                               Cmic=Cmic.dna.init_M,
                               Niter = 10000)
mend_substrates_c<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.prot_H"),
                               ColM = c("time", "r"),
                               Cmic=Cmic.prot.init_H,
                               Niter = 10000)
mend_substrates_d<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                               Vars = c("Time", "r", "E", "Cmic.prot_L"),
                               ColM = c("time", "r"),
                               Cmic=Cmic.prot.init_L,
                               Niter = 10000)

#for each separately 
mend_each_a<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                         Vars = c("Time", "r", "E", "Cmic.dna_H"),
                         ColM = c("time", "r"),
                         Cmic=Cmic.dna.init_H,
                         Niter = 10000)
mend_each_b<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                         Vars = c("Time", "r", "E", "Cmic.dna_M"),
                         ColM = c("time", "r"),
                         Cmic=Cmic.dna.init_M,
                         Niter = 10000)
mend_each_c<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                         Vars = c("Time", "r", "E", "Cmic.prot_H"),
                         ColM = c("time", "r"),
                         Cmic=Cmic.prot.init_H,
                         Niter = 10000)
mend_each_d<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                         Vars = c("Time", "r", "E", "Cmic.prot_L"),
                         ColM = c("time", "r"),
                         Cmic=Cmic.prot.init_L,
                         Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mend_all_a$ll
mend_all_b$ll
mend_all_c$ll
mend_all_d$ll

mend_substrates_a$ll
mend_substrates_b$ll
mend_substrates_c$ll
mend_substrates_d$ll

mend_structures_a$ll
mend_structures_b$ll
mend_structures_c$ll
mend_structures_d$ll

mend_each_a$ll
mend_each_b$ll
mend_each_c$ll
mend_each_d$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model##################################
###############################################################################################
source("../mmem_function.R")


#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mmem_all_a<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.dna.init_H,
                          Niter = 10000)
mmem_all_b<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.dna.init_M,
                          Niter = 10000)
mmem_all_c<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.prot.init_H,
                          Niter = 10000)
mmem_all_d<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.prot.init_L,
                          Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mmem_structures_a<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
mmem_structures_b<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
mmem_structures_c<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
mmem_structures_d<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for different substrates
mmem_substrates_a<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
mmem_substrates_b<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
mmem_substrates_c<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
mmem_substrates_d<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for each separately 
mmem_each_a<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
mmem_each_b<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
mmem_each_c<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
mmem_each_d<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mmem_all_a$ll
mmem_all_b$ll
mmem_all_c$ll
mmem_all_d$ll

mmem_substrates_a$ll
mmem_substrates_b$ll
mmem_substrates_c$ll
mmem_substrates_d$ll

mmem_structures_a$ll
mmem_structures_b$ll
mmem_structures_c$ll
mmem_structures_d$ll

mmem_each_a$ll
mmem_each_b$ll
mmem_each_c$ll
mmem_each_d$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model by Andrew##################################
###############################################################################################
source("../andrew_function.R")


#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
an_all_a<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.dna.init_H,
                          Niter = 10000)
an_all_b<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.dna.init_M,
                          Niter = 10000)
an_all_c<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.prot.init_H,
                          Niter = 10000)
an_all_d<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                          ColM = c("time", "r", "r1", "r2"),
                          Cmic=Cmic.prot.init_L,
                          Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
an_structures_a<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
an_structures_b<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
an_structures_c<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
an_structures_d<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "r", "r1", "r2"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for different substrates
an_substrates_a<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
an_substrates_b<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
an_substrates_c<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
an_substrates_d<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "r"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for each separately 
an_each_a<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
an_each_b<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
an_each_c<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
an_each_d<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time", "r"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

an_all_a$ll
an_all_b$ll
an_all_c$ll
an_all_d$ll

an_substrates_a$ll
an_substrates_b$ll
an_substrates_c$ll
an_substrates_d$ll

an_structures_a$ll
an_structures_b$ll
an_structures_c$ll
an_structures_d$ll

an_each_a$ll
an_each_b$ll
an_each_c$ll
an_each_d$ll

####################################################################################################
####################################################################################################
####################################################################################################
#All models are further calibrated against the microbial biomass because models calibrated against  
#the respiration rate do not predict microbial biomass well
###############################################################################################
#####################################Monod growth##############################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and d for lowest protein
monod_all_mica<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                                   Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                   ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                   Cmic=Cmic.dna.init_H,
                                   Niter = 10000)
monod_all_micb<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                                   Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                   ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                   Cmic=Cmic.dna.init_M,
                                   Niter = 10000)
monod_all_micc<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                                   Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                   ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                   Cmic=Cmic.prot.init_H,
                                   Niter = 10000)
monod_all_micd<-monod_growth_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                                   Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                   ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                   Cmic=Cmic.prot.init_L,
                                   Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
monod_structures_mica<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
monod_structures_micb<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
monod_structures_micc<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
monod_structures_micd<-monod_growth_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

#for different substrates
monod_substrates_mica<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                          ColM = c("time", "Cmic"),
                                          Cmic=Cmic.dna.init_H,
                                          Niter = 10000)
monod_substrates_micb<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                          ColM = c("time", "Cmic"),
                                          Cmic=Cmic.dna.init_M,
                                          Niter = 10000)
monod_substrates_micc<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                          ColM = c("time", "Cmic"),
                                          Cmic=Cmic.prot.init_H,
                                          Niter = 10000)
monod_substrates_micd<-monod_growth_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                          ColM = c("time", "Cmic"),
                                          Cmic=Cmic.prot.init_L,
                                          Niter = 10000)

#for each separately 
monod_each_mica<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "Cmic"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
monod_each_micb<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "Cmic"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
monod_each_micc<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "Cmic"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
monod_each_micd<-monod_growth_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "Cmic"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

monod_all_mica$ll
monod_all_micb$ll
monod_all_micc$ll
monod_all_micd$ll

monod_substrates_mica$ll
monod_substrates_micb$ll
monod_substrates_micc$ll
monod_substrates_micd$ll

monod_structures_mica$ll
monod_structures_micb$ll
monod_structures_micc$ll
monod_structures_micd$ll

monod_each_mica$ll
monod_each_micb$ll
monod_each_micc$ll
monod_each_micd$ll

###############################################################################################
#####################################Microbial - enzyme model##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mem_all_mica<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                        Vars = c("Time", "r", "E", "Cmic.dna_H"),
                        ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                        Cmic=Cmic.dna.init_H,
                        Niter = 10000)
mem_all_micb<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                        Vars = c("Time", "r", "E", "Cmic.dna_M"),
                        ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                        Cmic=Cmic.dna.init_M,
                        Niter = 10000)
mem_all_micc<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                        Vars = c("Time", "r", "E", "Cmic.prot_H"),
                        ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                        Cmic=Cmic.prot.init_H,
                        Niter = 10000)
mem_all_micd<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                        Vars = c("Time", "r", "E", "Cmic.prot_L"),
                        ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                        Cmic=Cmic.prot.init_L,
                        Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mem_structures_mica<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.dna_H"),
                               ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                               Cmic=Cmic.dna.init_H,
                               Niter = 10000)
mem_structures_micb<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.dna_M"),
                               ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                               Cmic=Cmic.dna.init_M,
                               Niter = 10000)
mem_structures_micc<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.prot_H"),
                               ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                               Cmic=Cmic.prot.init_H,
                               Niter = 10000)
mem_structures_micd<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.prot_L"),
                               ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                               Cmic=Cmic.prot.init_L,
                               Niter = 10000)

#for different substrates
mem_substrates_mica<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.dna_H"),
                               ColM = c("time", "Cmic"),
                               Cmic=Cmic.dna.init_H,
                               Niter = 10000)
mem_substrates_micb<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.dna_M"),
                               ColM = c("time", "Cmic"),
                               Cmic=Cmic.dna.init_M,
                               Niter = 10000)
mem_substrates_micc<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.prot_H"),
                               ColM = c("time", "Cmic"),
                               Cmic=Cmic.prot.init_H,
                               Niter = 10000)
mem_substrates_micd<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                               Vars = c("Time", "r", "E", "Cmic.prot_L"),
                               ColM = c("time", "Cmic"),
                               Cmic=Cmic.prot.init_L,
                               Niter = 10000)

#for each separately 
mem_each_mica<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.dna_H"),
                         ColM = c("time", "Cmic"),
                         Cmic=Cmic.dna.init_H,
                         Niter = 10000)
mem_each_micb<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.dna_M"),
                         ColM = c("time", "Cmic"),
                         Cmic=Cmic.dna.init_M,
                         Niter = 10000)
mem_each_micc<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.prot_H"),
                         ColM = c("time", "Cmic"),
                         Cmic=Cmic.prot.init_H,
                         Niter = 10000)
mem_each_micd<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.prot_L"),
                         ColM = c("time", "Cmic"),
                         Cmic=Cmic.prot.init_L,
                         Niter = 10000)

#stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mem_all_mica$ll
mem_all_micb$ll
mem_all_micc$ll
mem_all_micd$ll

mem_substrates_mica$ll
mem_substrates_micb$ll
mem_substrates_micc$ll
mem_substrates_micd$ll

mem_structures_mica$ll
mem_structures_micb$ll
mem_structures_micc$ll
mem_structures_micd$ll

mem_each_mica$ll
mem_each_micb$ll
mem_each_micc$ll
mem_each_micd$ll

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mend_all_mica<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.dna.init_H,
                          Niter = 10000)
mend_all_micb<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.dna.init_M,
                          Niter = 10000)
mend_all_micc<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.prot.init_H,
                          Niter = 10000)
mend_all_micd<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.prot.init_L,
                          Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_structures_mica<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
mend_structures_micb<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
mend_structures_micc<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
mend_structures_micd<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for different substrates
mend_substrates_mica<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
mend_substrates_micb<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
mend_substrates_micc<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
mend_substrates_micd<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for each separately 
mend_each_mica<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
mend_each_micb<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
mend_each_micc<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
mend_each_micd<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)

stopImplicitCluster()-1

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mend_all_mica$ll
mend_all_micb$ll
mend_all_micc$ll
mend_all_micd$ll

mend_substrates_mica$ll
mend_substrates_micb$ll
mend_substrates_micc$ll
mend_substrates_micd$ll

mend_structures_mica$ll
mend_structures_micb$ll
mend_structures_micc$ll
mend_structures_micd$ll

mend_each_mica$ll
mend_each_micb$ll
mend_each_micc$ll
mend_each_micd$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mmem_all_mica<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.dna.init_H,
                          Niter = 10000)
mmem_all_micb<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.dna.init_M,
                          Niter = 10000)
mmem_all_micc<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.prot.init_H,
                          Niter = 10000)
mmem_all_micd<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.prot.init_L,
                          Niter = 10000)



no_cors<-detectCores()
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mmem_structures_mica<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
mmem_structures_micb<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
mmem_structures_micc<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
mmem_structures_micd<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for different substrates
mmem_substrates_mica<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
mmem_substrates_micb<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
mmem_substrates_micc<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
mmem_substrates_micd<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for each separately 
mmem_each_mica<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
mmem_each_micb<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
mmem_each_micc<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
mmem_each_micd<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mmem_all_mica$ll
mmem_all_micb$ll
mmem_all_micc$ll
mmem_all_micd$ll

mmem_substrates_mica$ll
mmem_substrates_micb$ll
mmem_substrates_micc$ll
mmem_substrates_micd$ll

mmem_structures_mica$ll
mmem_structures_micb$ll
mmem_structures_micc$ll
mmem_structures_micd$ll

mmem_each_mica$ll
mmem_each_micb$ll
mmem_each_micc$ll
mmem_each_micd$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model by Andrew##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
an_all_mica<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.dna.init_H,
                          Niter = 10000)
an_all_micb<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.dna.init_M,
                          Niter = 10000)
an_all_micc<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.prot.init_H,
                          Niter = 10000)
an_all_micd<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                          ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                          Cmic=Cmic.prot.init_L,
                          Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
an_structures_mica<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
an_structures_micb<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
an_structures_micc<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
an_structures_micd<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "Cmic", "Cmic1", "Cmic2"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for different substrates
an_substrates_mica<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.dna.init_H,
                                 Niter = 10000)
an_substrates_micb<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.dna.init_M,
                                 Niter = 10000)
an_substrates_micc<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.prot.init_H,
                                 Niter = 10000)
an_substrates_micd<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                 Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                 ColM = c("time", "Cmic"),
                                 Cmic=Cmic.prot.init_L,
                                 Niter = 10000)

#for each separately 
an_each_mica<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
an_each_micb<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
an_each_micc<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
an_each_micd<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time", "Cmic"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

an_all_mica$ll
an_all_micb$ll
an_all_micc$ll
an_all_micd$ll

an_substrates_mica$ll
an_substrates_micb$ll
an_substrates_micc$ll
an_substrates_micd$ll

an_structures_mica$ll
an_structures_micb$ll
an_structures_micc$ll
an_structures_micd$ll

an_each_mica$ll
an_each_micb$ll
an_each_micc$ll
an_each_micd$ll

####################################################################################################
####################################################################################################
####################################################################################################
#All models are further calibrated against the enzymes

###############################################################################################
#####################################Microbial - enzyme model##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mem_all_Ea<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time", "E", "E1", "E2"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
mem_all_Eb<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time", "E", "E1", "E2"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
mem_all_Ec<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time", "E", "E1", "E2"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
mem_all_Ed<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time", "E", "E1", "E2"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mem_structures_Ea<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time", "E", "E1", "E2"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
mem_structures_Eb<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time", "E", "E1", "E2"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
mem_structures_Ec<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time", "E", "E1", "E2"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
mem_structures_Ed<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time", "E", "E1", "E2"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for different substrates
mem_substrates_Ea<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time", "E"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
mem_substrates_Eb<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time", "E"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
mem_substrates_Ec<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time", "E"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
mem_substrates_Ed<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time", "E"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for each separately 
mem_each_Ea<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                            Vars = c("Time", "r", "E", "Cmic.dna_H"),
                            ColM = c("time", "E"),
                            Cmic=Cmic.dna.init_H,
                            Niter = 10000)
mem_each_Eb<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                            Vars = c("Time", "r", "E", "Cmic.dna_M"),
                            ColM = c("time", "E"),
                            Cmic=Cmic.dna.init_M,
                            Niter = 10000)
mem_each_Ec<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                            Vars = c("Time", "r", "E", "Cmic.prot_H"),
                            ColM = c("time", "E"),
                            Cmic=Cmic.prot.init_H,
                            Niter = 10000)
mem_each_Ed<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                            Vars = c("Time", "r", "E", "Cmic.prot_L"),
                            ColM = c("time", "E"),
                            Cmic=Cmic.prot.init_L,
                            Niter = 10000)

#stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mem_all_Ea$ll
mem_all_Eb$ll
mem_all_Ec$ll
mem_all_Ed$ll

mem_substrates_Ea$ll
mem_substrates_Eb$ll
mem_substrates_Ec$ll
mem_substrates_Ed$ll

mem_structures_Ea$ll
mem_structures_Eb$ll
mem_structures_Ec$ll
mem_structures_Ed$ll

mem_each_Ea$ll
mem_each_Eb$ll
mem_each_Ec$ll
mem_each_Ed$ll

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mend_all_Ea<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.dna_H"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.dna.init_H,
                             Niter = 10000)
mend_all_Eb<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.dna_M"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.dna.init_M,
                             Niter = 10000)
mend_all_Ec<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.prot_H"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.prot.init_H,
                             Niter = 10000)
mend_all_Ed<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.prot_L"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.prot.init_L,
                             Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_structures_Ea<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
mend_structures_Eb<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
mend_structures_Ec<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
mend_structures_Ed<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

#for different substrates
mend_substrates_Ea<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
mend_substrates_Eb<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
mend_substrates_Ec<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
mend_substrates_Ed<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

#for each separately 
mend_each_Ea<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.dna_H"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.dna.init_H,
                              Niter = 10000)
mend_each_Eb<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.dna_M"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.dna.init_M,
                              Niter = 10000)
mend_each_Ec<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.prot_H"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.prot.init_H,
                              Niter = 10000)
mend_each_Ed<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.prot_L"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.prot.init_L,
                              Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mend_all_Ea$ll
mend_all_Eb$ll
mend_all_Ec$ll
mend_all_Ed$ll

mend_substrates_Ea$ll
mend_substrates_Eb$ll
mend_substrates_Ec$ll
mend_substrates_Ed$ll

mend_structures_Ea$ll
mend_structures_Eb$ll
mend_structures_Ec$ll
mend_structures_Ed$ll

mend_each_Ea$ll
mend_each_Eb$ll
mend_each_Ec$ll
mend_each_Ed$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mmem_all_Ea<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.dna_H"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.dna.init_H,
                             Niter = 10000)
mmem_all_Eb<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.dna_M"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.dna.init_M,
                             Niter = 10000)
mmem_all_Ec<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.prot_H"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.prot.init_H,
                             Niter = 10000)
mmem_all_Ed<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.prot_L"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.prot.init_L,
                             Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mmem_structures_Ea<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
mmem_structures_Eb<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
mmem_structures_Ec<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
mmem_structures_Ed<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

#for different substrates
mmem_substrates_Ea<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
mmem_substrates_Eb<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
mmem_substrates_Ec<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
mmem_substrates_Ed<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

#for each separately 
mmem_each_Ea<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.dna_H"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.dna.init_H,
                              Niter = 10000)
mmem_each_Eb<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.dna_M"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.dna.init_M,
                              Niter = 10000)
mmem_each_Ec<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.prot_H"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.prot.init_H,
                              Niter = 10000)
mmem_each_Ed<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.prot_L"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.prot.init_L,
                              Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mmem_all_Ea$ll
mmem_all_Eb$ll
mmem_all_Ec$ll
mmem_all_Ed$ll

mmem_substrates_Ea$ll
mmem_substrates_Eb$ll
mmem_substrates_Ec$ll
mmem_substrates_Ed$ll

mmem_structures_Ea$ll
mmem_structures_Eb$ll
mmem_structures_Ec$ll
mmem_structures_Ed$ll

mmem_each_Ea$ll
mmem_each_Eb$ll
mmem_each_Ec$ll
mmem_each_Ed$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model by Andrew##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
an_all_Ea<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.dna_H"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.dna.init_H,
                             Niter = 10000)
an_all_Eb<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.dna_M"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.dna.init_M,
                             Niter = 10000)
an_all_Ec<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.prot_H"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.prot.init_H,
                             Niter = 10000)
an_all_Ed<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                             Vars = c("Time", "r", "E", "Cmic.prot_L"),
                             ColM = c("time", "E", "E1", "E2"),
                             Cmic=Cmic.prot.init_L,
                             Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
an_structures_Ea<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
an_structures_Eb<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
an_structures_Ec<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
an_structures_Ed<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "E", "E1", "E2"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

#for different substrates
an_substrates_Ea<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.dna.init_H,
                                    Niter = 10000)
an_substrates_Eb<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.dna.init_M,
                                    Niter = 10000)
an_substrates_Ec<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.prot.init_H,
                                    Niter = 10000)
an_substrates_Ed<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                    Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                    ColM = c("time", "E"),
                                    Cmic=Cmic.prot.init_L,
                                    Niter = 10000)

#for each separately 
an_each_Ea<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.dna_H"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.dna.init_H,
                              Niter = 10000)
an_each_Eb<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.dna_M"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.dna.init_M,
                              Niter = 10000)
an_each_Ec<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.prot_H"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.prot.init_H,
                              Niter = 10000)
an_each_Ed<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                              Vars = c("Time", "r", "E", "Cmic.prot_L"),
                              ColM = c("time", "E"),
                              Cmic=Cmic.prot.init_L,
                              Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

an_all_Ea$ll
an_all_Eb$ll
an_all_Ec$ll
an_all_Ed$ll

an_substrates_Ea$ll
an_substrates_Eb$ll
an_substrates_Ec$ll
an_substrates_Ed$ll

an_structures_Ea$ll
an_structures_Eb$ll
an_structures_Ec$ll
an_structures_Ed$ll

an_each_Ea$ll
an_each_Eb$ll
an_each_Ec$ll
an_each_Ed$ll

####################################################################################################
####################################################################################################
####################################################################################################
#All models are further calibrated against everything

###############################################################################################
#####################################Microbial - enzyme model##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mem_all_alla<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.dna_H"),
                         ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                         Cmic=Cmic.dna.init_H,
                         Niter = 10000)
mem_all_allb<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.dna_M"),
                         ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                         Cmic=Cmic.dna.init_M,
                         Niter = 10000)
mem_all_allc<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.prot_H"),
                         ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                         Cmic=Cmic.prot.init_H,
                         Niter = 10000)
mem_all_alld<-mem_function(data=dat, SUB = TRUE, FACT = 4, Mean=FALSE,
                         Vars = c("Time", "r", "E", "Cmic.prot_L"),
                         ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                         Cmic=Cmic.prot.init_L,
                         Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mem_structures_alla<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                Cmic=Cmic.dna.init_H,
                                Niter = 10000)
mem_structures_allb<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                Cmic=Cmic.dna.init_M,
                                Niter = 10000)
mem_structures_allc<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                Cmic=Cmic.prot.init_H,
                                Niter = 10000)
mem_structures_alld<-mem_function(data=dat, SUB = TRUE, FACT = 2, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                Cmic=Cmic.prot.init_L,
                                Niter = 10000)

#for different substrates
mem_substrates_alla<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                ColM = c("time", "r","E", "Cmic"),
                                Cmic=Cmic.dna.init_H,
                                Niter = 10000)
mem_substrates_allb<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                ColM = c("time", "r","E", "Cmic"),
                                Cmic=Cmic.dna.init_M,
                                Niter = 10000)
mem_substrates_allc<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                ColM = c("time", "r","E", "Cmic"),
                                Cmic=Cmic.prot.init_H,
                                Niter = 10000)
mem_substrates_alld<-mem_function(data=dat, SUB = FALSE, FACT = 1, Mean=FALSE,
                                Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                ColM = c("time", "r","E", "Cmic"),
                                Cmic=Cmic.prot.init_L,
                                Niter = 10000)

#for each separately 
mem_each_alla<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_H"),
                          ColM = c("time", "E"),
                          Cmic=Cmic.dna.init_H,
                          Niter = 10000)
mem_each_allb<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.dna_M"),
                          ColM = c("time", "r","E", "Cmic"),
                          Cmic=Cmic.dna.init_M,
                          Niter = 10000)
mem_each_allc<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_H"),
                          ColM = c("time", "r","E", "Cmic"),
                          Cmic=Cmic.prot.init_H,
                          Niter = 10000)
mem_each_alld<-mem_function(data=dat, SUB = FALSE, FACT = 3, Mean=FALSE,
                          Vars = c("Time", "r", "E", "Cmic.prot_L"),
                          ColM = c("time", "r","E", "Cmic"),
                          Cmic=Cmic.prot.init_L,
                          Niter = 10000)

#stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mem_all_alla$ll
mem_all_allb$ll
mem_all_allc$ll
mem_all_alld$ll

mem_substrates_alla$ll
mem_substrates_allb$ll
mem_substrates_allc$ll
mem_substrates_alld$ll

mem_structures_alla$ll
mem_structures_allb$ll
mem_structures_allc$ll
mem_structures_alld$ll

mem_each_alla$ll
mem_each_allb$ll
mem_each_allc$ll
mem_each_alld$ll

###############################################################################################
###########################################MEND model##########################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mend_all_alla<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
mend_all_allb<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
mend_all_allc<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
mend_all_alld<-mend_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mend_structures_alla<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
mend_structures_allb<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
mend_structures_allc<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
mend_structures_alld<-mend_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for different substrates
mend_substrates_alla<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
mend_substrates_allb<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
mend_substrates_allc<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
mend_substrates_alld<-mend_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for each separately 
mend_each_alla<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.dna_H"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.dna.init_H,
                            Niter = 10000)
mend_each_allb<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.dna_M"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.dna.init_M,
                            Niter = 10000)
mend_each_allc<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.prot_H"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.prot.init_H,
                            Niter = 10000)
mend_each_alld<-mend_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.prot_L"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.prot.init_L,
                            Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mend_all_alla$ll
mend_all_allb$ll
mend_all_allc$ll
mend_all_alld$ll

mend_substrates_alla$ll
mend_substrates_allb$ll
mend_substrates_allc$ll
mend_substrates_alld$ll

mend_structures_alla$ll
mend_structures_allb$ll
mend_structures_allc$ll
mend_structures_alld$ll

mend_each_alla$ll
mend_each_allb$ll
mend_each_allc$ll
mend_each_alld$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
mmem_all_alla<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
mmem_all_allb<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
mmem_all_allc<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
mmem_all_alld<-mmem_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
mmem_structures_alla<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
mmem_structures_allb<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
mmem_structures_allc<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
mmem_structures_alld<-mmem_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for different substrates
mmem_substrates_alla<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
mmem_substrates_allb<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
mmem_substrates_allc<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
mmem_substrates_alld<-mmem_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for each separately 
mmem_each_alla<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.dna_H"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.dna.init_H,
                            Niter = 10000)
mmem_each_allb<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.dna_M"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.dna.init_M,
                            Niter = 10000)
mmem_each_allc<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.prot_H"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.prot.init_H,
                            Niter = 10000)
mmem_each_alld<-mmem_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.prot_L"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.prot.init_L,
                            Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

mmem_all_alla$ll
mmem_all_allb$ll
mmem_all_allc$ll
mmem_all_alld$ll

mmem_substrates_alla$ll
mmem_substrates_allb$ll
mmem_substrates_allc$ll
mmem_substrates_alld$ll

mmem_structures_alla$ll
mmem_structures_allb$ll
mmem_structures_allc$ll
mmem_structures_alld$ll

mmem_each_alla$ll
mmem_each_allb$ll
mmem_each_allc$ll
mmem_each_alld$ll

###############################################################################################
###########################Metabolic Microbial - enzyme model by Andrew##################################
###############################################################################################
#across all treatments
#a stands for highest DNA, b for median, c for highest protein and c for lowest protein
an_all_alla<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_H"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.dna.init_H,
                           Niter = 10000)
an_all_allb<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.dna_M"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.dna.init_M,
                           Niter = 10000)
an_all_allc<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_H"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.prot.init_H,
                           Niter = 10000)
an_all_alld<-andrew_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE,
                           Vars = c("Time", "r", "E", "Cmic.prot_L"),
                           ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                           Cmic=Cmic.prot.init_L,
                           Niter = 10000)



no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
an_structures_alla<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
an_structures_allb<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
an_structures_allc<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
an_structures_alld<-andrew_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time","r", "r1", "r2", "E", "E1", "E2", "Cmic", "Cmic1", "Cmic2"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for different substrates
an_substrates_alla<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_H"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.dna.init_H,
                                  Niter = 10000)
an_substrates_allb<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.dna_M"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.dna.init_M,
                                  Niter = 10000)
an_substrates_allc<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_H"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.prot.init_H,
                                  Niter = 10000)
an_substrates_alld<-andrew_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE,
                                  Vars = c("Time", "r", "E", "Cmic.prot_L"),
                                  ColM = c("time", "r","E", "Cmic"),
                                  Cmic=Cmic.prot.init_L,
                                  Niter = 10000)

#for each separately 
an_each_alla<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.dna_H"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.dna.init_H,
                            Niter = 10000)
an_each_allb<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.dna_M"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.dna.init_M,
                            Niter = 10000)
an_each_allc<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.prot_H"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.prot.init_H,
                            Niter = 10000)
an_each_alld<-andrew_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE,
                            Vars = c("Time", "r", "E", "Cmic.prot_L"),
                            ColM = c("time", "r","E", "Cmic"),
                            Cmic=Cmic.prot.init_L,
                            Niter = 10000)

stopImplicitCluster()

#Models comparison absed on Likelihood ratio test
#this should be correct way
#1-pchisq(-2*(-decay_substrates_ll[1]--decay_all_ll[1]), df=2)

# decay_all_results$likelihood
# decay_structures_results$likelihood
# decay_substrates_results$likelihood
# decay_unique_results$likelihood

an_all_alla$ll
an_all_allb$ll
an_all_allc$ll
an_all_alld$ll

an_substrates_alla$ll
an_substrates_allb$ll
an_substrates_allc$ll
an_substrates_alld$ll

an_structures_alla$ll
an_structures_allb$ll
an_structures_allc$ll
an_structures_alld$ll

an_each_alla$ll
an_each_allb$ll
an_each_allc$ll
an_each_alld$ll

###########################################################################################
#########################################DEB###############################################
###########################################################################################
source("../deb_function.R")

#across all treatments
deb_all<-deb_function(data=dat, SUB = TRUE, FACT = 4, Mean=TRUE, Niter = 10000)
deb_all$ll

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

#for different structures
deb_structures<-deb_function(data=dat, SUB = TRUE, FACT = 2, Mean=TRUE, Niter = 10000)

#for different substrates
deb_substrates<-deb_function(data=dat, SUB = FALSE, FACT = 1, Mean=TRUE, Niter = 10000)

#for each separately 
deb_each<-deb_function(data=dat, SUB = FALSE, FACT = 3, Mean=TRUE, Niter = 10000)

stopImplicitCluster()




deb_structures$ll
deb_substrates$ll
deb_each$ll
