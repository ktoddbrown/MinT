#LIBRARIES
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
library(ggplot2)
library(foreach)
library(doParallel)
library(openxlsx)
library(DEoptim)
library(gridExtra)
library(kableExtra)
library(lmerTest)
library(psycho)
library(ggeffects)
library(rcompanion)
library(vegan)
library(stringr)
###############################################################################################
#GGPLOT THEME
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
###############################################################################################
#DATA
##respiration rates in umol/ml/h
resp<-read.csv(file=c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/MinT/data_checked_respiration_raw.csv"))

##arrange the data
resp.ordered<-resp[order(resp$Sample, resp$Day), ]

##biological data
biology<-read.csv(file=c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/MinT/data_checked_biology_raw.csv"))

##arrange the data in the same way
biology.ordered<-biology[order(biology$Sample, biology$Day), -c(2,3)]
biology.ordered$Sample<-paste0("CSub3_", biology.ordered$Sample)

##merging both data sets
m0<-merge(resp.ordered, biology.ordered, by.x=c("Sample", "Day"), by.y = c("Sample", "Day"), all=T)

##removing rows with NA
m0<-m0[!is.na(m0$Time), ]

##recalculating ptoteins and DNA concentration to umol C - basis
m0$Protinc<-m0$Prot.in*0.46/12.01/4
m0$Protoutc<-m0$Prot.out*0.46/12.01/4
m0$DNAc<-m0$DNA*0.51/12.01/4
m0$DNA.initc<-m0$DNA.init*0.51/12.01/4

##Figures
m0 %>% filter(Substrate=="Celluloze" | Substrate=="Mix") %>% 
  group_by(Day, Structure, Substrate) %>%
  summarize(y=mean(Protinc, na.rm=T),
            y.sd=sd(Protinc, na.rm=T)) %>%
  ggplot(aes(Day, y))+geom_point(cex=6, pch=21, aes(fill=Structure))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))+
  theme_min+facet_wrap(~Substrate, scales="free")

m0 %>% filter(Substrate!="Free") %>% 
  group_by(Day, Structure, Substrate) %>%
  summarize(y=mean(Protoutc, na.rm=T),
            y.sd=sd(Protoutc, na.rm=T)) %>%
  ggplot(aes(Day, y))+geom_point(cex=6, pch=21, aes(fill=Structure))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))+
  theme_min+facet_grid(~Substrate)

m0 %>% filter(Substrate=="Celluloze" | Substrate=="Mix") %>% 
  group_by(Day, Structure, Substrate) %>%
  summarize(y=mean(DNAc, na.rm=T),
            y.sd=sd(DNAc, na.rm=T)) %>%
  ggplot(aes(Day, y))+geom_point(cex=6, pch=21, aes(fill=Structure))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))+
  theme_min+facet_wrap(~Substrate, scales="free")

##Bacteria
bac<-as.data.frame(read.csv(file=c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/MinT/Community_analysis/OTU_bacteria.txt"), header=T,
                            sep="\t"))
###without the taxonomy 
bacr<-bac[, -175]

###transpose and change colnames
bacr<-as.data.frame(t(bacr[, -1]))
colnames(bacr)<-as.character(bac[, 1])


###samples IDs and labels
bac_env<-rownames(bacr)
for(i in 1:length(bac_env)){
  bac_env[i]<-sub("[.]","_", sub("[.]","_", bac_env[i]))
}
rownames(bacr)<-bac_env
bac_env<-as.data.frame(bac_env)
colnames(bac_env)<-c("Sample")

biology.ordered$Sample<-str_sub(biology.ordered$Sample, 7, -1)
biology.ordered$Sample<-gsub("[W]", "", biology.ordered$Sample)
biology.ordered[substr(biology.ordered$Sample, 1, 2)=="MG", "Sample"]<-paste0(substr(biology.ordered[substr(biology.ordered$Sample, 1, 2)=="MG", "Sample"], 1, 1),
                                                    str_sub(biology.ordered[substr(biology.ordered$Sample, 1, 2)=="MG", "Sample"], 3, -1))

m0$Sample<-str_sub(m0$Sample, 7, -1)
m0$Sample<-gsub("[W]", "", m0$Sample)
m0[substr(m0$Sample, 1, 2)=="MG", "Sample"]<-paste0(substr(m0[substr(m0$Sample, 1, 2)=="MG", "Sample"], 1, 1),
                                                                              str_sub(m0[substr(m0$Sample, 1, 2)=="MG", "Sample"], 3, -1))


biology.ordered<-merge(biology.ordered, m0[, c("Sample", "Day", "Structure", "Substrate", "Time")],
                       by=c("Sample", "Day"), all.y = F)

bac_env$Day<-numeric(length = nrow(bac_env))
bac_env$Substrate<-character(length = nrow(bac_env))
bac_env$Structure<-character(length = nrow(bac_env))
bac_env$Time<-numeric(length = nrow(bac_env))

for(i in 1:nrow(bac_env)){
  tryCatch({
    bac_env[i, "Day"]<-as.numeric(biology.ordered[biology.ordered$Sample==bac_env$Sample[i], "Day"])
    bac_env[i, "Substrate"]<-as.character(biology.ordered[biology.ordered$Sample==bac_env$Sample[i], "Substrate"])
    bac_env[i, "Structure"]<-as.character(biology.ordered[biology.ordered$Sample==bac_env$Sample[i], "Structure"])
    bac_env[i, "Time"]<-as.numeric(biology.ordered[biology.ordered$Sample==bac_env$Sample[i], "Time"])
  }, error = function(e){print("Not found")})
  
}

###Initial community structure
bacr_i<-as.data.frame(bacr[substr(rownames(bacr),1,1)=="I", ])

###All community structure without initials and NAs
bacr$Time<-bac_env$Time

bac_env<-bac_env[bac_env$Time>0, ]
bacr_f<-bacr[bacr$Time>0, ]
bacr_f<-bacr_f[, -2134]

#Fungi
fung<-as.data.frame(read.csv(file=c("C:/Users/cape159/Documents/pracovni/data_statistika/minT/MinT/Community_analysis/OTU_fungi.txt"), header=T,
                            sep="\t"))
###without the taxonomy 
fungr<-fung[, -170]

###transpose and change colnames
fungr<-as.data.frame(t(fungr[, -1]))
colnames(fungr)<-as.character(fung[, 1])


###samples IDs and labels
fun_env<-rownames(fungr)
for(i in 1:length(fun_env)){
  fun_env[i]<-sub("[.]","_", sub("[.]","_", fun_env[i]))
}
rownames(fungr)<-fun_env
fun_env<-as.data.frame(fun_env)
colnames(fun_env)<-c("Sample")

fun_env$Day<-numeric(length = nrow(fun_env))
fun_env$Substrate<-character(length = nrow(fun_env))
fun_env$Structure<-character(length = nrow(fun_env))
fun_env$Time<-numeric(length = nrow(fun_env))

for(i in 1:nrow(fun_env)){
  tryCatch({
    fun_env[i, "Day"]<-as.numeric(biology.ordered[biology.ordered$Sample==fun_env$Sample[i], "Day"])
    fun_env[i, "Substrate"]<-as.character(biology.ordered[biology.ordered$Sample==fun_env$Sample[i], "Substrate"])
    fun_env[i, "Structure"]<-as.character(biology.ordered[biology.ordered$Sample==fun_env$Sample[i], "Structure"])
    fun_env[i, "Time"]<-as.numeric(biology.ordered[biology.ordered$Sample==fun_env$Sample[i], "Time"])
  }, error = function(e){print("Not found")})
  
}

###Initial community structure
fungr_i<-as.data.frame(fungr[substr(rownames(fungr),1,1)=="I", ])

###All community structure without initials and NAs
fungr$Time<-fun_env$Time

fun_env<-fun_env[fun_env$Time>0, ]
fungr_f<-fungr[fungr$Time>0, ]
fungr_f<-fungr_f[, -625]

#Statistical analysis
##Does time significantly affect microbial community structure?
#Bacteria across all 
bac.norm<-decostand(bacr_f, method=c("total"))
bac.dist<-vegdist(bac.norm, method = "bray")
bac.ad<-adonis(bac.dist~Substrate+Structure/Substrate+Day, bac_env)

#For Celluloze only
bac.normC<-bac.norm[bac_env$Substrate=="Celluloze", ]
bac_envC<-bac_env[bac_env$Substrate=="Celluloze", ]

bac.distC<-vegdist(bac.normC, method = "bray")
bac.adC<-adonis(bac.distC~Structure+Day, bac_envC)

#Fungi across all 
fun.norm<-decostand(fungr_f, method=c("total"))
fun.dist<-vegdist(fun.norm, method = "bray")
fun.ad<-adonis(fun.dist~Substrate+Structure/Substrate+Day, fun_env)

#For Celluloze only
func.normC<-fun.norm[fun_env$Substrate=="Celluloze", ]
fun_envC<-fun_env[fun_env$Substrate=="Celluloze", ]

fun.distC<-vegdist(func.normC, method = "bray")
fun.adC<-adonis(fun.distC~Structure+Day, fun_envC)


##############################################################################################
#Celluloze - Mixed glass
##How the distances changes with the time

###Bacteria
####Subset the data
bacCB<-bacr_f[(bac_env$Substrate=="Celluloze" & bac_env$Structure=="Mixed glass"), ]
####Add initials
bacCB<-rbind(bacr_i, bacCB)
####Environmentals
bac_envCB<-bac_env[(bac_env$Substrate=="Celluloze" & bac_env$Structure=="Mixed glass"), 
                   c("Sample", "Time", "Day")]
####Add initials
bac_envCB<-rbind(data.frame(Sample=rep("I", 3),
                            Time=rep(0, 3),
                            Day=rep(0, 3)), bac_envCB)

####Remove zeros
z<-numeric()
for(i in 1:ncol(bacCB)){
  if(sum(bacCB[, i])>0){
    z<-append(z, i)
  }else{}
}
bacCB<-bacCB[, z]

bacCB.norm<-decostand(bacCB, method=c("total"))

####Distance matrix
#####Manual method
bac_envCB$D<-numeric(length = nrow(bac_envCB))
for(i in 1:nrow(bacCB.norm)){
  bac_envCB[i, "D"]<-sqrt(sum((bacCB.norm[i, ]-bacCB.norm[1, ]+
                                 bacCB.norm[i, ]-bacCB.norm[2, ]+
                                 bacCB.norm[i, ]-bacCB.norm[3, ])^2))
}

bac_envCB %>% group_by(Day) %>%
  summarise(x=mean(Time),
            y=mean(D),
            y.sd=sd(D)) %>%
  ggplot(aes(x, y)) + geom_point(cex=6, pch=21, fill="grey")+
  theme_min+geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))+
  stat_function(fun=function(x){1.946*x/(80.642+x+x^2/31.159)})

summary(nlsLM(D~a*Time/(b+Time+Time^2/c), data=bac_envCB,
            start = list(a=0.5, b=10, c=0.2)))

####Scaling
bac_envCB$Dsc<-bac_envCB$D/bac_envCB$D[1:3]

bac_envCB %>% group_by(Day) %>%
  summarise(x=mean(Time),
            y=mean(Dsc),
            y.sd=sd(Dsc)) %>%
  ggplot(aes(x, y)) + geom_point(cex=6, pch=21, fill="grey")+
  theme_min+geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))+
  stat_function(fun=function(x){4.221*x/(30.700+x+x^2/106.021)+1})

summary(nlsLM(Dsc~a*Time/(b+Time+Time^2/c)+1, data=bac_envCB,
              start = list(a=3, b=50, c=80)))
        
###Fungi
####Subset the data
funCB<-fungr_f[(fun_env$Substrate=="Celluloze" & fun_env$Structure=="Mixed glass"), ]
####Add initials
funCB<-rbind(fungr_i, funCB)
####Environmentals
fun_envCB<-fun_env[(fun_env$Substrate=="Celluloze" & fun_env$Structure=="Mixed glass"), 
                   c("Sample", "Time", "Day")]
####Add initials
fun_envCB<-rbind(data.frame(Sample=rep("I", 3),
                            Time=rep(0, 3),
                            Day=rep(0, 3)), fun_envCB)

####Remove zeros
z<-numeric()
for(i in 1:ncol(funCB)){
  if(sum(funCB[, i])>0){
    z<-append(z, i)
  }else{}
}
funCB<-funCB[, z]

funCB.norm<-decostand(funCB, method=c("total"))

####Distance matrix
#####Manual method
fun_envCB$D<-numeric(length = nrow(fun_envCB))
for(i in 1:nrow(funCB.norm)){
  fun_envCB[i, "D"]<-sqrt(sum((funCB.norm[i, ]-funCB.norm[1, ]+
                                 funCB.norm[i, ]-funCB.norm[2, ]+
                                 funCB.norm[i, ]-funCB.norm[3, ])^2))
}

fun_envCB %>% group_by(Day) %>%
  summarise(x=mean(Time),
            y=mean(D),
            y.sd=sd(D)) %>%
  ggplot(aes(x, y)) + geom_point(cex=6, pch=21, fill="grey")+
  theme_min+geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))

####Scaling
fun_envCB$Dsc<-fun_envCB$D/mean(fun_envCB$D[1:3])

fun_envCB %>% group_by(Day) %>%
  summarise(x=mean(Time),
            y=mean(Dsc),
            y.sd=sd(Dsc)) %>%
  ggplot(aes(x, y)) + geom_point(cex=6, pch=21, fill="grey")+
  theme_min+geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd))+
  geom_point(cex=6, data=bac_envCB, aes(Time, Dsc))

#####Modeling
m0CB<-subset(m0, Substrate=="Celluloze" & Structure=="Mixed glass")

ggplot(m0CB, aes(Time, r)) + geom_point(cex=6) + theme_min
ggplot(m0CB, aes(Time, Protinc)) + geom_point(cex=6) + theme_min
ggplot(m0CB, aes(Time, DNAc)) + geom_point(cex=6) + theme_min
ggplot(m0CB, aes(Time, Protoutc)) + geom_point(cex=6) + theme_min


####Round the time
m0CB$Time2<-round(m0CB$Time, 0)

####First all model parameters are constant
source("C:/Users/cape159/Documents/pracovni/data_statistika/minT/MinT/Community_analysis/R_Functions/DB_constant.R")

out_const<-DB_constant(m0CB)
out_const$pars
out_const$fit$Gfit

ggplot(out_const$fit$Yhat, aes(time, obs))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~variable, scales="free")+theme_min+
  geom_line(aes(time, yhat))
