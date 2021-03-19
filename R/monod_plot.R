monod_plot<-function(ep){
  #model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      Protc=fp*Cmic
      DNAc=fd*Cmic
      
      dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
      dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
      
      return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C), Protc=Protc, DNAc=DNAc))
      
    })
  }
  #define names of parameters
  parnames<-c("Vmax", "Km", "CUE", "k", "fp", "fd")
  
  #Broth
  par<-as.numeric(ep[ep$Name=="Broth", c(2:7)])
  names(par)<-parnames
  yhatB<-as.data.frame(ode(y=c(Cmic=2.01*0.51/12.01/4/par[["fd"]], C=25), 
                           parms=par, deriv, times=seq(0, 120)))
  yhatB$Struc<-"BROTH"
  #Glass
  par<-as.numeric(ep[ep$Name=="Glass", c(2:7)])
  names(par)<-parnames
  yhatG<-as.data.frame(ode(y=c(Cmic=2.01*0.51/12.01/4/par[["fd"]], C=25), 
                             parms=par, deriv, times=seq(0, 120)))
  yhatG$Struc<-"GLASS"
  #Wool
  par<-as.numeric(ep[ep$Name=="Wool", c(2:7)])
  names(par)<-parnames
  yhatW<-as.data.frame(ode(y=c(Cmic=2.01*0.51/12.01/4/par[["fd"]], C=25), 
                           parms=par, deriv, times=seq(0, 120)))
  yhatW$Struc<-"WOOL"
  
  return(rbind(yhatB, yhatG, yhatW))
}