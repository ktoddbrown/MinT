two_pool_plot<-function(ep){
  #model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #Carbon uptake
      Cu=Vmax*C*S/(Km+C)
      #Respiration
      r=f*R*(1-Yu)
      
      Protc=fpr*R+fps*S
      DNAc=fds*S
      
      dR<-Cu-f*R
      dS<-f*R*Yu - k*S
      dC<--Cu + k*S
      
      
      return(list(c(dR, dS, dC), r=r, Protc=Protc, DNAc=DNAc))
      
    })
  }
  #define names of parameters
  parnames<-c("Vmax","Km", "f","k", "Yu", "fpr", "fps", "fds")
  
  #Broth
  par<-as.numeric(ep[ep$Name=="Broth", c(2:9)])
  names(par)<-parnames
  Sinit=2.01*0.51/12.01/4/par[["fds"]]
  Rinit=0.04012019*Sinit
  yhatB<-as.data.frame(ode(y=c(R=Rinit, S=Sinit, C=25), parms=par, deriv, 
                              times=seq(0, 120)))
  yhatB$Struc<-"BROTH"
  #Glass
  par<-as.numeric(ep[ep$Name=="Glass", c(2:9)])
  names(par)<-parnames
  Sinit=2.01*0.51/12.01/4/par[["fds"]]
  Rinit=0.04012019*Sinit
  yhatG<-as.data.frame(ode(y=c(R=Rinit, S=Sinit, C=25), parms=par, deriv, 
                           times=seq(0, 120)))
  yhatG$Struc<-"GLASS"
  #Wool
  par<-as.numeric(ep[ep$Name=="Wool", c(2:9)])
  names(par)<-parnames
  Sinit=2.01*0.51/12.01/4/par[["fds"]]
  Rinit=0.04012019*Sinit
  yhatW<-as.data.frame(ode(y=c(R=Rinit, S=Sinit, C=25), parms=par, deriv, 
                           times=seq(0, 120)))
  yhatW$Struc<-"WOOL"
  
  return(rbind(yhatB, yhatG, yhatW))
}