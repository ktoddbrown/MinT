mem_init<-function(data, FACT, Mean, Niter, DNAci){
  
  #SUB - needs to be calibrated in matrix
  #FACT 1 = Substrate, 2=Structure, 3=Both, 4=No
  #Vars = names of variables to calibrated the model against
  #ColM = colnames for these varaibles in matrix
  #ColV = colnames for these varaibles in normal setting
  #VarsCmic = which biomass to use for likelihood calculation (it should follow the same notation as Vars)
  #Cmic = initial microbial biomass C
  
  #if factor needs to be included, its levels are defined
  if(FACT==1){
    
    dat<-setDT(data)[, id := .GRP, by = .(Substrate)]
    
  }else{
    
    if(FACT==2){
      
      dat<-setDT(data)[, id := .GRP, by = .(Structure)]
    }else{
      
      if(FACT==3){
        
        dat<-setDT(data)[, id := .GRP, by = .(Substrate,Structure)]
      }else{
        dat<-data
        dat$id<-c(1)
        
      }
    }
  }
  
    #mem model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dCmic<--kmic*Cmic+CUE*Vmax*Cmic*C/(Km+C)-pe*Cmic
        dC<-kmic*Cmic-Vmax*Cmic*C/(Km+C)+ke*E
        dE<-pe*Cmic-ke*E
        
        #DNA and Protein content in biomass
        DNAc=fd*Cmic
        Protc=fp*Cmic
        
        return(list(c(dCmic, dC, dE), DNAc=DNAc, Protc=Protc, r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    
    #create cost function
    estim<-function(data){
      
      Obs_dat<-select(data, c("Time", "r","E", "DNAc", "Protc"))
      colnames(Obs_dat)<-c("time", "r","E", "DNAc", "Protc")
      mtrue<-rbind(Obs_dat, data.frame(time=0, r=NA, E=NA, DNAc=DNAci, Protc=NA))  
      
      cost_function<-function(pars){
        
        
        out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
        cost<-modCost(model = out, obs = mtrue)
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, kmic=0.01, ke=0.01, pe=0.01, fd=0.05, fp=0.5, Cmic_0=0.12),
                   lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, kmic=1e-5, ke=1e-6, pe=1e-6, fd=0, fp=0, Cmic_0=1e-6),
                   upper=c(Vmax=1e4, Km=1e4, CUE=0.999, kmic=1e5, ke=1e6, pe=1e6, fd=1, fp=1, Cmic_0=25),niter=Niter)
      
      return(res)
      
    }
    
    #parameter estimation
    if(FACT==4){
      
      res<-vector("list", length = 1)
      res[[1]]<-estim(data=dat)
      
    }else{
      
      res<-foreach(i=unique(dat$id), .combine=list, .multicombine = TRUE,
                   .packages=c("FME", "dplyr")) %dopar% {
                     
                     estim(data=dat[dat$id==i,])
                     
                   }
    }
    
    
    
    
    #parameters extraction
    parameters<-data.frame(Vmax=vector("numeric", length = length(unique(dat$id))), 
                           Km=vector("numeric", length = length(unique(dat$id))), 
                           CUE=vector("numeric", length = length(unique(dat$id))), 
                           kmic=vector("numeric", length = length(unique(dat$id))),
                           ke=vector("numeric", length = length(unique(dat$id))),
                           pe=vector("numeric", length = length(unique(dat$id))),
                           fd=vector("numeric", length = length(unique(dat$id))),
                           fp=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0=vector("numeric", length = length(unique(dat$id))),
                           Vmax.l=vector("numeric", length = length(unique(dat$id))), 
                           Km.l=vector("numeric", length = length(unique(dat$id))), 
                           CUE.l=vector("numeric", length = length(unique(dat$id))), 
                           kmic.l=vector("numeric", length = length(unique(dat$id))),
                           ke.l=vector("numeric", length = length(unique(dat$id))),
                           pe.l=vector("numeric", length = length(unique(dat$id))),
                           fd.l=vector("numeric", length = length(unique(dat$id))),
                           fp.l=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0.l=vector("numeric", length = length(unique(dat$id))),
                           Vmax.u=vector("numeric", length = length(unique(dat$id))), 
                           Km.u=vector("numeric", length = length(unique(dat$id))), 
                           CUE.u=vector("numeric", length = length(unique(dat$id))), 
                           kmic.u=vector("numeric", length = length(unique(dat$id))),
                           ke.u=vector("numeric", length = length(unique(dat$id))),
                           pe.u=vector("numeric", length = length(unique(dat$id))),
                           fd.u=vector("numeric", length = length(unique(dat$id))),
                           fp.u=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0.u=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    
    #median
    for(i in unique(dat$id)){
      for(n in 1:((ncol(parameters)-1)/3)){
        
        parameters[i, n]<-res[[i]]$bestpar[n]
      }
    }
    #lower quartile
    for(i in unique(dat$id)){
      for(n in 1:((ncol(parameters)-1)/3)){
        
        parameters[i, n+((ncol(parameters)-1)/3)]<-summary(res[[i]])[5,n]
      }
    }
    #lower quartile
    for(i in unique(dat$id)){
      for(n in 1:((ncol(parameters)-1)/3)){
        
        parameters[i, n+2*((ncol(parameters)-1)/3)]<-summary(res[[i]])[7,n]
      }
    }
    
   
    labeles<-dat[,c("Time", "Substrate", "Structure")]
    
    #OvP for respiration
    obs_r<-numeric()
    mod_r<-numeric()
    time_r<-numeric()
    
    
    cost_r<-function(pars, data){
      
      Obs_datr<-select(data, c("Time", "r"))
      colnames(Obs_datr)<-c("time", "r")
      
      out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_datr)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      
      obs_r<-append(obs_r, cost_r(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$obs)
      mod_r<-append(mod_r, cost_r(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$mod)
      time_r<-append(time_r, cost_r(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$x)
      
    }
    
    OvP_r<-data.frame(obs_r, mod_r, time_r)
    
    #add Substrate and structure
    OvP_r<-merge(OvP_r, labeles, by.x="time_r", by.y="Time")
    
    #logLik calculation
    mu_r<-mean(obs_r)
    variance_r<-sd(obs_r)^2
    
    
    ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
    
    
    SSmodel_r<-sum((obs_r-mod_r)^2)
    SSdata_r<-sum((obs_r-mean(obs_r))^2)
    
    rsq_r=1-(SSmodel_r/SSdata_r)
    
    likelihood_r<-c(logLik=ll_r, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_r, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_r)
    
    
    #OvP for DNAc
    obs_DNAc<-numeric()
    mod_DNAc<-numeric()
    time_DNAc<-numeric()
    
    
    cost_DNAc<-function(pars, data){
      
      Obs_datd<-select(data, c("Time", "r","E", "DNAc", "Protc"))
      colnames(Obs_datd)<-c("time", "r","E", "DNAc", "Protc")
      #mtrue<-rbind(Obs_dat, data.frame(time=0, r=NA, E=NA, DNAc=DNAci, Protc=NA))  
      Obs_datd<-select(Obs_datd, c("time", "DNAc"))
      
      out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_datd)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_DNAc<-append(obs_DNAc, cost_DNAc(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$obs)
      mod_DNAc<-append(mod_DNAc, cost_DNAc(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$mod)
      time_DNAc<-append(time_DNAc, cost_DNAc(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$x)
      
    }
    
    OvP_DNAc<-data.frame(obs_DNAc, mod_DNAc,time_DNAc)
    
    #add Substrate and structure
    OvP_DNAc<-merge(OvP_DNAc, labeles, by.x="time_DNAc", by.y="Time")
    
    #logLik calculation
    mu_DNAc<-mean(obs_DNAc)
    variance_DNAc<-sd(obs_DNAc)^2
    
    
    ll_DNAc<--1*sum((obs_DNAc-mod_DNAc)^2)/2/variance_DNAc
    
    
    SSmodel_DNAc<-sum((obs_DNAc-mod_DNAc)^2)
    SSdata_DNAc<-sum((obs_DNAc-mean(obs_DNAc))^2)
    
    rsq_DNAc=1-(SSmodel_DNAc/SSdata_DNAc)
    
    likelihood_DNAc<-c(logLik=ll_DNAc, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_DNAc, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_DNAc)
    
    #OvP for Protc
    obs_Protc<-numeric()
    mod_Protc<-numeric()
    time_Protc<-numeric()
    
    cost_Protc<-function(pars, data){
      
      
      Obs_datp<-select(data, c("Time", "Protc"))
      colnames(Obs_datp)<-c("time", "Protc")
      
      out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_datp)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_Protc<-append(obs_Protc, cost_Protc(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$obs)
      mod_Protc<-append(mod_Protc, cost_Protc(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$mod)
      time_Protc<-append(time_Protc, cost_Protc(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$x)
    }
    
    OvP_Protc<-data.frame(obs_Protc, mod_Protc,time_Protc)
    
    #add Substrate and structure
    OvP_Protc<-merge(OvP_Protc, labeles, by.x="time_Protc", by.y="Time")
    
    #logLik calculation
    mu_Protc<-mean(obs_Protc)
    variance_Protc<-sd(obs_Protc)^2
    
    
    ll_Protc<--1*sum((obs_Protc-mod_Protc)^2)/2/variance_Protc
    
    
    SSmodel_Protc<-sum((obs_Protc-mod_Protc)^2)
    SSdata_Protc<-sum((obs_Protc-mean(obs_Protc))^2)
    
    rsq_Protc=1-(SSmodel_Protc/SSdata_Protc)
    
    likelihood_Protc<-c(logLik=ll_Protc, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_Protc, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_Protc)
    
    
    #OvP for enzymes
    obs_E<-numeric()
    mod_E<-numeric()
    time_E<-numeric()
    
    cost_E<-function(pars, data){
      
      Obs_date<-select(data, c("Time", "E"))
      colnames(Obs_date)<-c("time", "E")
      
      out<-ode(y=c(Cmic=pars[["Cmic_0"]], C=25, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_date)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_E<-append(obs_E, cost_E(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$obs)
      mod_E<-append(mod_E, cost_E(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$mod)
      time_E<-append(time_E, cost_E(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$x)
      
    }
    
    OvP_E<-data.frame(obs_E, mod_E,time_E)
    
    #add Substrate and structure
    OvP_E<-merge(OvP_E, labeles, by.x="time_E", by.y="Time")
    
    #logLik calculation
    mu_E<-mean(obs_E)
    variance_E<-sd(obs_E)^2
    
    
    ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
    
    
    SSmodel_E<-sum((obs_E-mod_E)^2)
    SSdata_E<-sum((obs_E-mean(obs_E))^2)
    
    rsq_E=1-(SSmodel_E/SSdata_E)
    
    likelihood_E<-c(logLik=ll_E, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_E, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_E)
    
    al<-list()
    
    if(FACT==1){
      
      parameters$Substrate<-ddply(dat, .(id,Substrate), summarize, mean(r))[,2]
      
    }else{
      
      if(FACT==4){
        
        
      }else{
        
        if(FACT==3){
          
          parameters$Substrate<-ddply(dat, .(id,Substrate, Structure), summarize, mean(r))[,2]
          parameters$Structure<-ddply(dat, .(id,Substrate, Structure), summarize, mean(r))[,3]
          
        }else{
          
          parameters$Structure<-ddply(dat, .(id,Structure), summarize, mean(r))[,2]  
        }
      }
    }
    
    al$parameters<-parameters
    al$OvP_r<-OvP_r
    al$ll<-rbind(likelihood_r, likelihood_DNAc, likelihood_Protc, likelihood_E)
    al$OvP_DNAc<-OvP_DNAc
    al$OvP_Protc<-OvP_Protc
    al$OvP_E<-OvP_E
    
    return(al)
    
  
  
}