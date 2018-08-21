mem_function<-function(data, FACT, Vars, ColM, Mean, Cmic, Niter){
  
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
        
        return(list(c(dCmic, dC, dE), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    
    #create cost function
    estim<-function(data){
      
      Obs_dat<-select(data, Vars)
      colnames(Obs_dat)<-c("time", "r","E", "Cmic")
      mtrue<-select(Obs_dat, ColM)  
      
      cost_function<-function(pars){
        
        
        out<-as.data.frame(ode(y=c(Cmic=Cmic, C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
        cost<-modCost(model = out, obs = mtrue)
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, kmic=0.01, ke=0.01, pe=0.01),
                   lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, kmic=1e-5, ke=1e-6, pe=1e-6),
                   upper=c(Vmax=1e4, Km=1e4, CUE=0.999, kmic=1e5, ke=1e6, pe=1e6),niter=Niter)
      
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
                           Vmax.l=vector("numeric", length = length(unique(dat$id))), 
                           Km.l=vector("numeric", length = length(unique(dat$id))), 
                           CUE.l=vector("numeric", length = length(unique(dat$id))), 
                           kmic.l=vector("numeric", length = length(unique(dat$id))),
                           ke.l=vector("numeric", length = length(unique(dat$id))),
                           pe.l=vector("numeric", length = length(unique(dat$id))),
                           Vmax.u=vector("numeric", length = length(unique(dat$id))), 
                           Km.u=vector("numeric", length = length(unique(dat$id))), 
                           CUE.u=vector("numeric", length = length(unique(dat$id))), 
                           kmic.u=vector("numeric", length = length(unique(dat$id))),
                           ke.u=vector("numeric", length = length(unique(dat$id))),
                           pe.u=vector("numeric", length = length(unique(dat$id))),
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
        
        parameters[i, n+((ncol(parameters)-1)/3)*2]<-summary(res[[i]])[7,n]
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
      
      out<-as.data.frame(ode(y=c(Cmic=Cmic, C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
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
    
    
    #OvP for biomass
    obs_Cmic<-numeric()
    mod_Cmic<-numeric()
    time_Cmic<-numeric()
    
    cost_Cmic<-function(pars, data){
      
      Obs_dat<-select(data, Vars[c(1,4)])
      
      colnames(Obs_dat)<-c("time", "Cmic")
      
      out<-as.data.frame(ode(y=c(Cmic=Cmic, C=25, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$obs)
      mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$mod)
      time_Cmic<-append(time_Cmic, cost_Cmic(pars=res[[i]]$bestpar, data=dat[dat$id==i, ])$residuals$x)
    }
    
    OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
    
    #add Substrate and structure
    OvP_Cmic<-merge(OvP_Cmic, labeles, by.x="time_Cmic", by.y="Time")
    
    #logLik calculation
    mu_Cmic<-mean(obs_Cmic)
    variance_Cmic<-sd(obs_Cmic)^2
    
    
    ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
    
    
    SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
    SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
    
    rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
    
    likelihood_Cmic<-c(logLik=ll_Cmic, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_Cmic, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_Cmic)
    
    #OvP for enzymes
    obs_E<-numeric()
    mod_E<-numeric()
    time_E<-numeric()
    
    cost_E<-function(pars, data){
      
      Obs_date<-select(data, c("Time", "E"))
      colnames(Obs_date)<-c("time", "E")
      
      out<-ode(y=c(Cmic=Cmic, C=25, E=0), parms=pars, times=seq(0,130), func=deriv)
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
    al$ll<-rbind(likelihood_r, likelihood_Cmic, likelihood_E)
    al$OvP_Cmic<-OvP_Cmic
    al$OvP_E<-OvP_E
    
    return(al)
    
  
  
}