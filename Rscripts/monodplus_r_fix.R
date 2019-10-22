monodplus_r_fix<-function(data, FACT){
  
  #FACT 1 = Substrate, 2=Structure, 3=Both, 4=No
  
  if(FACT==1){
    
    dat<-data
    dat$id<-ifelse(d$Substrate=="Glucose", 1, 2)
    
  }else{
    
    if(FACT==2){
      
      dat<-data
      dat<-data
      ids<-dat %>% group_by(Structure) %>% summarise(id=n())
      ids$id<-seq(1:nrow(ids))
      
      dat<-merge(dat, ids, by.x=c("Structure"), by.y=c("Structure"))
      
    }else{
      
      if(FACT==3){
        
        dat<-data
        ids<-dat %>% group_by(Substrate, Structure) %>% summarise(id=n())
        ids$id<-seq(1:nrow(ids))
        
        dat<-merge(dat, ids, by.x=c("Substrate","Structure"), by.y=c("Substrate", "Structure"))
        
        
      }else{
        
        dat<-data
        dat$id<-c(1)
        
      }
    }
  }
  
  #mend model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #equations
      #C uptake
      Cu=Vmax*Cmic*C/(Km+C)
      #respiration
      r=mr*Cmic+(1-CUE)*Cu
      
      #Protc=Cmic*fp
      #DNAc=Cmic*fd
      
      #states
      dCmic<-CUE*Cu-mr*Cmic-k*Cmic
      dC<--Cu + k*Cmic
      
      return(list(c(dCmic, dC), r=r))
      
    })
  }
  #define names of parameters
  parnames<-c("Vmax", "Km", "CUE", "mr", "k")
  
  #parameters estimation function
  estim<-function(odeset){
    
    #defining cost function
    cost<-function(x){
      
      par<-x[1:length(parnames)]
      
      names(par)<-parnames
      
      #first, pars dependent output from ode is matched with measured values
      yhat_all<-as.data.frame(ode(y=c(Cmic=0.09125118, C=25), parms=par, deriv, times=sort(odeset$Time)))
      
      #select time and the measured variables 
      yhat<-select(yhat_all, c("time", "r"))
      colnames(yhat)<-c("time", "value")
      
      #add the measured data to a data frame
      yhat$obs<-c(odeset[order(odeset$Time), c("r")])
      
      #now, the root mean square error is calculated
      NRMSE<-with(yhat, sum((((value-obs))^2), na.rm = T))
      
      return(NRMSE)
      
    }
    
    #defining goodness of fit function 
    rsq_ode<-function(x){
      
      par<-x[1:length(parnames)]
      
      names(par)<-parnames
      
      #first, pars dependent output from ode is matched with measured values
      yhat_all<-as.data.frame(ode(y=c(Cmic=0.09125118, C=25), parms=par, deriv, times=sort(odeset$Time)))
      
      #select time and the measured variables 
      yhat<-select(yhat_all, c("time", "r"))
      colnames(yhat)<-c("time", "value")
      
      #add the measured data to a data frame
      yhat$obs<-c(odeset[order(odeset$Time), c("r")])
      
      #rsquared calculation for each variable
      SSres=with(yhat, sum(((obs-value)^2), na.rm = T))
      SStot=with(yhat, sum(((obs-mean(obs, na.rm = T))^2), na.rm = T))
      ll=with(yhat, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
      R2<-1-SSres/SStot
      N<-length(x)
      AIC<-2*N-2*ll
      
      Gfit<-c(R2=R2, N=N, AIC=AIC, ll=ll)
      
      rsq_out<-list(Yhat=yhat, Gfit=Gfit)
      
      return(rsq_out)
      
    }
    
    #approximate parameter estimation is done by MCMC method
    par_mcmc<-modMCMC(f=cost, p=c(Vmax=0.1, Km=3, CUE=0.5, mr=0.01, k=1e-3), 
                      lower=c(Vmax=1e-3, Km=1e-3, CUE=0, mr=1e-5, k=1e-6),
                      upper=c(Vmax=10, Km=100, CUE=1, mr=10, k=10), niter=10000)
    
    #lower and upper limits for parameters are extracted
    pl<-summary(par_mcmc)["min",]
    pu<-summary(par_mcmc)["max",]
    
    #these limits are used to find global optimum by DEoptim
    opt_par<-DEoptim(fn=cost, lower=pl, upper=pu, 
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
                                 trace=FALSE, strategy=3, NP=250))
    
    #global optimum parameters are further used in MCMC to find parameters distribution
    par_prof<-modMCMC(f=cost, p=opt_par$optim$bestmem, 
                      lower=pl,
                      upper=pu, niter=5000)
    
    #goodness of fit
    fit<-rsq_ode(opt_par$optim$bestmem)
    
    #best parameters
    p<-opt_par$optim$bestmem
    names(p)<-parnames
    
    #return list with opt_par and par_prof
    estim_out<-list(pars=p, par_prof=par_prof, fit=fit)
    
    return(estim_out)
    
  }
  
  
  #parameter estimation
  if(FACT==4){
    
    res<-vector("list", length = 1)
    res[[1]]<-estim(odeset=dat)
    
  }else{
    
    res<-foreach(i=unique(dat$id), .combine=list, .multicombine = TRUE,
                 .packages=c("FME", "dplyr", "DEoptim", "reshape")) %dopar% {
                   
                   estim(odeset=dat[dat$id==i,])
                   
                 }
  }
  
  
  #calculation of overall goodness of fit from individual results
  if(FACT==4){
    
    res$goodness<-res[[1]]$fit$Gfit
    res$OvP<-res[[1]]$fit$Yhat
    
  }else{
    
    if(FACT==1){
      
      res$OvP<-rbind(res[[1]]$fit$Yhat, res[[2]]$fit$Yhat)
      
      #rsquared calculation for each variable
      SSres=with(res$OvP, sum(((obs-value)^2), na.rm = T))
      SStot=with(res$OvP, sum(((obs-mean(obs, na.rm = T))^2), na.rm = T))
      ll=with(res$OvP, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
      R2<-1-SSres/SStot
      N<-length(parnames)*2+1
      AIC<-2*N-2*ll
      
      res$goodness<-c(R2=R2, N=N, AIC=AIC, ll=ll)
      
    }else{
      
      if(FACT==2){
        
        res$OvP<-rbind(res[[1]]$fit$Yhat, res[[2]]$fit$Yhat, res[[3]]$fit$Yhat)
        
        #rsquared calculation for each variable
        SSres=with(res$OvP, sum(((obs-value)^2), na.rm = T))
        SStot=with(res$OvP, sum(((obs-mean(obs, na.rm = T))^2), na.rm = T))
        ll=with(res$OvP, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
        R2<-1-SSres/SStot
        N<-length(parnames)*3+1
        AIC<-2*N-2*ll
        
        res$goodness<-c(R2=R2, N=N, AIC=AIC, ll=ll)
        
      }else{
        
        res$OvP<-rbind(res[[1]]$fit$Yhat, res[[2]]$fit$Yhat, res[[3]]$fit$Yhat, res[[4]]$fit$Yhat, res[[5]]$fit$Yhat, res[[6]]$fit$Yhat)
        
        #rsquared calculation for each variable
        SSres=with(res$OvP, sum(((obs-value)^2), na.rm = T))
        SStot=with(res$OvP, sum(((obs-mean(obs, na.rm = T))^2), na.rm = T))
        ll=with(res$OvP, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
        R2<-1-SSres/SStot
        N<-length(parnames)*6+1
        AIC<-2*N-2*ll
        
        res$goodness<-c(R2=R2, N=N, AIC=AIC, ll=ll)
      }
      
    }
    
  }
  
  return(res)
  
}