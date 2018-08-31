andrew_init_a<-function(data, FACT){
  
  #SUB - needs to be calibrated in matrix
  #FACT 1 = Substrate, 2=Structure, 3=Both, 4=No
  #DNAci = initial DNA concentration
  
  #if factor needs to be included, its levels are defined
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
        
        #specific growth rate
        u=umax*C/(Km+C+Ki*C^2)
        
        #respiration rate
        r=m*Cmic+(1/Ymic-1)*Cmic*u+(1/Ye-1)*(Cmic*u*alfa)
        
        
        dCmic<-Cmic*u-kmic*Cmic
        dC<--(1/Ymic)*Cmic*u-(1/Ye)*(Cmic*u*alfa)-m*Cmic+kmic*Cmic+ke*E
        dE<-Cmic*u*alfa-ke*E
        
        #DNA and Protein content in biomass
        DNAc=fd*Cmic
        # Protc=fp*Cmic
        
        return(list(c(dCmic, dC, dE), DNAc=DNAc, r=r))
        
      })
    }
    #define names of parameters
    parnames<-c("umax", "Km", "Ki", "Ymic", "Ye", "alfa", "m", "kmic", "ke", "Cmic_0", "fd")
    
    #parameters estimation function
    estim<-function(odeset){
      
      #defining cost function
      cost<-function(x){
        
        par<-x[1:length(parnames)]
        
        names(par)<-parnames
        
        #first, pars dependent output from ode is matched with measured values
        yhat_all<-as.data.frame(ode(y=c(Cmic=par[["Cmic_0"]], C=25, E=0), parms=par, deriv, times=sort(odeset$Time)))
        
        #select time and the measured variables 
        yhat<-select(yhat_all, c("time", "DNAc", "r", "E"))
        
        #reformat to long format data frame
        Yhat<-melt(yhat, id.vars = "time")
        
        #add the measured data to a data frame
        Yhat$obs<-c(odeset[order(odeset$Time), "DNAc"], odeset[order(odeset$Time), c("r")], odeset[order(odeset$Time), "E"])
        
        #add the weighting factor
        #I want to have the weighting factor to be proportional to mean of the given variable 
        yweights<-Yhat %>% group_by(variable) %>% summarise(weights=mean(value, na.rm = T))
        
        #match with the Yhat data frame
        Yhat<-merge(Yhat, yweights, by.x=c("variable"), by.y=c("variable"))
        
        #now, the root mean square error is calculated
        NRMSE<-as.numeric(Yhat %>% group_by(variable) %>% summarise(NRMSE=sum((((value-obs)/mean(weights))^2), na.rm = T)) %>% 
                            summarise(NRMSE=sum(NRMSE)) )
        
        return(NRMSE)
        
      }
      
      #defining goodness of fit function 
      rsq_ode<-function(x){
        
        par<-x[1:length(parnames)]
        
        names(par)<-parnames
        
        #first, pars dependent output from ode is matched with measured values
        yhat_all<-as.data.frame(ode(y=c(Cmic=par[["Cmic_0"]], C=25, E=0), parms=par, deriv, times=sort(odeset$Time)))
        
        #select time and the measured variables 
        yhat<-select(yhat_all, c("time", "DNAc", "r", "E"))
        
        #reformat to long format data frame
        Yhat<-melt(yhat, id.vars = "time")
        
        #add the measured data to a data frame
        Yhat$obs<-c(odeset[order(odeset$Time), "DNAc"], odeset[order(odeset$Time), c("r")], odeset[order(odeset$Time), "E"])
        Yhat$Substrate<-rep(odeset[order(odeset$Time), "Substrate"], times=3)
        Yhat$Structure<-rep(odeset[order(odeset$Time), "Structure"], times=3)
        
        #rsquared calculation for each variable
        Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T), 
                                                        SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                        ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
        Gfit$R2<-with(Gfit, 1-SSres/SStot)
        Gfit$N<-c(11)
        Gfit$AIC<-with(Gfit, 2*N-2*ll)
        
        rsq_out<-list(Yhat=Yhat, Gfit=Gfit)
        
        return(rsq_out)
        
      }
      
      
      #approximate parameter estimation is done by MCMC method
      par_mcmc<-modMCMC(f=cost, p=c(0.1, 3, 3, 0.8, 0.5, 0.1, 0.01, 0.01, 0.01, 0.1, 0.004), 
                          lower=c(1e-3, 1e-3,1e-5, 0.01, 0.01, 1e-5, 1e-5, 1e-5, 1e-5, 0.01, 0.004),
                          upper=c(1, 10, 10, 1, 1, 1, 1, 0.1, 0.1, 2, 0.5), niter=10000)
      
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
        res$goodness<-res$OvP %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T), 
                                                               SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                               ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
        res$goodness$R2<-with(res$goodness, 1-SSres/SStot)
        res$goodness$N<-rep(length(parnames)*2, times=3)
        res$goodness$AIC<-with(res$goodness, 2*N-2*ll)
        
        
      }else{
        
        if(FACT==2){
          
          res$OvP<-rbind(res[[1]]$fit$Yhat, res[[2]]$fit$Yhat, res[[3]]$fit$Yhat)
          
          #rsquared calculation for each variable
          res$goodness<-res$OvP %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T), 
                                                                 SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                 ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
          res$goodness$R2<-with(res$goodness, 1-SSres/SStot)
          res$goodness$N<-rep(length(parnames)*3, times=3)
          res$goodness$AIC<-with(res$goodness, 2*N-2*ll)
          
        }else{
          
          res$OvP<-rbind(res[[1]]$fit$Yhat, res[[2]]$fit$Yhat, res[[3]]$fit$Yhat, res[[4]]$fit$Yhat, res[[5]]$fit$Yhat, res[[6]]$fit$Yhat)
          
          #rsquared calculation for each variable
          res$goodness<-res$OvP %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T), 
                                                                 SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                 ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
          res$goodness$R2<-with(res$goodness, 1-SSres/SStot)
          res$goodness$N<-rep(length(parnames)*6, times=3)
          res$goodness$AIC<-with(res$goodness, 2*N-2*ll)
          
        }
        
      }
      
    }
    
  return(res)
  
}