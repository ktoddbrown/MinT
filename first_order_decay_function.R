###first order decay
first_order_function<-function(data, SUB, FACT, Vars, Niter){
  
  #FACT 1 = Substrate, 2=Structure, 3=Both, 4=No
  
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
  #############################################################################################
  #if parameters are estimated across different substrates, different initial C conc
  #has to be aknowledged
  
  if(SUB==TRUE){
    
    #function I need
    first_order<-function(X, pars, t){
      
      #first order decay function
      deriv<-function(time, state, pars){
        
        with(as.list(c(state, pars)),{
          
          dC<--k*C
          
          return(list(c(dC), r=k*C))
          
        })
      }
      
      #function to apply acros 3 different initial substrate concentrations
      base_function<-function(Ci){
        
        
        out<-ode(y=c(C=Ci), parms=pars, times=t, func=deriv)
        return(as.data.frame(out))}
      
      #results
      res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
      
      return(res)
      
    }
    
    #parameter estimation function
    estim<-function(data){
      
      #Obs_decay<-dat[,c("Substrate", "r", "Time")]
      Obs_decay<-select(data, Vars)
      
      #m1<-merge(Obs_decay[Obs_decay$Substrate=="Glucose", c("r", "Time")], 
                #Obs_decay[Obs_decay$Substrate=="Cellobiose", c("r", "Time")],all = T, by="Time")
      
      m1<-merge(select(filter(Obs_decay, Substrate=="Glucose"), Vars[-1]),
                select(filter(Obs_decay, Substrate=="Cellobiose"), Vars[-1]),all = T, by="Time")
        
      m2<-merge(m1, 
                select(filter(Obs_decay, Substrate=="Cellobiose"), Vars[-1]),all = T, by="Time")
      
      #m2<-merge(m1, 
                #Obs_decay[Obs_decay$Substrate=="Mix", c("r", "Time")],all = T, by="Time")
      
      #Glucose = r, Cellobiose = r1, Mix = r2 
      colnames(m2)<-c("time", "r", "r1", "r2")
      
      #create cost function
      cost_decay_SUB<-function(pars){
        
        out<-first_order(X=c(25+0.1244, 25+0.1244, 16.5+0.1244), pars = pars, t=seq(0,130))
        cost<-modCost(model = out, obs = m2)
        
        return(cost)
        
      }
      
      mcmc<-modMCMC(f=cost_decay_SUB, p=c(k=0.01), niter=Niter)
      
      return(mcmc)
      
    }
    
    #parameter estimation
    if(FACT==4){
      res=vector("list", length = 1)
      res[[1]]<-estim(data)
    }else{
      res<-foreach(i=unique(dat$id), .combine=list, .multicombine = TRUE,
                   .packages=c("FME", "dplyr")) %dopar% {estim(data=dat[dat$id==i,])}
    }
    
    
    #parameters extraction
    #data frame with no values is created and populated
    parameters<-data.frame(k50=vector("numeric", length = length(unique(dat$id))),
                           k25=vector("numeric", length = length(unique(dat$id))),
                           k75=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    for(i in unique(dat$id)){
      
      parameters[i, "k50"]<-summary(res[[i]])[6,1]
      parameters[i, "k25"]<-summary(res[[i]])[5,1]
      parameters[i, "k75"]<-summary(res[[i]])[7,1]
      
    }
    
    if(FACT==2){
      
      parameters$Structure<-ddply(dat, .(id,Structure), summarize, mean(r))[,2]
      
    }else{
      
      
    }
    
    #combining results
    cost_decay_SUB2<-function(pars, data){
      
      Obs_decay<-select(data, Vars)
      
      m1<-merge(select(filter(Obs_decay, Substrate=="Glucose"), Vars[-1]),
                select(filter(Obs_decay, Substrate=="Cellobiose"), Vars[-1]),all = T, by="Time")
      
      m2<-merge(m1, 
                select(filter(Obs_decay, Substrate=="Cellobiose"), Vars[-1]),all = T, by="Time")
      
      #Glucose = r, Cellobiose = r1, Mix = r2 
      colnames(m2)<-c("time", "r", "r1", "r2")
      
      out<-first_order(X=c(25+0.1244, 25+0.1244, 16.5+0.1244), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2)
      
      return(cost)
      
    }
    
    obs<-numeric()
    mod<-numeric()
    
    for(i in unique(dat$id)){
      
      obs<-append(obs, cost_decay_SUB2(pars = c(k=parameters[i,"k50"]), data=dat[dat$id==i,])$residuals$obs)
      mod<-append(mod, cost_decay_SUB2(pars = c(k=parameters[i,"k50"]), data=dat[dat$id==i,])$residuals$mod)
    }
    
    OvP<-data.frame(obs, mod)
    
    #logLik calculation
    mu<-mean(obs)
    variance<-sd(obs)^2
    
    
    ll<--1*sum((obs-mod)^2)/2/variance
    
    
    SSmodel<-sum((obs-mod)^2)
    SSdata<-sum((obs-mean(obs))^2)
    
    rsq=1-(SSmodel/SSdata)
    
    likelihood<-c(logLik=ll, npar=length(unique(dat$id)), rsq=rsq, AIC=2*length(unique(dat$id))-2*ll)
    
    al<-list()
    al$parameters<-parameters
    al$OvP<-OvP
    al$likelihood<-likelihood
    
    
    
    return(al)
    
    
    
  }else{
    
    #first order decay function
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dC<--k*C
        
        return(list(c(dC), r=k*C))
        
      })
    }
    
    #parameter estimation function
    estim<-function(data){
      
      Obs_dat<-select(data, Vars)
      colnames(Obs_dat)<-c("time", "r")
      cinit<-as.numeric(dat[1,"DOCinit"])+0.1244
      
      
      #cost function
      cost_function<-function(pars){
        
        out<-as.data.frame(ode(y=c(C=cinit), parms = pars, t=seq(0,130), func = deriv))
        cost<-modCost(model = out, obs = Obs_dat)
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(k=0.01), niter=Niter)
      
      
      return(res)
      
    }
    
    #parameter estimation
    res<-foreach(i=unique(dat$id), .combine=list, .multicombine = TRUE,
                 .packages=c("FME", "dplyr")) %dopar% {
                   
                   estim(data=dat[dat$id==i,])
                   
                 }
    
    
    #parameters extraction
    #data frame with no values is created and populated
    parameters<-data.frame(k50=vector("numeric", length = length(unique(dat$id))),
                           k25=vector("numeric", length = length(unique(dat$id))),
                           k75=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    for(i in unique(dat$id)){
      
      parameters[i, "k50"]<-summary(res[[i]])[6,1]
      parameters[i, "k25"]<-summary(res[[i]])[5,1]
      parameters[i, "k75"]<-summary(res[[i]])[7,1]
      
    }
    
    if(FACT==1){
      
      parameters$Substrate<-ddply(dat, .(id,Substrate), summarize, mean(r))[,2]
      
    }else{
      
      parameters$Substrate<-ddply(dat, .(id,Substrate, Structure), summarize, mean(r))[,2]
      parameters$Structure<-ddply(dat, .(id,Substrate, Structure), summarize, mean(r))[,3]
    }
    
    #OvP
    
    obs<-numeric()
    mod<-numeric()
    
    cost_function2<-function(pars, data){
      
      Obs_dat<-select(data, Vars)
      colnames(Obs_dat)<-c("time", "r")
      cinit<-as.numeric(data[1,"DOCinit"])+0.1244
      
      out<-as.data.frame(ode(y=c(C=cinit), parms = pars, t=seq(0,130), func = deriv))
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      
      obs<-append(obs, cost_function2(pars = c(k=as.numeric(parameters[i,"k50"])), data=dat[dat$id==i,])$residuals$obs)
      mod<-append(mod, cost_function2(pars = c(k=as.numeric(parameters[i,"k50"])), data=dat[dat$id==i,])$residuals$mod)
    }
    
    
    OvP<-data.frame(obs, mod)
    
    #logLik calculation
    mu<-mean(obs)
    variance<-sd(obs)^2
    
    
    ll<--1*sum((obs-mod)^2)/2/variance
    
    
    SSmodel<-sum((obs-mod)^2)
    SSdata<-sum((obs-mean(obs))^2)
    
    rsq=1-(SSmodel/SSdata)
    
    likelihood<-c(logLik=ll, npar=length(unique(dat$id)), rsq=rsq, AIC=2*length(unique(dat$id))-2*ll)
    
    al<-list()
    al$parameters<-parameters
    al$OvP<-OvP
    al$likelihood<-likelihood
    
    
    
    return(al)
    
  }
  
}