andrew_function<-function(data, SUB, FACT, Vars, ColM, Mean, Cmic, Niter){
  
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
  
  ############################################################################################
  #if parameters are estimated across different substrates, different initial C conc
  #has to be acknowledged
  
  if(SUB==TRUE){
    
    #this is the mmem model fitted across different substrates 
    #with different initial carbon concentrations
    mmem<-function(X, pars, t){
      
      #mmem model
      deriv<-function(time, state, pars){
        
        with(as.list(c(state, pars)),{
          
          #specific growth rate
          
          u=umax*C/(Km+C+Ki*C^2)
          
          #respiration rate
          r=m*Cmic+(1/Ymic-1)*Cmic*u+(1/Ye-1)*(Cmic*u*alfa)
          
          
          dCmic<-Cmic*u-kmic*Cmic
          dC<--(1/Ymic)*Cmic*u-(1/Ye)*(Cmic*u*alfa)-m*Cmic+kmic*Cmic+ke*E
          dE<-Cmic*u*alfa-ke*E
          
          return(list(c(dCmic, dC, dE), r=r))
          
        })
      }
      
      #function to apply acros 3 different initial substrate concentrations
      base_function<-function(Ci){
        
        
        out<-ode(y=c(Cmic=Cmic, C=Ci, E=0), parms=pars, times=t, func=deriv)
        return(as.data.frame(out))}
      
      #results
      res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
      
      return(res)
      
    }
    
    #create cost function
    estim<-function(data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), Vars)
      colnames(gl)<-c("time", "r", "E" ,"Cmic")
      #round the time to one digit
      gl$time<-round(gl$time, 1)
      
      gl2<-select(filter(data, Substrate!="Glucose" ), Vars[1])
      colnames(gl2)<-"time"
      #round the time to one digit
      gl2$time<-round(gl2$time, 1)
      
      gl2$r<-NA
      gl2$E<-NA
      gl2$Cmic<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), Vars)
      colnames(cel)<-c("time", "r1", "E1" ,"Cmic1")
      #round the time to one digit
      cel$time<-round(cel$time, 1)
      
      cel2<-select(filter(data, Substrate!="Cellobiose" ), Vars[1])
      colnames(cel2)<-"time"
      #round the time to one digit
      cel2$time<-round(cel2$time, 1)
      
      
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$Cmic1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), Vars)
      colnames(mix)<-c("time", "r2", "E2" ,"Cmic2")
      #round the time to one digit
      mix$time<-round(mix$time, 1)
      
      
      mix2<-select(filter(data, Substrate!="Mix" ), Vars[1])
      colnames(mix2)<-"time"
      #round the time to one digit
      mix2$time<-round(mix2$time, 1)
      
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$Cmic2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-select(m1, ColM)
      
      
      cost_function<-function(pars){
        
        
        out<-mmem(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
        if(Mean==FALSE){
          cost<-modCost(model = out, obs = mtrue)
        }else{
          cost<-modCost(model = out, obs = mtrue, weight="mean")
        }
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(umax=0.1, Km=3, Ki=0.08, Ymic=0.8, Ye=0.8, alfa=2, m=0.001, kmic=0.01, ke=0.01),
                   lower=c(umax=0.001, Km=0.03, Ki=0.0008, Ymic=0.008, Ye=0.008, alfa=0.02, m=0.00001, kmic=0.0001, ke=0.0001),
                   upper=c(umax=100, Km=50, Ki=8, Ymic=1, Ye=1, alfa=20,  m=10, kmic=10, ke=10),niter=Niter)
      
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
    parameters<-data.frame(umax=vector("numeric", length = length(unique(dat$id))), 
                           Km=vector("numeric", length = length(unique(dat$id))), 
                           Ki=vector("numeric", length = length(unique(dat$id))),
                           Ymic=vector("numeric", length = length(unique(dat$id))),
                           Ye=vector("numeric", length = length(unique(dat$id))),
                           alfa=vector("numeric", length = length(unique(dat$id))),
                           #beta=vector("numeric", length = length(unique(dat$id))),
                           m=vector("numeric", length = length(unique(dat$id))),
                           kmic=vector("numeric", length = length(unique(dat$id))),
                           ke=vector("numeric", length = length(unique(dat$id))),
                           umax.l=vector("numeric", length = length(unique(dat$id))), 
                           Km.l=vector("numeric", length = length(unique(dat$id))), 
                           Ki.l=vector("numeric", length = length(unique(dat$id))),
                           Ymic.l=vector("numeric", length = length(unique(dat$id))),
                           Ye.l=vector("numeric", length = length(unique(dat$id))),
                           alfa.l=vector("numeric", length = length(unique(dat$id))),
                           #beta.l=vector("numeric", length = length(unique(dat$id))),
                           m.l=vector("numeric", length = length(unique(dat$id))),
                           kmic.l=vector("numeric", length = length(unique(dat$id))),
                           ke.l=vector("numeric", length = length(unique(dat$id))),
                           umax.u=vector("numeric", length = length(unique(dat$id))), 
                           Km.u=vector("numeric", length = length(unique(dat$id))), 
                           Ki.u=vector("numeric", length = length(unique(dat$id))),
                           Ymic.u=vector("numeric", length = length(unique(dat$id))),
                           Ye.u=vector("numeric", length = length(unique(dat$id))),
                           alfa.u=vector("numeric", length = length(unique(dat$id))),
                           #beta.u=vector("numeric", length = length(unique(dat$id))),
                           m.u=vector("numeric", length = length(unique(dat$id))),
                           kmic.u=vector("numeric", length = length(unique(dat$id))),
                           ke.u=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    
    #median
    for(i in unique(dat$id)){
      for(n in 1:9){
        
        parameters[i, n]<-summary(res[[i]])[6,n]
      }
    }
    #lower quartile
    for(i in unique(dat$id)){
      for(n in 1:9){
        
        parameters[i, n+9]<-summary(res[[i]])[5,n]
      }
    }
    #lower quartile
    for(i in unique(dat$id)){
      for(n in 1:9){
        
        parameters[i, n+18]<-summary(res[[i]])[7,n]
      }
    }
    
    if(FACT==2){
      
      parameters$Structure<-ddply(dat, .(id,Structure), summarize, mean(r))[,2]
      
    }else{
      
      
    }
    
    #OvP for respiration
    obs_r<-numeric()
    mod_r<-numeric()
    
    cost_r<-function(pars, data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), Vars)
      colnames(gl)<-c("time", "r", "E" ,"Cmic")
      #round the time to one digit
      gl$time<-round(gl$time, 1)
      
      gl2<-select(filter(data, Substrate!="Glucose" ), Vars[1])
      colnames(gl2)<-"time"
      #round the time to one digit
      gl2$time<-round(gl2$time, 1)
      
      gl2$r<-NA
      gl2$E<-NA
      gl2$Cmic<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), Vars)
      colnames(cel)<-c("time", "r1", "E1" ,"Cmic1")
      #round the time to one digit
      cel$time<-round(cel$time, 1)
      
      cel2<-select(filter(data, Substrate!="Cellobiose" ), Vars[1])
      colnames(cel2)<-"time"
      #round the time to one digit
      cel2$time<-round(cel2$time, 1)
      
      
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$Cmic1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), Vars)
      colnames(mix)<-c("time", "r2", "E2" ,"Cmic2")
      #round the time to one digit
      mix$time<-round(mix$time, 1)
      
      
      mix2<-select(filter(data, Substrate!="Mix" ), Vars[1])
      colnames(mix2)<-"time"
      #round the time to one digit
      mix2$time<-round(mix2$time, 1)
      
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$Cmic2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-select(m1, c("time", "r", "r1", "r2"))
      
      out<-mmem(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = mtrue)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      
      obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
      mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
      
    }
    
    OvP_r<-data.frame(obs_r, mod_r)
    
    #logLik calculation
    mu_r<-mean(obs_r)
    variance_r<-sd(obs_r)^2
    
    
    ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
    
    
    SSmodel_r<-sum((obs_r-mod_r)^2)
    SSdata_r<-sum((obs_r-mean(obs_r))^2)
    
    rsq_r=1-(SSmodel_r/SSdata_r)
    
    likelihood_r<-c(logLik=ll_r, npar=9*length(unique(dat$id)), rsq=rsq_r, AIC=2*9*length(unique(dat$id))-2*ll_r)
    
    #OvP for biomass
    obs_Cmic<-numeric()
    mod_Cmic<-numeric()
    
    cost_Cmic<-function(pars, data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), Vars)
      colnames(gl)<-c("time", "r", "E" ,"Cmic")
      #round the time to one digit
      gl$time<-round(gl$time, 1)
      
      gl2<-select(filter(data, Substrate!="Glucose" ), Vars[1])
      colnames(gl2)<-"time"
      #round the time to one digit
      gl2$time<-round(gl2$time, 1)
      
      gl2$r<-NA
      gl2$E<-NA
      gl2$Cmic<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), Vars)
      colnames(cel)<-c("time", "r1", "E1" ,"Cmic1")
      #round the time to one digit
      cel$time<-round(cel$time, 1)
      
      cel2<-select(filter(data, Substrate!="Cellobiose" ), Vars[1])
      colnames(cel2)<-"time"
      #round the time to one digit
      cel2$time<-round(cel2$time, 1)
      
      
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$Cmic1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), Vars)
      colnames(mix)<-c("time", "r2", "E2" ,"Cmic2")
      #round the time to one digit
      mix$time<-round(mix$time, 1)
      
      
      mix2<-select(filter(data, Substrate!="Mix" ), Vars[1])
      colnames(mix2)<-"time"
      #round the time to one digit
      mix2$time<-round(mix2$time, 1)
      
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$Cmic2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-select(m1, c("time", "Cmic", "Cmic1", "Cmic2"))
      
      out<-mmem(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = mtrue)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
      mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
    
    #logLik calculation
    mu_Cmic<-mean(obs_Cmic)
    variance_Cmic<-sd(obs_Cmic)^2
    
    
    ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
    
    
    SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
    SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
    
    rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
    
    likelihood_Cmic<-c(logLik=ll_Cmic, npar=9*length(unique(dat$id)), rsq=rsq_Cmic, AIC=2*9*length(unique(dat$id))-2*ll_Cmic)
    
    #OvP for enzymes
    obs_E<-numeric()
    mod_E<-numeric()
    
    cost_E<-function(pars, data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), Vars)
      colnames(gl)<-c("time", "r", "E" ,"Cmic")
      #round the time to one digit
      gl$time<-round(gl$time, 1)
      
      gl2<-select(filter(data, Substrate!="Glucose" ), Vars[1])
      colnames(gl2)<-"time"
      #round the time to one digit
      gl2$time<-round(gl2$time, 1)
      
      gl2$r<-NA
      gl2$E<-NA
      gl2$Cmic<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), Vars)
      colnames(cel)<-c("time", "r1", "E1" ,"Cmic1")
      #round the time to one digit
      cel$time<-round(cel$time, 1)
      
      cel2<-select(filter(data, Substrate!="Cellobiose" ), Vars[1])
      colnames(cel2)<-"time"
      #round the time to one digit
      cel2$time<-round(cel2$time, 1)
      
      
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$Cmic1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), Vars)
      colnames(mix)<-c("time", "r2", "E2" ,"Cmic2")
      #round the time to one digit
      mix$time<-round(mix$time, 1)
      
      
      mix2<-select(filter(data, Substrate!="Mix" ), Vars[1])
      colnames(mix2)<-"time"
      #round the time to one digit
      mix2$time<-round(mix2$time, 1)
      
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$Cmic2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-select(m1, c("time", "E", "E1", "E2"))
      
      out<-mmem(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
      mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_E<-data.frame(obs_E, mod_E)
    
    #logLik calculation
    mu_E<-mean(obs_E)
    variance_E<-sd(obs_E)^2
    
    
    ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
    
    
    SSmodel_E<-sum((obs_E-mod_E)^2)
    SSdata_E<-sum((obs_E-mean(obs_E))^2)
    
    rsq_E=1-(SSmodel_E/SSdata_E)
    
    likelihood_E<-c(logLik=ll_E, npar=9*length(unique(dat$id)), rsq=rsq_E, AIC=2*9*length(unique(dat$id))-2*ll_E)
    
    al<-list()
    
    al$parameters<-parameters
    al$OvP_r<-OvP_r
    al$ll<-rbind(likelihood_r, likelihood_Cmic, likelihood_E)
    al$OvP_Cmic<-OvP_Cmic
    al$OvP_E<-OvP_E
    
    
    return(al)
    
    
  }else{
    
    
    #mmem model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        #specific growth rate
        
        u=umax*C/(Km+C+Ki*C^2)
        
        #respiration rate
        r=m*Cmic+(1/Ymic-1)*Cmic*u+(1/Ye-1)*(Cmic*u*alfa)
        
        
        dCmic<-Cmic*u-kmic*Cmic
        dC<--(1/Ymic)*Cmic*u-(1/Ye)*(Cmic*u*alfa)-m*Cmic+kmic*Cmic+ke*E
        dE<-Cmic*u*alfa-ke*E
        
        return(list(c(dCmic, dC, dE), r=r))
        
      })
    }
    
    
    #create cost function
    estim<-function(data){
      
      Obs_dat<-select(data, Vars)
      colnames(Obs_dat)<-c("time", "r","E", "Cmic")
      #rounding the time
      Obs_dat$time<-round(Obs_dat$time,1)
      
      mtrue<-select(Obs_dat, ColM)  
      
      cinit<-as.numeric(data[1, "DOCinit"])
      
      cost_function<-function(pars){
        
        
        out<-as.data.frame(ode(y=c(Cmic=Cmic, C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
        if(Mean==FALSE){
          cost<-modCost(model = out, obs = mtrue)
        }else{
          cost<-modCost(model = out, obs = mtrue, weight="mean")
        }
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(umax=0.1, Km=3, Ki=0.08, Ymic=0.8, Ye=0.8, alfa=2,  m=0.001, kmic=0.01, ke=0.01),
                   lower=c(umax=0.001, Km=0.03, Ki=0.0008, Ymic=0.008, Ye=0.008, alfa=0.02,  m=0.00001, kmic=0.0001, ke=0.0001),
                   upper=c(umax=100, Km=50, Ki=8, Ymic=1, Ye=1, alfa=20, m=10, kmic=10, ke=10),niter=Niter)
      
      return(res)
      
    }
    
    #parameter estimation
    res<-foreach(i=unique(dat$id), .combine=list, .multicombine = TRUE,
                 .packages=c("FME", "dplyr")) %dopar% {
                   
                   estim(data=dat[dat$id==i,])
                   
                 }
    
    
    
    
    #parameters extraction
    parameters<-data.frame(umax=vector("numeric", length = length(unique(dat$id))), 
                           Km=vector("numeric", length = length(unique(dat$id))), 
                           Ki=vector("numeric", length = length(unique(dat$id))),
                           Ymic=vector("numeric", length = length(unique(dat$id))),
                           Ye=vector("numeric", length = length(unique(dat$id))),
                           alfa=vector("numeric", length = length(unique(dat$id))),
                           beta=vector("numeric", length = length(unique(dat$id))),
                           m=vector("numeric", length = length(unique(dat$id))),
                           kmic=vector("numeric", length = length(unique(dat$id))),
                           ke=vector("numeric", length = length(unique(dat$id))),
                           umax.l=vector("numeric", length = length(unique(dat$id))), 
                           Km.l=vector("numeric", length = length(unique(dat$id))), 
                           Ki.l=vector("numeric", length = length(unique(dat$id))),
                           Ymic.l=vector("numeric", length = length(unique(dat$id))),
                           Ye.l=vector("numeric", length = length(unique(dat$id))),
                           alfa.l=vector("numeric", length = length(unique(dat$id))),
                           #beta.l=vector("numeric", length = length(unique(dat$id))),
                           m.l=vector("numeric", length = length(unique(dat$id))),
                           kmic.l=vector("numeric", length = length(unique(dat$id))),
                           ke.l=vector("numeric", length = length(unique(dat$id))),
                           umax.u=vector("numeric", length = length(unique(dat$id))), 
                           Km.u=vector("numeric", length = length(unique(dat$id))), 
                           Ki.u=vector("numeric", length = length(unique(dat$id))),
                           Ymic.u=vector("numeric", length = length(unique(dat$id))),
                           Ye.u=vector("numeric", length = length(unique(dat$id))),
                           alfa.u=vector("numeric", length = length(unique(dat$id))),
                           #beta.u=vector("numeric", length = length(unique(dat$id))),
                           m.u=vector("numeric", length = length(unique(dat$id))),
                           kmic.u=vector("numeric", length = length(unique(dat$id))),
                           ke.u=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    
    #median
    for(i in unique(dat$id)){
      for(n in 1:9){
        
        parameters[i, n]<-summary(res[[i]])[6,n]
      }
    }
    #lower quartile
    for(i in unique(dat$id)){
      for(n in 1:9){
        
        parameters[i, n+9]<-summary(res[[i]])[5,n]
      }
    }
    #lower quartile
    for(i in unique(dat$id)){
      for(n in 1:9){
        
        parameters[i, n+18]<-summary(res[[i]])[7,n]
      }
    }
    
    if(FACT==1){
      
      parameters$Substrate<-ddply(dat, .(id,Substrate), summarize, mean(r))[,2]
      
    }else{
      
      parameters$Substrate<-ddply(dat, .(id,Substrate, Structure), summarize, mean(r))[,2]
      parameters$Structure<-ddply(dat, .(id,Substrate, Structure), summarize, mean(r))[,3]
    }
    
    #parameters[i, "Vmax"]<-summary(res[[i]])[6,1]
    #parameters[i, "Km"]<-summary(res[[i]])[6,2]
    #parameters[i, "CUE"]<-summary(res[[i]])[6,3]
    #parameters[i, "k"]<-summary(res[[i]])[6,4]
    
    #parameters[i, "Vmax.l"]<-summary(res[[i]])[5,1]
    #parameters[i, "Km.l"]<-summary(res[[i]])[5,2]
    #parameters[i, "CUE.l"]<-summary(res[[i]])[5,3]
    #parameters[i, "k.l"]<-summary(res[[i]])[5,4]
    
    #parameters[i, "Vmax.u"]<-summary(res[[i]])[7,1]
    #parameters[i, "Km.u"]<-summary(res[[i]])[7,2]
    #parameters[i, "CUE.u"]<-summary(res[[i]])[7,3]
    #parameters[i, "k.u"]<-summary(res[[i]])[7,4]
    
    #OvP for respiration
    obs_r<-numeric()
    mod_r<-numeric()
    
    cost_r<-function(pars, data){
      
      Obs_dat<-select(data, Vars[1:2])
      
      colnames(Obs_dat)<-c("time", "r")
      
      #rounding the time
      Obs_dat$time<-round(Obs_dat$time,1)
      
      cinit<-as.numeric(data[1,"DOCinit"])
      
      out<-out<-as.data.frame(ode(y=c(Cmic=Cmic, C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      
      obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
      mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
      
    }
    
    OvP_r<-data.frame(obs_r, mod_r)
    
    #logLik calculation
    mu_r<-mean(obs_r)
    variance_r<-sd(obs_r)^2
    
    
    ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
    
    
    SSmodel_r<-sum((obs_r-mod_r)^2)
    SSdata_r<-sum((obs_r-mean(obs_r))^2)
    
    rsq_r=1-(SSmodel_r/SSdata_r)
    
    likelihood_r<-c(logLik=ll_r, npar=9*length(unique(dat$id)), rsq=rsq_r, AIC=2*9*length(unique(dat$id))-2*ll_r)
    
    #OvP for biomass
    obs_Cmic<-numeric()
    mod_Cmic<-numeric()
    
    cost_Cmic<-function(pars, data){
      
      
      Obs_decay<-select(data, Vars[c(1,4)])
      
      colnames(Obs_dat)<-c("time", "Cmic")
      #rounding the time
      Obs_dat$time<-round(Obs_dat$time,1)
      
      cinit<-as.numeric(data[1, "DOCinit"])
      
      out<-as.data.frame(ode(y=c(Cmic=Cmic, C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
      mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
    
    #logLik calculation
    mu_Cmic<-mean(obs_Cmic)
    variance_Cmic<-sd(obs_Cmic)^2
    
    
    ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
    
    
    SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
    SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
    
    rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
    
    likelihood_Cmic<-c(logLik=ll_Cmic, npar=9*length(unique(dat$id)), rsq=rsq_Cmic, AIC=2*9*length(unique(dat$id))-2*ll_Cmic)
    
    #OvP for enzymes
    obs_E<-numeric()
    mod_E<-numeric()
    
    cost_E<-function(pars, data){
      
      
      Obs_decay<-select(data, Vars[c(1,3)])
      
      colnames(Obs_dat)<-c("time", "E")
      #rounding the time
      Obs_dat$time<-round(Obs_dat$time,1)
      
      Ci<-as.numeric(data[1, "DOCinit"])
      
      out<-ode(y=c(Cmic=Cmic, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
      mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_E<-data.frame(obs_E, mod_E)
    
    #logLik calculation
    mu_E<-mean(obs_E)
    variance_E<-sd(obs_E)^2
    
    
    ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
    
    
    SSmodel_E<-sum((obs_E-mod_E)^2)
    SSdata_E<-sum((obs_E-mean(obs_E))^2)
    
    rsq_E=1-(SSmodel_E/SSdata_E)
    
    likelihood_E<-c(logLik=ll_E, npar=9*length(unique(dat$id)), rsq=rsq_E, AIC=2*9*length(unique(dat$id))-2*ll_E)
    
    
    
    al<-list()
    
    al$parameters<-parameters
    al$OvP_r<-OvP_r
    al$ll<-rbind(likelihood_r, likelihood_Cmic, likelihood_E)
    al$OvP_Cmic<-OvP_Cmic
    al$OvP_E<-OvP_E
    
    return(al)
    
  }
  
}