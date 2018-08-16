mend_init_a<-function(data, SUB, FACT, Mean, Niter, DNAci){
  
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
    
    #this is the mend model fitted across different substrates 
    #with different initial carbon concentrations
    mend<-function(X, pars, t){
      
      #mend model
      deriv<-function(time, state, pars){
        
        with(as.list(c(state, pars)),{
          
          #equations
          #C uptake
          F1=1/CUE*(Vmax+mr)*Cmic*C/(Km+C)
          #growth respiration
          F4=(1/CUE-1)*Vmax*Cmic*C/(Km+C)
          #maintnance respiration
          F5=(1/CUE-1)*mr*Cmic*C/(Km+C)
          #microbial mortality
          F8=(1-pe)*mr*Cmic
          #enzyme production
          F9=pe*mr*Cmic
          #enzyme decay
          F10=ke*E
          
          #states
          dCmic<-F1-(F4+F5)-F8-F9
          dC<--F1+F8+F10
          dE<-F9-F10
          
          #DNA and Protein content in biomass
          DNAc=fd*Cmic
          #Protc=fp*Cmic
          
          return(list(c(dCmic, dC, dE), r=F4+F5, DNAc=DNAc))
          
        })
      }
      
      #function to apply acros 3 different initial substrate concentrations
      base_function<-function(Ci){
        
        
        out<-ode(y=c(Cmic=pars[["Cmic_0"]], C=Ci, E=0), parms=pars, times=t, func=deriv)
        return(as.data.frame(out))}
      
      #results
      res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
      
      return(res)
      
    }
    
    #create cost function
    estim<-function(data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(gl)<-c("time", "r", "E" ,"DNAc")
      gl2<-select(filter(data, Substrate!="Glucose" ), "Time")
      colnames(gl2)<-"time"
      gl2$r<-NA
      gl2$E<-NA
      gl2$DNAc<-NA
      #gl2$Protc<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(cel)<-c("time", "r1", "E1" ,"DNAc1")
      cel2<-select(filter(data, Substrate!="Cellobiose" ), "Time")
      colnames(cel2)<-"time"
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$DNAc1<-NA
      #cel2$Protc1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), c("Time", "r", "E" ,"DNAc"))
      colnames(mix)<-c("time", "r2", "E2" ,"DNAc2")
      mix2<-select(filter(data, Substrate!="Mix" ), "Time")
      colnames(mix2)<-"time"
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$DNAc2<-NA
      #mix2$Protc2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-rbind(m1, data.frame(time=0, r=NA, E=0, DNAc=DNAci,
                                  r1=NA, E1=0, DNAc1=DNAci, 
                                  r2=NA, E2=0, DNAc2=DNAci))
      
      
      cost_function<-function(pars){
        
        
        out<-mend(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
        
        if(Mean==FALSE){
        cost<-modCost(model = out, obs = mtrue)
        }else{
          cost<-modCost(model = out, obs = mtrue, weight="mean")
        }
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, mr=0.01, ke=0.01, pe=0.5, fd=0.05, Cmic_0=0.12),
                   lower=c(Vmax=1e-4, Km=1e-4, CUE=0, mr=1e-5, ke=1e-6, pe=0, fd=0,  Cmic_0=1e-6),
                   upper=c(Vmax=1e4, Km=1e4, CUE=1, mr=1e5, ke=1e6, pe=1, fd=1, Cmic_0=25),niter=Niter)
      
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
                           mr=vector("numeric", length = length(unique(dat$id))),
                           ke=vector("numeric", length = length(unique(dat$id))),
                           pe=vector("numeric", length = length(unique(dat$id))),
                           fd=vector("numeric", length = length(unique(dat$id))),
                           #fp=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0=vector("numeric", length = length(unique(dat$id))),
                           Vmax.l=vector("numeric", length = length(unique(dat$id))), 
                           Km.l=vector("numeric", length = length(unique(dat$id))), 
                           CUE.l=vector("numeric", length = length(unique(dat$id))), 
                           mr.l=vector("numeric", length = length(unique(dat$id))),
                           ke.l=vector("numeric", length = length(unique(dat$id))),
                           pe.l=vector("numeric", length = length(unique(dat$id))),
                           fd.l=vector("numeric", length = length(unique(dat$id))),
                           #fp.l=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0.l=vector("numeric", length = length(unique(dat$id))),
                           Vmax.u=vector("numeric", length = length(unique(dat$id))), 
                           Km.u=vector("numeric", length = length(unique(dat$id))), 
                           CUE.u=vector("numeric", length = length(unique(dat$id))), 
                           mr.u=vector("numeric", length = length(unique(dat$id))),
                           ke.u=vector("numeric", length = length(unique(dat$id))),
                           pe.u=vector("numeric", length = length(unique(dat$id))),
                           fd.u=vector("numeric", length = length(unique(dat$id))),
                           #fp.u=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0.u=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    
    #median
    for(i in unique(dat$id)){
      for(n in 1:((ncol(parameters)-1)/3)){
        
        parameters[i, n]<-summary(res[[i]])[1,n]
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
    
    if(FACT==2){
      
      parameters$Structure<-ddply(dat, .(id,Structure), summarize, mean(r))[,2]
      
    }else{
      
      
    }
    
    #OvP for respiration
    obs_r<-numeric()
    mod_r<-numeric()
    
    cost_r<-function(pars, data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(gl)<-c("time", "r", "E" ,"DNAc")
      gl2<-select(filter(data, Substrate!="Glucose" ), "Time")
      colnames(gl2)<-"time"
      gl2$r<-NA
      gl2$E<-NA
      gl2$DNAc<-NA
      #gl2$Protc<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(cel)<-c("time", "r1", "E1" ,"DNAc1")
      cel2<-select(filter(data, Substrate!="Cellobiose" ), "Time")
      colnames(cel2)<-"time"
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$DNAc1<-NA
      #cel2$Protc1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), c("Time", "r", "E" ,"DNAc"))
      colnames(mix)<-c("time", "r2", "E2" ,"DNAc2")
      mix2<-select(filter(data, Substrate!="Mix" ), "Time")
      colnames(mix2)<-"time"
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$DNAc2<-NA
      #mix2$Protc2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-rbind(m1, data.frame(time=0, r=NA, E=0, DNAc=DNAci,
                                  r1=NA, E1=0, DNAc1=DNAci, 
                                  r2=NA, E2=0, DNAc2=DNAci))
      
      mtrue<-select(mtrue, c("time", "r", "r1", "r2"))
      
      out<-mend(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = mtrue)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      
      obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
      mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
      
    }
    
    OvP_r<-data.frame(obs_r, mod_r)
    
    #logLik calculation
    mu_r<-mean(obs_r)
    variance_r<-sd(obs_r)^2
    
    
    ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
    
    
    SSmodel_r<-sum((obs_r-mod_r)^2)
    SSdata_r<-sum((obs_r-mean(obs_r))^2)
    
    rsq_r=1-(SSmodel_r/SSdata_r)
    
    likelihood_r<-c(logLik=ll_r, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_r, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_r)
    
    #OvP for DNA
    obs_DNAc<-numeric()
    mod_DNAc<-numeric()
    
    cost_DNAc<-function(pars, data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(gl)<-c("time", "r", "E" ,"DNAc")
      gl2<-select(filter(data, Substrate!="Glucose" ), "Time")
      colnames(gl2)<-"time"
      gl2$r<-NA
      gl2$E<-NA
      gl2$DNAc<-NA
      #gl2$Protc<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(cel)<-c("time", "r1", "E1" ,"DNAc1")
      cel2<-select(filter(data, Substrate!="Cellobiose" ), "Time")
      colnames(cel2)<-"time"
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$DNAc1<-NA
      #cel2$Protc1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), c("Time", "r", "E" ,"DNAc"))
      colnames(mix)<-c("time", "r2", "E2" ,"DNAc2")
      mix2<-select(filter(data, Substrate!="Mix" ), "Time")
      colnames(mix2)<-"time"
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$DNAc2<-NA
      #mix2$Protc2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-rbind(m1, data.frame(time=0, r=NA, E=0, DNAc=DNAci,
                                  r1=NA, E1=0, DNAc1=DNAci, 
                                  r2=NA, E2=0, DNAc2=DNAci))
      
      
      mtrue<-select(mtrue, c("time", "DNAc", "DNAc1", "DNAc2"))
      
      out<-mend(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = mtrue)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_DNAc<-append(obs_DNAc, cost_DNAc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
      mod_DNAc<-append(mod_DNAc, cost_DNAc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_DNAc<-data.frame(obs_DNAc, mod_DNAc)
    
    #logLik calculation
    mu_DNAc<-mean(obs_DNAc)
    variance_DNAc<-sd(obs_DNAc)^2
    
    
    ll_DNAc<--1*sum((obs_DNAc-mod_DNAc)^2)/2/variance_DNAc
    
    
    SSmodel_DNAc<-sum((obs_DNAc-mod_DNAc)^2)
    SSdata_DNAc<-sum((obs_DNAc-mean(obs_DNAc))^2)
    
    rsq_DNAc=1-(SSmodel_DNAc/SSdata_DNAc)
    
    likelihood_DNAc<-c(logLik=ll_DNAc, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), 
                       rsq=rsq_DNAc, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_DNAc)
    
    # #OvP for Protein
    # obs_Protc<-numeric()
    # mod_Protc<-numeric()
    # 
    # cost_Protc<-function(pars, data){
    #   
    #   ###Glucose
    #   gl<-select(filter(data, Substrate=="Glucose" ), c("Time", "r", "E" ,"DNAc", "Protc"))
    #   colnames(gl)<-c("time", "r", "E" ,"DNAc", "Protc")
    #   gl2<-select(filter(data, Substrate!="Glucose" ), "Time")
    #   colnames(gl2)<-"time"
    #   gl2$r<-NA
    #   gl2$E<-NA
    #   gl2$DNAc<-NA
    #   gl2$Protc<-NA
    #   m1<-rbind(gl, gl2)
    #   m1<-m1[order(m1$time),]
    #   
    #   ###Cellobiose
    #   cel<-select(filter(data, Substrate=="Cellobiose" ), c("Time", "r", "E" ,"DNAc", "Protc"))
    #   colnames(cel)<-c("time", "r1", "E1" ,"DNAc1", "Protc1")
    #   cel2<-select(filter(data, Substrate!="Cellobiose" ), "Time")
    #   colnames(cel2)<-"time"
    #   cel2$r1<-NA
    #   cel2$E1<-NA
    #   cel2$DNAc1<-NA
    #   cel2$Protc1<-NA
    #   m2<-rbind(cel, cel2)
    #   m2<-m2[order(m2$time),]
    #   
    #   ###Mix
    #   mix<-select(filter(data, Substrate=="Mix" ), c("Time", "r", "E" ,"DNAc", "Protc"))
    #   colnames(mix)<-c("time", "r2", "E2" ,"DNAc2", "Protc2")
    #   mix2<-select(filter(data, Substrate!="Mix" ), "Time")
    #   colnames(mix2)<-"time"
    #   mix2$r2<-NA
    #   mix2$E2<-NA
    #   mix2$DNAc2<-NA
    #   mix2$Protc2<-NA
    #   m3<-rbind(mix, mix2)
    #   m3<-m3[order(m3$time),]
    #   
    #   m1[,c(6:9)]<-m2[,c(2:5)]
    #   m1[,c(10:13)]<-m3[,c(2:5)]
    #   
    #   mtrue<-rbind(m1, data.frame(time=0, r=NA, E=0, DNAc=DNAci, Protc=NA,
    #                               r1=NA, E1=0, DNAc1=DNAci, Protc1=NA,
    #                               r2=NA, E2=0, DNAc2=DNAci, Protc2=NA))
    #   
    #   
    #   mtrue<-select(mtrue, c("time", "Protc", "Protc1", "Protc2"))
    #   
    #   out<-mend(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
    #   cost<-modCost(model = out, obs = mtrue)
    #   
    #   return(cost)
    #   
    # }
    # 
    # for(i in unique(dat$id)){
    #   obs_Protc<-append(obs_Protc, cost_Protc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
    #   mod_Protc<-append(mod_Protc, cost_Protc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
    # }
    # 
    # OvP_Protc<-data.frame(obs_Protc, mod_Protc)
    # 
    # #logLik calculation
    # mu_Protc<-mean(obs_Protc)
    # variance_Protc<-sd(obs_Protc)^2
    # 
    # 
    # ll_Protc<--1*sum((obs_Protc-mod_Protc)^2)/2/variance_Protc
    # 
    # 
    # SSmodel_Protc<-sum((obs_Protc-mod_Protc)^2)
    # SSdata_Protc<-sum((obs_Protc-mean(obs_Protc))^2)
    # 
    # rsq_Protc=1-(SSmodel_Protc/SSdata_Protc)
    # 
    # likelihood_Protc<-c(logLik=ll_Protc, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), 
    #                    rsq=rsq_Protc, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_Protc)
    # 
    
    #OvP for enzymes
    obs_E<-numeric()
    mod_E<-numeric()
    
    cost_E<-function(pars, data){
      
      ###Glucose
      gl<-select(filter(data, Substrate=="Glucose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(gl)<-c("time", "r", "E" ,"DNAc")
      gl2<-select(filter(data, Substrate!="Glucose" ), "Time")
      colnames(gl2)<-"time"
      gl2$r<-NA
      gl2$E<-NA
      gl2$DNAc<-NA
      #gl2$Protc<-NA
      m1<-rbind(gl, gl2)
      m1<-m1[order(m1$time),]
      
      ###Cellobiose
      cel<-select(filter(data, Substrate=="Cellobiose" ), c("Time", "r", "E" ,"DNAc"))
      colnames(cel)<-c("time", "r1", "E1" ,"DNAc1")
      cel2<-select(filter(data, Substrate!="Cellobiose" ), "Time")
      colnames(cel2)<-"time"
      cel2$r1<-NA
      cel2$E1<-NA
      cel2$DNAc1<-NA
      #cel2$Protc1<-NA
      m2<-rbind(cel, cel2)
      m2<-m2[order(m2$time),]
      
      ###Mix
      mix<-select(filter(data, Substrate=="Mix" ), c("Time", "r", "E" ,"DNAc"))
      colnames(mix)<-c("time", "r2", "E2" ,"DNAc2")
      mix2<-select(filter(data, Substrate!="Mix" ), "Time")
      colnames(mix2)<-"time"
      mix2$r2<-NA
      mix2$E2<-NA
      mix2$DNAc2<-NA
      #mix2$Protc2<-NA
      m3<-rbind(mix, mix2)
      m3<-m3[order(m3$time),]
      
      m1[,c(5:7)]<-m2[,c(2:4)]
      m1[,c(8:10)]<-m3[,c(2:4)]
      
      mtrue<-rbind(m1, data.frame(time=0, r=NA, E=0, DNAc=DNAci,
                                  r1=NA, E1=0, DNAc1=DNAci, 
                                  r2=NA, E2=0, DNAc2=DNAci))
      
      
      mtrue<-select(mtrue, c("time", "E", "E1", "E2"))
      
      out<-mend(X=c(25, 25, 16.5), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
      mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_E<-data.frame(obs_E, mod_E)
    
    #logLik calculation
    mu_E<-mean(obs_E)
    variance_E<-sd(obs_E)^2
    
    
    ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
    
    
    SSmodel_E<-sum((obs_E-mod_E)^2)
    SSdata_E<-sum((obs_E-mean(obs_E))^2)
    
    rsq_E=1-(SSmodel_E/SSdata_E)
    
    likelihood_E<-c(logLik=ll_E, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), rsq=rsq_E, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_E)
    
    al<-list()
    
    al$parameters<-parameters
    al$OvP_r<-OvP_r
    al$ll<-rbind(likelihood_r, likelihood_DNAc, likelihood_E)
    al$OvP_DNAc<-OvP_DNAc
    #al$OvP_Protc<-OvP_Protc
    al$OvP_E<-OvP_E
    
    return(al)
    
    
  }else{
    
    
    #mend model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        #equations
        #C uptake
        F1=1/CUE*(Vmax+mr)*Cmic*C/(Km+C)
        #growth respiration
        F4=(1/CUE-1)*Vmax*Cmic*C/(Km+C)
        #maintnance respiration
        F5=(1/CUE-1)*mr*Cmic*C/(Km+C)
        #microbial mortality
        F8=(1-pe)*mr*Cmic
        #enzyme production
        F9=pe*mr*Cmic
        #enzyme decay
        F10=ke*E
        
        #states
        dCmic<-F1-(F4+F5)-F8-F9
        dC<--F1+F8+F10
        dE<-F9-F10
        
        #DNA and Protein content in biomass
        DNAc=fd*Cmic
        #Protc=fp*Cmic
        
        return(list(c(dCmic, dC, dE), r=F4+F5, DNAc=DNAc))
        
      })
    }
    
    
    #create cost function
    estim<-function(data){
      
      Obs_dat<-select(data, c("Time", "r","E", "DNAc"))
      colnames(Obs_dat)<-c("time", "r","E", "DNAc")
      mtrue<-rbind(Obs_dat, data.frame(time=0, r=NA, E=NA, DNAc=DNAci))  
      
      cinit<-as.numeric(data[1, "DOCinit"])
      
      cost_function<-function(pars){
        
        
        out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
        
        if(Mean==FALSE){
          cost<-modCost(model = out, obs = mtrue)
        }else{
          cost<-modCost(model = out, obs = mtrue, weight="mean")
        }
        
        return(cost)
        
      }
      
      res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, mr=0.01, ke=0.01, pe=0.5, fd=0.05, Cmic_0=0.12),
                   lower=c(Vmax=1e-4, Km=1e-4, CUE=0, mr=1e-5, ke=1e-6, pe=0, fd=0, Cmic_0=1e-6),
                   upper=c(Vmax=1e4, Km=1e4, CUE=1, mr=1e5, ke=1e6, pe=1, fd=1, Cmic_0=25),niter=Niter)
      
      
      return(res)
      
    }
    
    #parameter estimation
    res<-foreach(i=unique(dat$id), .combine=list, .multicombine = TRUE,
                 .packages=c("FME", "dplyr")) %dopar% {
                   
                   estim(data=dat[dat$id==i,])
                   
                 }
    
    
    
    
    #parameters extraction
    #parameters extraction
    parameters<-data.frame(Vmax=vector("numeric", length = length(unique(dat$id))), 
                           Km=vector("numeric", length = length(unique(dat$id))), 
                           CUE=vector("numeric", length = length(unique(dat$id))), 
                           mr=vector("numeric", length = length(unique(dat$id))),
                           ke=vector("numeric", length = length(unique(dat$id))),
                           pe=vector("numeric", length = length(unique(dat$id))),
                           fd=vector("numeric", length = length(unique(dat$id))),
                           #fp=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0=vector("numeric", length = length(unique(dat$id))),
                           Vmax.l=vector("numeric", length = length(unique(dat$id))), 
                           Km.l=vector("numeric", length = length(unique(dat$id))), 
                           CUE.l=vector("numeric", length = length(unique(dat$id))), 
                           mr.l=vector("numeric", length = length(unique(dat$id))),
                           ke.l=vector("numeric", length = length(unique(dat$id))),
                           pe.l=vector("numeric", length = length(unique(dat$id))),
                           fd.l=vector("numeric", length = length(unique(dat$id))),
                           #fp.l=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0.l=vector("numeric", length = length(unique(dat$id))),
                           Vmax.u=vector("numeric", length = length(unique(dat$id))), 
                           Km.u=vector("numeric", length = length(unique(dat$id))), 
                           CUE.u=vector("numeric", length = length(unique(dat$id))), 
                           mr.u=vector("numeric", length = length(unique(dat$id))),
                           ke.u=vector("numeric", length = length(unique(dat$id))),
                           pe.u=vector("numeric", length = length(unique(dat$id))),
                           fd.u=vector("numeric", length = length(unique(dat$id))),
                           #fp.u=vector("numeric", length = length(unique(dat$id))),
                           Cmic_0.u=vector("numeric", length = length(unique(dat$id))),
                           ID=unique(dat$id))
    
    #median
    for(i in unique(dat$id)){
      for(n in 1:((ncol(parameters)-1)/3)){
        
        parameters[i, n]<-summary(res[[i]])[1,n]
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
      
      Obs_dat<-select(data, "Time", "r")
      
      colnames(Obs_dat)<-c("time", "r")
      cinit<-as.numeric(data[1,"DOCinit"])
      
      out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      
      obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
      mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
      
    }
    
    OvP_r<-data.frame(obs_r, mod_r)
    
    #logLik calculation
    mu_r<-mean(obs_r)
    variance_r<-sd(obs_r)^2
    
    
    ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
    
    
    SSmodel_r<-sum((obs_r-mod_r)^2)
    SSdata_r<-sum((obs_r-mean(obs_r))^2)
    
    rsq_r=1-(SSmodel_r/SSdata_r)
    
    likelihood_r<-c(logLik=ll_r, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), 
                    rsq=rsq_r, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_r)
    
    #OvP for DNA
    obs_DNAc<-numeric()
    mod_DNAc<-numeric()
    
    cost_DNAc<-function(pars, data){
      
      Obs_dat<-select(data, c("Time", "r","E", "DNAc", "Protc"))
      colnames(Obs_dat)<-c("time", "r","E", "DNAc", "Protc")
      mtrue<-rbind(Obs_dat, data.frame(time=0, r=NA, E=NA, DNAc=DNAci, Protc=NA))  
      
      mtrue<-select(mtrue, c("time", "DNAc"))
      
      out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = mtrue)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_DNAc<-append(obs_DNAc, cost_DNAc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
      mod_DNAc<-append(mod_DNAc, cost_DNAc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_DNAc<-data.frame(obs_DNAc, mod_DNAc)
    
    #logLik calculation
    mu_DNAc<-mean(obs_DNAc)
    variance_DNAc<-sd(obs_DNAc)^2
    
    
    ll_DNAc<--1*sum((obs_DNAc-mod_DNAc)^2)/2/variance_DNAc
    
    
    SSmodel_DNAc<-sum((obs_DNAc-mod_DNAc)^2)
    SSdata_DNAc<-sum((obs_DNAc-mean(obs_DNAc))^2)
    
    rsq_DNAc=1-(SSmodel_DNAc/SSdata_DNAc)
    
    likelihood_DNAc<-c(logLik=ll_DNAc, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), 
                       rsq=rsq_DNAc, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_DNAc)
    
    # #OvP for Prot
    # obs_Protc<-numeric()
    # mod_Protc<-numeric()
    # 
    # cost_Protc<-function(pars, data){
    #   
    #   Obs_dat<-select(data, c("Time", "r","E", "DNAc", "Protc"))
    #   colnames(Obs_dat)<-c("time", "r","E", "DNAc", "Protc")
    #   mtrue<-rbind(Obs_dat, data.frame(time=0, r=NA, E=NA, DNAc=DNAci, Protc=NA))  
    #   
    #   mtrue<-select(mtrue, c("time", "Protc"))
    #   
    #   out<-as.data.frame(ode(y=c(Cmic=pars[["Cmic_0"]], C=cinit, E=0), parms=pars, times=seq(0,130), func=deriv))
    #   cost<-modCost(model = out, obs = mtrue)
    #   
    #   return(cost)
    #   
    # }
    # 
    # for(i in unique(dat$id)){
    #   obs_Protc<-append(obs_Protc, cost_Protc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
    #   mod_Protc<-append(mod_Protc, cost_Protc(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
    # }
    # 
    # OvP_Protc<-data.frame(obs_Protc, mod_Protc)
    # 
    # #logLik calculation
    # mu_Protc<-mean(obs_Protc)
    # variance_Protc<-sd(obs_Protc)^2
    # 
    # 
    # ll_Protc<--1*sum((obs_Protc-mod_Protc)^2)/2/variance_Protc
    # 
    # 
    # SSmodel_Protc<-sum((obs_Protc-mod_Protc)^2)
    # SSdata_Protc<-sum((obs_Protc-mean(obs_Protc))^2)
    # 
    # rsq_Protc=1-(SSmodel_Protc/SSdata_Protc)
    # 
    # likelihood_Protc<-c(logLik=ll_Protc, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), 
    #                    rsq=rsq_Protc, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_Protc)
    # 
    #OvP for enzymes
    obs_E<-numeric()
    mod_E<-numeric()
    
    cost_E<-function(pars, data){
      
      Obs_decay<-select(data, c("Time", "E"))
      
      colnames(Obs_dat)<-c("time", "E")
      
      Ci<-as.numeric(data[1, "DOCinit"])
      
      out<-ode(y=c(Cmic=pars[["Cmic_0"]], C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
    
    for(i in unique(dat$id)){
      obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$obs)
      mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[1,], data=dat[dat$id==i, ])$residuals$mod)
    }
    
    OvP_E<-data.frame(obs_E, mod_E)
    
    #logLik calculation
    mu_E<-mean(obs_E)
    variance_E<-sd(obs_E)^2
    
    
    ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
    
    
    SSmodel_E<-sum((obs_E-mod_E)^2)
    SSdata_E<-sum((obs_E-mean(obs_E))^2)
    
    rsq_E=1-(SSmodel_E/SSdata_E)
    
    likelihood_E<-c(logLik=ll_E, npar=((ncol(parameters)-1)/3)*length(unique(dat$id)), 
                    rsq=rsq_E, AIC=2*((ncol(parameters)-1)/3)*length(unique(dat$id))-2*ll_E)
    
    
    
    al<-list()
    
    al$parameters<-parameters
    al$OvP_r<-OvP_r
    al$ll<-rbind(likelihood_r, likelihood_DNAc, likelihood_E)
    al$OvP_DNAc<-OvP_DNAc
    #al$OvP_Protc<-OvP_Protc
    al$OvP_E<-OvP_E
    
    return(al)
    
  }
  
}