DB_constantMB<-function(data){
  #Define the model with constant parameters
  DB<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      #Define equations
      ##Enzymatic degradation - multiplicative equation
      Deg<-Enz * Cs * Ecat
      ##Carbon uptake
      Cu<-0.450799793*S*DOC/(17.77449941 + DOC)
      ##maintnance
      m=S*0.004769861
      ##Reserves available
      an=3.188932004*R-m
      ##Respiration rate
      r=pmax(chi*an*(1-0.6), 0) + pmax((1-chi)*an*(1-Yue), 0) + m
      
      ##Proteins and DNA
      Protinc=fpr*R+fps*S
      DNAc=fds*S
      
      #States
      dCs <- -Deg
      dDOC <- Deg - Cu
      dR <- Cu - 3.188932004*R
      dS <- pmax(0, chi*an*0.6, 0) + pmin(0, an)
      dEnz <- pmax(0, (1-chi)*an*Yue, 0)
      
      return(list(c(dCs, dDOC, dR, dS, dEnz), 
                  r=r, Protinc = Protinc, 
                  DNAc = DNAc))
    })
  }
  
  #Define parameters for constant model
  parnames=c("Ecat", "chi", "Yue", "fpr", "fps", "fds","RSinit")
  
  #Define estimation function
  estim<-function(odeset){
    #Cost function
    cost<-function(x){
      p<-x
      #p<-out_const$pars
      names(p)<-parnames
      
      #Initial conditions
      Sinit<-mean(odeset$DNA.initc, na.rm=T)/p[["fds"]]
      #Sinit<-mean(m0CB$DNA.initc, na.rm=T)/p[["fds"]]
      Rinit<-p[["RSinit"]]*Sinit
      
      #Simulation
      yhat_all <- as.data.frame(ode(y=c(Cs=25/3, DOC=25/3*2, R=Rinit, S=Sinit, Enz=0), func = DB, times = seq(0,120),
                  parms = p[-which(names(p)=="RSinit")]))
      
      #Filter the output
      yhat <- yhat_all[, c("Protinc", "DNAc", "Enz", "r", "time")]
      colnames(yhat)<-c("Protinc", "DNAc", "Protoutc", "r", "time")
      Yhat <- melt(yhat, id.vars=c("time"))
      colnames(Yhat)<-c("time", "variable", "yhat")
      
      #Observations
      obs <- odeset[, c("Protinc", "DNAc", "Protoutc", "r", "Time2")]
      Obs <- melt(obs, id.vars=c("Time2"))
      colnames(Obs)<-c("time", "variable", "obs")
      
      #Merge both
      Yhat <- merge(Yhat, Obs, by = c("time", "variable"))
      
      #Add weighting factor
      yweights<-Yhat %>% group_by(variable) %>% summarise(weights=mean(obs, na.rm = T))
      
      #match with the Yhat data frame
      Yhat<-merge(Yhat, yweights, by=c("variable"))
      
      #now, the root mean square error is calculated
      NRMSE<-as.numeric(Yhat %>% group_by(variable) %>% summarise(NRMSE=sum((((yhat-obs)/mean(weights))^2), na.rm = T)) %>%
                          summarise(NRMSE=sum(NRMSE)))
      
      return(NRMSE)
    }
    
    #Goodness of fit function
    good<-function(x){
      p<-x
      names(p)<-parnames
      
      #Initial conditions
      Sinit<-mean(odeset$DNA.initc, na.rm=T)/p[["fds"]]
      Rinit<-p[["RSinit"]]*Sinit
      
      #Simulation
      yhat_all <- as.data.frame(ode(y=c(Cs=25/3, DOC=25/3*2, R=Rinit, S=Sinit, Enz=0), func = DB, times = seq(0,120),
                                    parms = p[-which(names(p)=="RSinit")]))
      
      #Filter the output
      yhat <- yhat_all[, c("Protinc", "DNAc", "Enz", "r", "time")]
      colnames(yhat)<-c("Protinc", "DNAc", "Protoutc", "r", "time")
      Yhat <- melt(yhat, id.vars=c("time"))
      colnames(Yhat)<-c("time", "variable", "yhat")
      
      #Observations
      obs <- odeset[, c("Protinc", "DNAc", "Protoutc", "r", "Time2")]
      Obs <- melt(obs, id.vars=c("Time2"))
      colnames(Obs)<-c("time", "variable", "obs")
      
      #Merge both
      Yhat <- merge(Yhat, Obs, by = c("time", "variable"))
      
      #rsquared calculation for each variable
      Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-yhat)^2), na.rm = T),
                                                      SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                      ll=-sum(((obs-yhat)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
      Gfit$R2<-with(Gfit, 1-SSres/SStot)
      Gfit$N<-length(p)
      Gfit$AIC<-with(Gfit, 2*N-2*ll)
      
      rsq_out<-list(Yhat=Yhat, Gfit=Gfit)
      
      return(rsq_out)
    }
    
    #Estimate the parameters
    #approximate parameter estimation is done by MCMC method
    par_mcmc<-modMCMC(f=cost, p=c(Ecat=0.01, chi=0.5, Yue=0.6, fpr=0.8, 
                                  fps=0.4, fds=0.004, RSinit=0.1),
                      lower=c(Ecat=1e-5, chi=0, Yue=0, fpr=0, 
                              fps=0, fds=0, RSinit=0),
                      upper=c(Ecat=10, chi=1, Yue=0.9, fpr=1, 
                              fps=1, fds=1, RSinit=50), niter=10000)
    
    #lower and upper limits for parameters are extracted
    pl<-summary(par_mcmc)["min",]
    pu<-summary(par_mcmc)["max",]
    
    #these limits are used to find global optimum by DEoptim
    opt_par<-DEoptim(fn=cost, lower=pl, upper=pu,
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8,
                                 trace=FALSE, strategy=3, NP=250))
    
    #goodness of fit
    fit<-good(opt_par$optim$bestmem)
    
    #best parameters
    p<-opt_par$optim$bestmem
    names(p)<-parnames
    
    #return list with opt_par and par_prof
    estim_out<-list(pars=p, par_mcmc=par_mcmc, fit=fit)
    
    return(estim_out)
  }
  
  #Use the function
  res <- estim(data)
  
  return(res)
}