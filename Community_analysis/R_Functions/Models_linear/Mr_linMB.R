Mr_linMB<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    #Define time variable parameter
    Mr = I + B*ifelse(time<25, 1 + 0.66352*time, 17.6413*exp(-0.0025*time))
    
    #Define equations
    ##Enzymatic degradation - multiplicative equation
    Deg<-Enz * Cs * Ecat
    ##Carbon uptake
    Cu<-0.450799793*S*DOC/(17.77449941 + DOC)
    ##maintnance
    m=S*Mr
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