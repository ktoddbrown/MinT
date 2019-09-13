fps_lin<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    #Define time variable parameter
    fps = I + B*(4.221*x/(30.700+x+x^2/106.021)+1)
    
    #Define equations
    ##Enzymatic degradation - multiplicative equation
    Deg<-Enz * Cs * Ecat
    ##Carbon uptake
    Cu<-0.283475617*S*DOC/(2.512518101 + DOC)
    ##maintnance
    m=S*Mr
    ##Reserves available
    an=f*R-m
    ##Respiration rate
    r=pmin(m, f*R)+pmax(chi*an*(1-0.6), 0) + pmax((1-chi)*an*(1-Yue), 0)
    
    ##Proteins and DNA
    Protinc=fpr*R+fps*S
    DNAc=fds*S
    
    #States
    dCs <- -Deg + pmax(0, -j*(an/Mr))
    dDOC <- Deg - Cu + Enz*k + pmax(0, -(1 - j)*(an/Mr))
    dR <- Cu - f*R
    dS <- pmax(0, chi*an*0.6, 0) - pmax(0, -an/Mr)
    dEnz <- pmax(0, (1-chi)*an*Yue, 0) - Enz*k
    
    return(list(c(dCs, dDOC, dR, dS, dEnz), 
                r=r, Protinc = Protinc, 
                DNAc = DNAc))
  })
}