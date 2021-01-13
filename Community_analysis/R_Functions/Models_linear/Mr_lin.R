Mr_lin<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    #Define time variable parameter
    Mr = I + B * (4.221*time/(30.700+time+time^2/106.021)+1)
    
    #Define equations
    ##Enzymatic degradation - multiplicative equation
    Deg<-Enz * Cs * Ecat
    ##Carbon uptake
    Cu<-0.283475617*S*DOC/(2.512518101 + DOC)
    ##maintnance
    m=S*Mr
    ##Reserves available
    an=3.05208273377997*R-m
    ##Respiration rate
    r=pmax(chi*an*(1-0.6), 0) + pmax((1-chi)*an*(1-Yue), 0) + m
    
    ##Proteins and DNA
    Protinc=0.293821346*R+0.063372172*S
    DNAc=0.004533891*S
    
    #States
    dCs <- -Deg
    dDOC <- Deg - Cu
    dR <- Cu - 3.05208273377997*R
    dS <- pmax(0, chi*an*0.6, 0) + pmin(0, an)
    dEnz <- pmax(0, (1-chi)*an*Yue, 0)
    
    return(list(c(dCs, dDOC, dR, dS, dEnz), 
                r=r, Protinc = Protinc, 
                DNAc = DNAc))
  })
}