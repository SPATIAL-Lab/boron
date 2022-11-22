model{
############################################################################################
#    LIKELIHOOD FUNCTION: evaluates the data against modeled values
############################################################################################    
  
  # Gaussian precision for Mg/Ca measurements 
  mgcaf.p = 1/0.3^2 

 
  for(i in 1:n.steps){
    mgca[i] ~ dnorm(xmg[i] / xca[i], mgcaf.p)
  }
  
############################################################################################
#    ENVIRONMENTAL MODEL
############################################################################################    
 
  # Environmental time-dependent prior initial conditions 
        # .phi = temporal autocorrelation of a parameter
        # .eps = random walk error term (i.e., eqn 5 in Bowen et al., 2020 CotP)
        # .tau = error precision for dt = 1
        # .pc = error precision of random walk error term (i.e., .eps)

  
  # Salinity (ppt)  
  sal[1] ~ dnorm(sal.m, sal.p)  
  sal.phi ~ dbeta(2,2)         
  sal.eps[1] = 0                 
  sal.tau ~ dgamma(1e2, 5e-3) 
  
  # Temp in C
  tempC[1] ~ dnorm(tempC.m, tempC.p) 
  tempC.phi ~ dbeta(2,5) 
  tempC.eps[1] = 0 
  tempC.tau ~ dgamma(10, 10)
    
  # [Ca] (mmol kg^-1)
  xca[1] ~ dnorm(xca.m, xca.p)
  xca.phi ~ dbeta(2,2) 
  xca.eps[1] = 0 
  xca.tau ~ dgamma(1e2, 1e-2)

  # [Mg] (mmol kg^-1)
  xmg[1] ~ dnorm(xmg.m, xmg.p)      
  xmg.phi ~ dbeta(2,2)  
  xmg.eps[1] = 0 
  xmg.tau ~ dgamma(1e2, 7e-3)
    
  # [SO4] (mmol kg^-1)
  xso4[1] ~ dnorm(xso4.m, xso4.p)      
  xso4.phi ~ dbeta(2,2) 
  xso4.eps[1] = 0 
  xso4.tau ~ dgamma(1e2, 4e-4)
  
  # d11B of seawater (per mille SRM-951) 
  d11Bsw[1] ~ dnorm(d11Bsw.m, d11Bsw.p)    
  d11Bsw.phi ~ dbeta(2,5)  
  d11Bsw.eps[1] = 0 
  d11Bsw.tau ~ dgamma(4e3, 1e-3)
    
  # d18O of seawater (per mille SMOW) 
  d18Osw[1] ~ dnorm(d18Osw.m, d18Osw.p)    
  d18Osw.phi ~ dbeta(2,5) 
  d18Osw.eps[1] = 0 
  d18Osw.tau ~ dgamma(1e3, 1e-3)
  
  # pH
  pH[1] ~ dunif(pH.l, pH.u)   
  pH.phi ~ dbeta(2,5)
  pH.eps[1] = 0 
  pH.tau ~ dgamma(100, 1e-1) 
  
  # DIC (mol kg^-1) 
  dic[1] ~ dnorm(dic.m, dic.p)
  dic.phi ~ dbeta(2,2)
  dic.eps[1] = 0
  dic.tau ~ dgamma(1e2, 1e-2)
  
  
  # Environmental time-dependent variables 
  
  for (i in 2:n.steps){
    
  # Salinity (ppt)  
  sal.pc[i] = sal.tau*((1-sal.phi^2)/(1-sal.phi^(2*dt[i-1])))   
  sal.eps[i] ~ dnorm(sal.eps[i-1]*(sal.phi^dt[i-1]), sal.pc[i])
  sal[i] = sal[1] * (1 + sal.eps[i])
  
  # Temp in C
  tempC.pc[i] = tempC.tau*((1-tempC.phi^2)/(1-tempC.phi^(2*dt[i-1]))) 
  tempC.eps[i] ~ dnorm(tempC.eps[i-1]*(tempC.phi^dt[i-1]), tempC.pc[i])
  tempC[i] = tempC[i-1] + tempC.eps[i]
  
  # [Ca] (mmol kg^-1); linear decline in Cenozoic follows Holland et al. (2020) 
#  xca.tdep[i] = -0.0001938*dt[i-1]     
  xca.pc[i] = xca.tau*((1-xca.phi^2)/(1-xca.phi^(2*dt[i-1])))
  xca.eps[i] ~ dnorm(xca.eps[i-1]*(xca.phi^dt[i-1]), xca.pc[i])
  xca[i] = xca[1] * (1 + xca.eps[i])# + xca.tdep[i]
  
  # [Mg] (mmol kg^-1); imposing max realistic decline in Mg/Casw; decrease suggested by paired benthic Mg/Ca+d18O from ODP Site 1209 (data suggest even greater decline than what's imposed here; Mg decrease used here is derived from max observed shifts in Mg/Casw over Cenozoic as in Holland et al. (2020)
#  xmg.tdep[i] = -0.00274*dt[i-1]      
  xmg.pc[i] = xmg.tau*((1-xmg.phi^2)/(1-xmg.phi^(2*dt[i-1])))
  xmg.eps[i] ~ dnorm(xmg.eps[i-1]*(xmg.phi^dt[i-1]), xmg.pc[i])
  xmg[i] = xmg[1] * (1 + xmg.eps[i]) #+ xmg.tdep[i]
  
  # [SO4] (mmol kg^-1)
  xso4.pc[i] = xso4.tau*((1-xso4.phi^2)/(1-xso4.phi^(2*dt[i-1])))
  xso4.eps[i] ~ dnorm(xso4.eps[i-1]*(xso4.phi^dt[i-1]), xso4.pc[i])
  xso4[i] = xso4[1] * (1 + xso4.eps[i])
  
  # d11B of seawater (per mille SRM-951) 
  d11Bsw.pc[i] = d11Bsw.tau*((1-d11Bsw.phi^2)/(1-d11Bsw.phi^(2*dt[i-1])))
  d11Bsw.eps[i] ~ dnorm(d11Bsw.eps[i-1]*(d11Bsw.phi^dt[i-1]), d11Bsw.pc[i])
  d11Bsw[i] = d11Bsw[i-1] + d11Bsw.eps[i]
  
  # d18O of seawater (per mille SMOW) 
  d18Osw.pc[i] = d18Osw.tau*((1-d18Osw.phi^2)/(1-d18Osw.phi^(2*dt[i-1])))
  d18Osw.eps[i] ~ dnorm(d18Osw.eps[i-1]*(d18Osw.phi^dt[i-1]), d18Osw.pc[i])
  d18Osw[i] = d18Osw[i-1] + d18Osw.eps[i]
  
  # pH
  pH.pc[i] = pH.tau*((1-pH.phi^2)/(1-pH.phi^(2*dt[i-1])))
  pH.eps[i] ~ dnorm(pH.eps[i-1]*(pH.phi^dt[i-1]), pH.pc[i])
  pH[i] = pH[i-1] + pH.eps[i]
  
  # DIC (mol kg^-1) - make change in DIC temp dependent (C cycle model)?
  dic.pc[i] = dic.tau*((1-dic.phi^2)/(1-dic.phi^(2*dt[i-1])))
  #temp.dep[i] = (tempC[i]-tempC[i-1])*0.00005
  #dic.eps[i] ~ dnorm((dic.eps[i-1]+temp.dep[i])*(dic.phi^dt[i-1]), dic.pc[i]) 
  dic.eps[i] ~ dnorm(dic.eps[i-1]*(dic.phi^dt[i-1]), dic.pc[i])
  dic[i] = dic[1] * (1 + dic.eps[i]) #+ temp.dep[i]
  
  }

}


  