
model{
  ############################################################################################
  #    LIKELIHOOD FUNCTION: evaluates the data against modeled values
  ############################################################################################    

  # Gaussian precision for d18O measurements 
  d18Of.p = 1/0.08^2
  
  for (i in 1:length(ai.mgca)){
    mgcaf.data[i] ~ dnorm(mgcaf[ai.mgca[i]], mgcaf.p[i])
    mgcaf.p[i] = 1/(mgcafu.data[i])^2  # Gaussian precision for Mg/Ca measurements (1sd = 1.5%) based on long-term consistency standard reproducibility 
  }
  
  for (i in 1:length(ai.d18O)){
    d18Of.data[i] ~ dnorm(d18Of[ai.d18O[i]], d18Of.p)
  }
  
  
  ############################################################################################
  #    PROXY SYSTEM MODEL
  ############################################################################################    
  
  
  # INPUT VALUES FOR CALCULATIONS 
  
  # Set modern concentrations for Mg, Ca, and SO4
  xcam = 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm = 52.8171 # modern [Mg] (mmol kg^-1)
  mgcaswm <- xmgm/xcam # modern Mg/Ca of seawater
  
  # Index of diagenetic overprint
  indexop ~ dnorm((seccal/100), (1/(seccal.u/100)^2))
  
  # Sensitivity of seawater d18O to salinity 
  sw.sens ~ dnorm(0.558, 1/0.03^2) # regression slope and uncertainty (i.e., from standard deviation) of GEOSECS obs. reported in Charles and Fairbanks (1990) after Duplessey et al. (1991)
  
  # Correction for Mg/Ca salinity effect 
  salcorrco ~ dnorm(0.042, (1/0.004^2))
  
  # Set H, B and A distributions (i.e., Evans et al., 2012, 2016b)
  # Includes calibration uncertainty in these terms
  Hp ~ dnorm(Hp.m, Hp.p) # nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer)
  Hp.m = Hp.mean
  Hp.p = 1/Hp.sd^2
  Bmod ~ dnorm(Bmod.m, Bmod.p) # modern pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
  Bmod.m = Bmod.mean  
  Bmod.p = 1/Bmod.sd^2
  A ~ dnorm(A.m, A.p) # Exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans )
  A.m = A.mean    
  A.p = 1/A.sd^2
  
  R = 83.131 # constant (cm^3 bar mol^-1 K^-1)
  
  
  for (i in 1:length(ai.prox)){
 
    # DETERMINE FORAMINIFERAL D18O AND MG/CA 
    
    # Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
    d18Osw.sc[i] <- d18Osw[ai.prox[i]] + (sw.sens*(sal[ai.prox[i]]-35))
    d18Oswpdb[i] <- d18Osw.sc[i] -0.27
    d18Of.pr[i] <- d18Oswpdb[i] + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC[ai.prox[i]])))^0.5))/(2*0.09))
    d18Of[i] <- d18Of.pr[i]*(1-indexop) + indexop*d18Oseccal
    
    # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach (CURRENTLY NO CARB CHEM CORRECTION IN THIS VERSION - COMMENTED OUT)
    mgcasw[i] <- (xmg[ai.prox[i]]/xca[ai.prox[i]])     
    Bcorr[i] <- ((mgcasw[i]^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr[i] <- Bcorr[i]*(exp(A*tempC[ai.prox[i]]))
    mgcaf[i] <- mgca_corr[i] / (1-(sal[ai.prox[i]]-35)*salcorrco)
    
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
  sal.phi ~ dbeta(5,2)         
  sal.eps[1] = 0                 
  sal.tau ~ dgamma(1e2, 5e-3) 
  
  # Temp in C
  tempC[1] ~ dnorm(tempC.m, tempC.p) 
  tempC.phi ~ dbeta(5,2) 
  tempC.eps[1] = 0 
  tempC.tau ~ dgamma(10, 10)
  
  # [Ca] (mmol kg^-1)
  xca[1] ~ dnorm(xca.m, xca.p)
  xca.phi ~ dbeta(5,2) 
  xca.eps[1] = 0 
  xca.tau ~ dgamma(1e3, 1e-3)
  
  # [Mg] (mmol kg^-1)
  xmg[1] ~ dnorm(xmg.m, xmg.p)      
  xmg.phi ~ dbeta(5,2)  
  xmg.eps[1] = 0 
  xmg.tau ~ dgamma(1e3, 1e-5)
  
  # d18O of seawater (per mille SMOW) 
  d18Osw[1] ~ dnorm(d18Osw.m, d18Osw.p)
  d18Osw.phi ~ dbeta(5,2) 
  d18Osw.eps[1] = 0 
  d18Osw.tau ~ dgamma(1e3, 1e-3)
  
  
  # Environmental time-dependent variables 
  
  for (i in 2:n.steps){
    
    # Salinity (ppt)  
    sal.pc[i] = sal.tau*((1-sal.phi^2)/(1-sal.phi^(2*dt[i-1])))   
    sal.eps[i] ~ dnorm(sal.eps[i-1]*(sal.phi^dt[i-1]), sal.pc[i])T(-0.3, 0.3)
    sal[i] = sal[1] * (1 + sal.eps[i])

    # Temp in C
    tempC.pc[i] = tempC.tau*((1-tempC.phi^2)/(1-tempC.phi^(2*dt[i-1]))) 
    tempC.eps[i] ~ dnorm(tempC.eps[i-1]*(tempC.phi^dt[i-1]), tempC.pc[i])T(-12, 12)
    tempC[i] = tempC[1] + tempC.eps[i]
    
    # [Ca] (mmol kg^-1); linear decline in Cenozoic follows Holland et al. (2020) 
    xca.pc[i] = xca.tau*((1-xca.phi^2)/(1-xca.phi^(2*dt[i-1])))
    xca.eps[i] ~ dnorm(xca.eps[i-1]*(xca.phi^dt[i-1]), xca.pc[i])T(-0.3, 0.3)
    xca[i] = xca[1] * (1 + xca.eps[i]) 
    
    # [Mg] (mmol kg^-1); imposing max realistic decline in Mg/Casw; decrease suggested by paired benthic Mg/Ca+d18O from ODP Site 1209 (data suggest even greater decline than what's imposed here; Mg decrease used here is derived from max observed shifts in Mg/Casw over Cenozoic as in Holland et al. (2020)
    xmg.pc[i] = xmg.tau*((1-xmg.phi^2)/(1-xmg.phi^(2*dt[i-1])))
    xmg.eps[i] ~ dnorm(xmg.eps[i-1]*(xmg.phi^dt[i-1]), xmg.pc[i])T(-0.3, 0.3)
    xmg[i] = xmg[1] * (1 + xmg.eps[i]) 
    
    # d18O of seawater (per mille SMOW) 
    d18Osw.pc[i] = d18Osw.tau*((1-d18Osw.phi^2)/(1-d18Osw.phi^(2*dt[i-1])))
    d18Osw.eps[i] ~ dnorm(d18Osw.eps[i-1]*(d18Osw.phi^dt[i-1]), d18Osw.pc[i])T(-0.2, 0.2)
    d18Osw[i] = d18Osw[1] * (1 + d18Osw.eps[i])
  }
}


