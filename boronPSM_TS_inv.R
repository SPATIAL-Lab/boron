model{
############################################################################################
#    LIKELIHOOD FUNCTION: evaluates the data against modeled values
############################################################################################    
 
  for (i in 1:ai.d11B){
  d11Bf.data ~ dnorm(d11Bf[i], d11Bf.p)
  d11Bf.p = 1/d11Bfu.data^2
  }
  
  for (i in 1:ai.MgCa){
  mgcaf.data ~ dnorm(mgcaf[i], mgcaf.p)
  mgcaf.p = 1/0.03^2 
  }
  
  for (i in 1:ai.d18O){
  d18Of.data ~ dnorm(d18Of[i], d18Of.p)
  d18Of.p = 1/0.1^2
  }
  
  
############################################################################################
#    DETERMINE 'm' AND 'c' VALUES USING MODERN SPECIES CALIBRATION DATA (input in driver)
############################################################################################    
  
  # Average measurement precision for modern calibration points
  Bmeas.sd <- 0.2 
  
  for(i in 1:length(d11Bcb)){
    d11Bcb[i] ~ dnorm(d11Bcb.m[i], 1/Bmeas.sd^2) 
    d11Bcb.m[i] = d11Bcfo[i]*m + c
  }


############################################################################################
#    PROXY SYSTEM MODEL
############################################################################################    
  
  
  # INPUT VALUES FOR CALCULATIONS 
  
  # Set modern concentrations for Mg, Ca, and SO4
  xcam <- 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm <- 52.8171 # modern [Mg] (mmol kg^-1)
  xso4m <- 28.24  # modern [SO4] (mmol kg^-1)
  mgcaswm <- xmgm/xcam # modern Mg/Ca of seawater

  
  # Set fractionation factor: Klochko et al. (2006)
  alpha ~ dnorm(alpha.m, alpha.p)     
  alpha.m = 1.0272    #gaussian mean
  alpha.p = 1/0.0003^2   #gaussian precision
  epsilon <- (alpha - 1)*1000  # Compute epsilon from alpha

  for (i in 1:ai.prox){
    # Hollis et al., 2019 input: set H, B and A distributions (i.e., Evans et al., 2012, 2016b)
    # Includes calibration uncertainty in these terms
    Hp ~ dnorm(Hp.m, Hp.p) # nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer)
    Hp.m = Hp.mean
    Hp.p = 1/Hp.sd^2
    Bmod ~ dnorm(Bmod.m, Bmod.p) # modern pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
    Bmod.m = Bmod.mean  
    Bmod.p = 1/Bmod.sd^2
    A ~ dnorm(A.m, A.p) # Exponenital constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans )
    A.m = A.mean    
    A.p = 1/A.sd^2
    
    
    # CARB CHEM EQUILIBRIUM CONSTANT CALCULATIONS FOLLOWING ZEEBE AND TYRRELL (2019)
    
    # Calculate equil. constants using salinity and temp:
    temp <- tempC[ai.prox[i]]+273.15
    Ks1m_st <-exp(2.83655-2307.1266/temp-1.5529413*(log(temp))-((0.20760841+4.0484/temp)*sqrt(sal[ai.prox[i]]))+0.0846834*sal[ai.prox[i]]-0.00654208*(sal[ai.prox[i]]^1.5)+log(1-(0.001005*sal[ai.prox[i]])))
    Ks2m_st <- exp(-9.226508-3351.6106/temp-0.2005743*(log(temp))-((0.106901773+23.9722/temp)*sqrt(sal[ai.prox[i]]))+0.1130822*sal[ai.prox[i]]-0.00846934*(sal[ai.prox[i]]^1.5)+log(1-(0.001005*sal[ai.prox[i]])))
    logKsspcm_st <- ((-171.9065-0.077993*temp+2839.319/temp+71.595*(log(temp)/log(10))+(-0.77712+0.0028426*temp+178.34/temp)*(sal[ai.prox[i]]^0.5)-0.07711*sal[ai.prox[i]]+0.0041249*(sal[ai.prox[i]]^1.5)))
    Ksspcm_st <- 10^(logKsspcm_st)
    lnKsB_st <- ((-8966.9-2890.53*sal[ai.prox[i]]^0.5-77.942*sal[ai.prox[i]]+1.728*sal[ai.prox[i]]^1.5-0.0996*sal[ai.prox[i]]^2)/temp)+148.0248+137.1942*sal[ai.prox[i]]^0.5+1.62142*sal[ai.prox[i]]-(24.4344+25.085*sal[ai.prox[i]]^0.5+0.2474*sal[ai.prox[i]])*(log(temp))+(0.053105*sal[ai.prox[i]]^0.5*temp)
    KsB_st <- exp(lnKsB_st)
    Ksw_st <- exp(148.96502-13847.26/temp-23.6521*(log(temp))+(118.67/temp-5.977+1.0495*(log(temp)))*(sal[ai.prox[i]]^0.5)-0.01615*sal[ai.prox[i]])
    K0 <- exp(9345.17/temp-60.2409+23.3585*(log(temp/100))+sal[ai.prox[i]]*(0.023517-0.00023656*temp+0.0047036*((temp/100)^2)))
    
    # Adjust equil. constants for the effect of pressure (Millero 1995):
    delV1 <- (-25.50)+0.1271*tempC[ai.prox[i]]
    delV2 <- (-15.82)+(-0.0219*tempC[ai.prox[i]])
    delVspc <- (-48.76)+(0.5304*tempC[ai.prox[i]])
    delVB <- (-29.48)+0.1622*tempC[ai.prox[i]]+(2.608/1000)*tempC[ai.prox[i]]^2
    delVw <- (-25.60)+0.2324*tempC[ai.prox[i]]+(-3.6246/1000)*tempC[ai.prox[i]]^2
    
    delk1 <- (-3.08/1000)+(0.0877/1000)*tempC[ai.prox[i]]
    delk2 <- (1.13/1000)+(-0.1475/1000)*tempC[ai.prox[i]]
    delkspc <- (-11.76/1000)+(0.3692/1000)*tempC[ai.prox[i]]
    delkB <- -2.84/1000
    delkw <- (-5.13/1000)+(0.0794/1000)*tempC[ai.prox[i]]
    
    R <- 83.131 # constant (cm^3 bar mol^-1 K^-1)
    
    Ks1m <- (exp(-((delV1/(R*temp))*press)+((0.5*delk1)/(R*temp))*press^2))*Ks1m_st
    Ks2m <- (exp(-((delV2/(R*temp))*press)+((0.5*delk2)/(R*temp))*press^2))*Ks2m_st
    Ksspcm <- (exp(-((delVspc/(R*temp))*press)+((0.5*delkspc)/(R*temp))*press^2))*Ksspcm_st
    KsB <- (exp(-((delVB/(R*temp))*press)+((0.5*delkB)/(R*temp))*press^2))*KsB_st
    Ksw <- (exp(-((delVw/(R*temp))*press)+((0.5*delkw)/(R*temp))*press^2))*Ksw_st
    
    # K*1, K*2, and K*spc are corrected for past seawater [Ca], [Mg], and [SO4] following ZT19:
    
    # Define ZT19 table 2 sensitivity parameters (si_j)
    s1_ca <- 5/1000
    s1_mg <- 17/1000
    s1_so4 <- 208/1000
    s2_ca <- 157/1000
    s2_mg <- 420/1000
    s2_so4 <- 176/1000
    sspc_ca <- 185/1000
    sspc_mg <- 518/1000
    sspc_so4 <- 106/1000
    
    # Compute K*is for past seawater composition 
    Ks1 <- Ks1m*(1+(s1_ca*(xca[ai.prox[i]]/xcam-1)+s1_mg*(xmg[ai.prox[i]]/xmgm-1)+s1_so4*(xso4[ai.prox[i]]/xso4m-1)))
    Ks2 <- Ks2m*(1+(s2_ca*(xca[ai.prox[i]]/xcam-1)+s2_mg*(xmg[ai.prox[i]]/xmgm-1)+s2_so4*(xso4[ai.prox[i]]/xso4m-1)))
    Ksspc <- Ksspcm*(1+(sspc_ca*(xca[ai.prox[i]]/xcam-1)+sspc_mg*(xmg[ai.prox[i]]/xmgm-1)+sspc_so4*(xso4[ai.prox[i]]/xso4m-1)))
  
    
    # DETERMINE FORAMINIFERAL D18O, D11B, AND MG/CA 
    
    # GJB - In my experiments this is likely to be the biggest source
    # of instability...I think for certain combos of pco2 and dic the 
    # roots h1 and h2 are undefined?
    
    # Compute pH from [co2] and DIC
    fco2 <- pco2[ai.prox[i]]*0.9968
    co2 <- fco2*K0
    qa <- co2-dic[ai.prox[i]]
    qb <- co2*Ks1
    qc <- co2*Ks1*Ks2
    h1 <- (-qb+sqrt((qb^2)-(4*qa*qc)))/(2*qa)
    h2 <- (-qb-sqrt((qb^2)-(4*qa*qc)))/(2*qa)
    pH <- -(log(max(h1, h2))/log(10))                      
    
    # Compute d11Bforam from pH and d11Bsw
    pKsB <- -(log(KsB)/log(10))
    t1 <- 10^(pKsB-pH)
    d11Bb <- ((t1*epsilon)-(t1*d11Bsw[ai.prox[i]])-d11Bsw[ai.prox[i]])/(-((t1*alpha)+1))
    c.final <- c + c.correction # adjusts final c value if desired - correction specified in driver
    d11Bf <- m*d11Bb + (c.final)
     
    # Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
    sw.sens ~ dnorm(0.558, 1/0.03^2) # Uncertainty represents std dev of regression slope in GEOSECS obs. reported in Charles and Fairbanks (1990)
    d18Osw.sc <- d18Osw[ai.prox[i]] + (sw.sens*(sal-35)) 
    d18Oswpdb <- d18Osw.sc -0.27
    d18Of.pr <- d18Oswpdb + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC[ai.prox[i]])))^0.5))/(2*0.09))
    indexop ~ dnorm(indexop.m, indexop.p)
    indexop.m = seccal/100
    indexop.p = 1/(seccal.u/100)^2
    d18Of <- d18Of.pr + indexop*Dd18Oseccal
    
    # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
    mgcasw <- (xmg[ai.prox[i]]/xca[ai.prox[i]])     
    Bcorr <- ((mgcasw^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr <- Bcorr*(exp(A*tempC[ai.prox[i]]))
    pHcorrco ~ dnorm(pHcorrco.m, pHcorrco.p)
    pHcorrco.m = 0.70
    pHcorrco.p = 1/0.09^2
    mgca_sal <- mgca_corr / (1-(8.05-pH)*pHcorrco)
    salcorrco ~ dnorm(salcorrco.m, salcorrco.p)
    salcorrco.m = 0.042
    salcorrco.p = 1/0.004^2
    mgcaf <- mgca_sal / (1-(sal[ai.prox[i]]-35)*salcorrco)
  }
  
############################################################################################
#    ENVIRONMENTAL MODEL
############################################################################################    
 
  # Environmental time-dependent priors
    
  sal <- list()
  tempC <- list()
  xca <- list()
  xmg <- list()
  xso4 <- list()
  d11Bsw <- list()
  d18Osw <- list()
  pco2 <- list()
  dic <- list()
    
  for (i in 2:n.steps){
    
  # Salinity (ppt)  
  sal.m = 35  
  sal.p = 1/0.5^2    
  sal[1] ~ dnorm(sal.m, sal.p)     
  sal.pc = 1
  sal.sig[i] ~ dnorm(0, sal.pc)
  sal[i] = sal[i-1] + sal.sig[i]
  
  # Temp in C
  tempC.m = 30   
  tempC.p = 1/10^2 
  tempC[1] ~ dnorm(tempC.m, tempC.p)     
  tempC.pc = 1
  tempC.sig[i] ~ dnorm(0, tempC.pc)
  tempC[i] = tempC[i-1] + tempC.sig[i]
  
  # [Ca] (mmol kg^-1)
  xca.m = 17
  xca.p = 1/0.5^2
  xca[1] ~ dnorm(xca.m, xca.p)       
  xca.pc = 0.1
  xca.sig[i] ~ dnorm(0, xca.pc)
  xca[i] = xca[i-1] + xca.sig[i]
  
  # [Mg] (mmol kg^-1)
  xmg.m = 36
  xmg.p = 1/0.5^2
  xmg[1] ~ dnorm(xmg.m, xmg.p)      
  xmg.pc = 0.1
  xmg.sig[i] ~ dnorm(0, xmg.pc)
  xmg[i] = xmg[i-1] + xmg.sig[i]
  
  # [SO4] (mmol kg^-1)
  xso4.m = 14
  xso4.p = 1/0.5^2
  xso4[1] ~ dnorm(xso4.m, xso4.p)      
  xso4.pc = 0.1
  xso4.sig[i] ~ dnorm(0, xso4.pc)
  xso4[i] = xso4[i-1] + xso4.sig[i]
  
  # d11B of seawater (per mille SRM-951) 
  d11Bsw.m = 38.45
  d11Bsw.p = 1/0.5^2
  d11Bsw[1] ~ dnorm(d11Bsw.m, d11Bsw.p)    
  d11Bsw.pc = 0.1
  d11Bsw.sig[i] ~ dnorm(0, d11Bsw.pc)
  d11Bsw[i] = d11Bsw[i-1] + d11Bsw.sig[i]
  
  # d18O of seawater (per mille SMOW) 
  d18Osw.m = -1
  d18Osw.p = 1/0.5^2
  d18Osw[1] ~ dnorm(d18Osw.m, d18Osw.p)    
  d18Osw.pc = 0.1
  d18Osw.sig[i] ~ dnorm(0, d18Osw.pc)
  d18Osw[i] = d18Osw[i-1] + d18Osw.sig[i] 
  
  # Atmospheric pCO2 (atm)
  pco2.l = 300e-6   
  pco2.u = 2500e-6   
  pco2[1] ~ dunif(pco2.l, pco2.u)    
  pco2.pc = 10e-6
  pco2.sig[i] ~ dnorm(0, pco2.pc)
  pco2[i] = pco2[i-1] + pco2.sig[i]
  
  # DIC (mol kg^-1) - make change in DIC temp dependent (C cycle model)?
  dic.m = 0.0022
  dic.p = 1/0.0002^2
  dic[1] ~ dnorm(dic.m, dic.p)
  dic.pc = 0.00002
  dic.sig[i] ~ dnorm(0, dic.pc)
  dic[i] = dic[i-1] + dic.sig[i]
  # dic[i] = dic[i-1] + (tempC[i]-tempC[i-1])*0.00005 + dic.sig[i]
  
  sal[[i]] <- sal[i]
  tempC[[i]] <- tempC[i]
  xca[[i]] <- xca[i]
  xmg[[i]] <- xmg[i]
  xso4[[i]] <- xso4[i]
  d11Bsw[[i]] <- d11Bsw[i]
  d18Osw[[i]] <- d18Osw[i]
  pco2[[i]] <- pco2[i]
  dic[[i]] <- dic[i]
  }

  # Age index vector for prior time bins
  ai.env = ceiling((ages.max - ages) / ages.bin)  
  
  
  # Time independent priors 
  
  # 'm' boron vital effect calibration parm 
  m ~ dnorm(m.mean, 1/m.sd^2)  
  
  # 'c' boron vital effect calibration parm
  c ~ dnorm(c.mean, 1/c.sd^2)     
  
  # Pressure (bar)
  press ~ dnorm(press.m, press.p)    
  press.p = 1/press.sd^2 
  
  
}


  