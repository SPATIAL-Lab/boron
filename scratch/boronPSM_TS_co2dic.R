model{
############################################################################################
#    LIKELIHOOD FUNCTION: evaluates the data against modeled values
############################################################################################    
  
  # Gaussian precision for Mg/Ca measurements 
  mgcaf.p = 1/0.03^2 
  # Gaussian precision for d18O measurements 
  d18Of.p = 1/0.1^2
  
  for (i in 1:length(ai.d11B1)){
  d11Bf.data1[i] ~ dnorm(d11Bf.1[ai.d11B1[i]], d11Bf.p1[i])
  d11Bf.p1[i] = 1/d11Bfu.data1[i]^2
  }
  
  for (i in 1:length(ai.d11B2)){
  d11Bf.data2[i] ~ dnorm(d11Bf.2[ai.d11B2[i]], d11Bf.p2[i])
  d11Bf.p2[i] = 1/d11Bfu.data2[i]^2
  }
  
  for (i in 1:length(ai.mgca)){
  mgcaf.data[i] ~ dnorm(mgcaf[ai.mgca[i]], mgcaf.p)
  }
  
  for (i in 1:length(ai.d18O)){
  d18Of.data[i] ~ dnorm(d18Of[ai.d18O[i]], d18Of.p)
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

  
  # Define ZT19 table 2 sensitivity parameters (si_j) for changing sw chemistry on carb chem
  s1_ca <- 5/1000
  s1_mg <- 17/1000
  s1_so4 <- 208/1000
  s2_ca <- 157/1000
  s2_mg <- 420/1000
  s2_so4 <- 176/1000
  sspc_ca <- 185/1000
  sspc_mg <- 518/1000
  sspc_so4 <- 106/1000
 
 
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
  
  R <- 83.131 # constant (cm^3 bar mol^-1 K^-1)
  
    
  for (i in 1:length(ai.prox)){
    # CARB CHEM EQUILIBRIUM CONSTANT CALCULATIONS FOLLOWING ZEEBE AND TYRRELL (2019)
    
    # Calculate equil. constants using salinity and temp:
    temp[i] <- tempC[ai.prox[i]]+273.15
    Ks1m_st[i] <-exp(2.83655-2307.1266/temp[i]-1.5529413*(log(temp[i]))-((0.20760841+4.0484/temp[i])*sqrt(sal[ai.prox[i]]))+0.0846834*sal[ai.prox[i]]-0.00654208*(sal[ai.prox[i]]^1.5)+log(1-(0.001005*sal[ai.prox[i]])))
    Ks2m_st[i] <- exp(-9.226508-3351.6106/temp[i]-0.2005743*(log(temp[i]))-((0.106901773+23.9722/temp[i])*sqrt(sal[ai.prox[i]]))+0.1130822*sal[ai.prox[i]]-0.00846934*(sal[ai.prox[i]]^1.5)+log(1-(0.001005*sal[ai.prox[i]])))
    logKsspcm_st[i] <- ((-171.9065-0.077993*temp[i]+2839.319/temp[i]+71.595*(log(temp[i])/log(10))+(-0.77712+0.0028426*temp[i]+178.34/temp[i])*(sal[ai.prox[i]]^0.5)-0.07711*sal[ai.prox[i]]+0.0041249*(sal[ai.prox[i]]^1.5)))
    Ksspcm_st[i] <- 10^(logKsspcm_st[i])
    lnKsB_st[i] <- ((-8966.9-2890.53*sal[ai.prox[i]]^0.5-77.942*sal[ai.prox[i]]+1.728*sal[ai.prox[i]]^1.5-0.0996*sal[ai.prox[i]]^2)/temp[i])+148.0248+137.1942*sal[ai.prox[i]]^0.5+1.62142*sal[ai.prox[i]]-(24.4344+25.085*sal[ai.prox[i]]^0.5+0.2474*sal[ai.prox[i]])*(log(temp[i]))+(0.053105*sal[ai.prox[i]]^0.5*temp[i])
    KsB_st[i] <- exp(lnKsB_st[i])
    Ksw_st[i] <- exp(148.96502-13847.26/temp[i]-23.6521*(log(temp[i]))+(118.67/temp[i]-5.977+1.0495*(log(temp[i])))*(sal[ai.prox[i]]^0.5)-0.01615*sal[ai.prox[i]])
    K0[i] <- exp(9345.17/temp[i]-60.2409+23.3585*(log(temp[i]/100))+sal[ai.prox[i]]*(0.023517-0.00023656*temp[i]+0.0047036*((temp[i]/100)^2)))
    
    # Adjust equil. constants for the effect of pressure (Millero 1995):
    delV1[i] <- (-25.50)+0.1271*tempC[ai.prox[i]]
    delV2[i] <- (-15.82)+(-0.0219*tempC[ai.prox[i]])
    delVspc[i]<- (-48.76)+(0.5304*tempC[ai.prox[i]])
    delVB[i] <- (-29.48)+0.1622*tempC[ai.prox[i]]+(2.608/1000)*tempC[ai.prox[i]]^2
    delVw[i] <- (-25.60)+0.2324*tempC[ai.prox[i]]+(-3.6246/1000)*tempC[ai.prox[i]]^2
    
    delk1[i] <- (-3.08/1000)+(0.0877/1000)*tempC[ai.prox[i]]
    delk2[i] <- (1.13/1000)+(-0.1475/1000)*tempC[ai.prox[i]]
    delkspc[i] <- (-11.76/1000)+(0.3692/1000)*tempC[ai.prox[i]]
    delkB[i] <- -2.84/1000
    delkw[i] <- (-5.13/1000)+(0.0794/1000)*tempC[ai.prox[i]]
    
    Ks1m[i] <- (exp(-((delV1[i]/(R*temp[i]))*press)+((0.5*delk1[i])/(R*temp[i]))*press^2))*Ks1m_st[i]
    Ks2m[i] <- (exp(-((delV2[i]/(R*temp[i]))*press)+((0.5*delk2[i])/(R*temp[i]))*press^2))*Ks2m_st[i]
    Ksspcm[i] <- (exp(-((delVspc[i]/(R*temp[i]))*press)+((0.5*delkspc[i])/(R*temp[i]))*press^2))*Ksspcm_st[i]
    KsB[i] <- (exp(-((delVB[i]/(R*temp[i]))*press)+((0.5*delkB[i])/(R*temp[i]))*press^2))*KsB_st[i]
    Ksw[i] <- (exp(-((delVw[i]/(R*temp[i]))*press)+((0.5*delkw[i])/(R*temp[i]))*press^2))*Ksw_st[i]
    
    # K*1, K*2, and K*spc are corrected for past seawater [Ca], [Mg], and [SO4] following ZT19
    Ks1[i] <- Ks1m[i]*(1+(s1_ca*(xca[ai.prox[i]]/xcam-1)+s1_mg*(xmg[ai.prox[i]]/xmgm-1)+s1_so4*(xso4[ai.prox[i]]/xso4m-1)))
    Ks2[i] <- Ks2m[i]*(1+(s2_ca*(xca[ai.prox[i]]/xcam-1)+s2_mg*(xmg[ai.prox[i]]/xmgm-1)+s2_so4*(xso4[ai.prox[i]]/xso4m-1)))
    Ksspc[i] <- Ksspcm[i]*(1+(sspc_ca*(xca[ai.prox[i]]/xcam-1)+sspc_mg*(xmg[ai.prox[i]]/xmgm-1)+sspc_so4*(xso4[ai.prox[i]]/xso4m-1)))
  
    
    # DETERMINE FORAMINIFERAL D18O, D11B, AND MG/CA 
    
    # GJB - In my experiments this is likely to be the biggest source
    # of instability...I think for certain combos of pco2 and dic the 
    # roots h1 and h2 are undefined?
    
    # Compute pH from [co2] and DIC
    fco2[i] <- pco2[ai.prox[i]]*0.9968
    co2[i] <- fco2[i]*K0[i]
    qa[i] <- co2[i]-dic[ai.prox[i]]
    qb[i] <- co2[i]*Ks1[i]
    qc[i] <- co2[i]*Ks1[i]*Ks2[i]
    h1[i] <- (-qb[i]+sqrt((qb[i]^2)-(4*qa[i]*qc[i])))/(2*qa[i])
    h2[i] <- (-qb[i]-sqrt((qb[i]^2)-(4*qa[i]*qc[i])))/(2*qa[i])
    pH[i] <- -(log(max(h1[i], h2[i]))/log(10))                      
    
    # Compute d11Bforam from pH and d11Bsw
    pKsB[i] <- -(log(KsB[i])/log(10))
    t1[i] <- 10^(pKsB[i]-pH[i])
    d11Bb[i] <- ((t1[i]*epsilon)-(t1[i]*d11Bsw[ai.prox[i]])-d11Bsw[ai.prox[i]])/(-((t1[i]*alpha)+1))
    
    c.final1[i] <- c.1 + c.correction1 # adjusts final c value if desired - correction specified in driver
    d11Bf.1[i] <- m.1*d11Bb[i]+ (c.final1[i])
   
    c.final2[i] <- c.2 + c.correction2 # adjusts final c value if desired - correction specified in driver
    d11Bf.2[i] <- m.2*d11Bb[i]+ (c.final2[i])
    
    # Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
    sw.sens[i] ~ dnorm(0.558, 1/0.03^2) # Uncertainty represents std dev of regression slope in GEOSECS obs. reported in Charles and Fairbanks (1990)
    d18Osw.sc[i] <- d18Osw[ai.prox[i]] + (sw.sens[i]*(sal[ai.prox[i]]-35)) 
    d18Oswpdb[i] <- d18Osw.sc[i] -0.27
    d18Of.pr[i] <- d18Oswpdb[i] + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC[ai.prox[i]])))^0.5))/(2*0.09))
    indexop[i] ~ dnorm((seccal/100), (1/(seccal.u/100)^2))
    d18Of[i] <- d18Of.pr[i] + indexop[i]*Dd18Oseccal
    
    # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
    mgcasw[i] <- (xmg[ai.prox[i]]/xca[ai.prox[i]])     
    Bcorr[i] <- ((mgcasw[i]^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr[i] <- Bcorr[i]*(exp(A*tempC[ai.prox[i]]))
    pHcorrco[i] ~ dnorm(0.70, (1/0.09^2))
    mgca_sal[i] <- mgca_corr[i] / (1-(8.05-pH[i])*pHcorrco[i])
    salcorrco[i] ~ dnorm(0.042, (1/0.004^2))
    mgcaf[i] <- mgca_sal[i] / (1-(sal[ai.prox[i]]-35)*salcorrco[i])
    
  }
  
############################################################################################
#    ENVIRONMENTAL MODEL
############################################################################################    
 
  # Environmental time-dependent prior initial conditions ("-.m" = mean, "-.p" = precision, 
  # "-.l" = lower, "-.u" = upper) and change terms ("-.mc" = mean change, "-.pc" = change precision)
  
  # Salinity (ppt)  
  sal.m = 35  
  sal.p = 1/0.5^2    
  sal[1] ~ dnorm(sal.m, sal.p)     
  sal.pc = 1/0.001^2
  sal.mc = 1
  
  # Temp in C
  tempC.m = 30   
  tempC.p = 1/5^2 
  tempC[1] ~ dnorm(tempC.m, tempC.p)     
  tempC.pc = 1
  tempC.mc = 0
  
  # [Ca] (mmol kg^-1)
  xca.m = 21.5842 # calculated using Holland et al., 2020 @ 59 Ma
  xca.p = 1/0.5^2
  xca[1] ~ dnorm(xca.m, xca.p)       
  xca.pc = 1e6
  xca.mc = -0.1938*ages.bin
  
  # [Mg] (mmol kg^-1)
  xmg.m = 63.783  # 68.090 large decrease
  xmg.p = 1/0.5^2
  xmg[1] ~ dnorm(xmg.m, xmg.p)      
  xmg.pc = 1e6
  xmg.mc = -2.958*ages.bin  # -4.819*ages.bin large decrease
  
  # [SO4] (mmol kg^-1)
  xso4.m = 14
  xso4.p = 1/0.5^2
  xso4[1] ~ dnorm(xso4.m, xso4.p)      
  xso4.pc = 1e6
  xso4.mc = 0
  
  # d11B of seawater (per mille SRM-951) 
  d11Bsw.m = 38.45
  d11Bsw.p = 1/0.5^2
  d11Bsw[1] ~ dnorm(d11Bsw.m, d11Bsw.p)    
  d11Bsw.pc = 1/0.05^2
  d11Bsw.mc = 0 
  
  # d18O of seawater (per mille SMOW) 
  d18Osw.m = -1
  d18Osw.p = 1/0.5^2
  d18Osw[1] ~ dnorm(d18Osw.m, d18Osw.p)    
  d18Osw.pc = 1/0.01^2
  d18Osw.mc = 0
  
  # Atmospheric pCO2 (atm)
  pco2.m = 800e-6   
  pco2.p = 500e-6   
  pco2[1] ~ dnorm(pco2.m, pco2.p)I(200e-6,1000e-6) 
  pco2.pc = 1e6
  pco2.mc = 1
  
  # DIC (mol kg^-1) - make change in DIC temp dependent (C cycle model)?
  dic.m = 0.0022
  dic.p = 1/0.0001^2
  dic[1] ~ dnorm(dic.m, dic.p)I(0.002,0.0022)
  dic.pc = 1e8
  dic.mc = 1
  
  
  # Environmental time-dependent priors
  
  for (i in 2:n.steps){
    
    # Salinity (ppt)  
    sal.sig[i] ~ dnorm(sal.mc, sal.pc)I(0.95,1.01)
    sal[i] = sal[i-1] * sal.sig[i]
    
    # Temp in C
    tempC.sig[i] ~ dnorm(tempC.mc, tempC.pc)
    tempC[i] = tempC[i-1] + tempC.sig[i]
    
    # [Ca] (mmol kg^-1)
    xca.sig[i] ~ dnorm(xca.mc, xca.pc)
    xca[i] = xca[i-1] + xca.sig[i]
    
    # [Mg] (mmol kg^-1)
    xmg.sig[i] ~ dnorm(xmg.mc, xmg.pc)
    xmg[i] = xmg[i-1] + xmg.sig[i]
    
    # [SO4] (mmol kg^-1)
    xso4.sig[i] ~ dnorm(xso4.mc, xso4.pc)
    xso4[i] = xso4[i-1] + xso4.sig[i]
    
    # d11B of seawater (per mille SRM-951) 
    d11Bsw.sig[i] ~ dnorm(d11Bsw.mc, d11Bsw.pc)
    d11Bsw[i] = d11Bsw[i-1] + d11Bsw.sig[i]
    
    # d18O of seawater (per mille SMOW) 
    d18Osw.sig[i] ~ dnorm(d18Osw.mc, d18Osw.pc)
    d18Osw[i] = d18Osw[i-1] + d18Osw.sig[i] 
    
    # Atmospheric pCO2 (atm)
    pco2.sig[i] ~ dnorm(pco2.mc, pco2.pc)I(0.5,1.5)
    pco2[i] = pco2[i-1] * pco2.sig[i]
    
    # DIC (mol kg^-1) - make change in DIC temp dependent (C cycle model)?
    dic.sig[i] ~ dnorm(dic.mc, dic.pc)I(0.99,1.05)
    #dic[i] = dic[i-1] * dic.sig[i]
    dic[i] = dic[i-1]*dic.sig[i] + (tempC[i]-tempC[i-1])*0.00005 
    
  }
  
  # Time independent parameters 
  
  # 'm' boron vital effect calibration parm 
  m.1 ~ dnorm(m.Grub, 1/m.Grubu^2)  
  m.2 ~ dnorm(m.Tsac, 1/m.Tsacu^2)  
  
  # 'c' boron vital effect calibration parm
  c.1 ~ dnorm(c.Grub, 1/c.Grubu^2)     
  c.2 ~ dnorm(c.Tsac, 1/c.Tsacu^2) 
  
  # Pressure (bar)
  press ~ dnorm(10, press.p)    
  press.p = 1/2^2 
  
}



  