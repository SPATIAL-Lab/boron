model{
############################################################################################
#    LIKELIHOOD FUNCTION: evaluates the data against modeled values
############################################################################################    

  # Gaussian precision for d18O measurements (1sd = 0.8 per mille based on long-term standard reproducibility)
  d18Of.p = 1/0.08^2
  
  for (i in 1:length(ai.d11B1)){
  d11Bf.data1[i] ~ dnorm(d11Bf.1[ai.d11B1[i]], d11Bf.p1[i])
  d11Bf.p1[i] = 1/d11Bfu.data1[i]^2 # Gaussian precision for d11Bf (species #1) measurements from se of replicate analyses 
  }
  
  for (i in 1:length(ai.d11B2)){
  d11Bf.data2[i] ~ dnorm(d11Bf.2[ai.d11B2[i]], d11Bf.p2[i])
  d11Bf.p2[i] = 1/d11Bfu.data2[i]^2 # Gaussian precision for d11Bf (species #2) measurements from se of replicate analyses 
  }
  
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
  xcam <- 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm <- 52.8171 # modern [Mg] (mmol kg^-1)
  xso4m <- 28.24  # modern [SO4] (mmol kg^-1)
  mgcaswm <- xmgm/xcam # modern Mg/Ca of seawater

  # Set fractionation factor: Klochko et al. (2006)
  alpha ~ dnorm(alpha.m, alpha.p)     
  alpha.m = 1.0272    #gaussian mean
  alpha.p = 1/0.0003^2   #gaussian precision
  epsilon <- (alpha - 1)*1000  # Compute epsilon from alpha
  
  # Index of diagenetic overprint
  indexop ~ dnorm((seccal/100), (1/(seccal.u/100)^2))
  
  # Sensitivity of seawater d18O to salinity 
  sw.sens ~ dnorm(0.558, 1/0.03^2) # regression slope and uncertainty (i.e., from standard deviation) of GEOSECS obs. reported in Charles and Fairbanks (1990) after Duplessey et al. (1991)

  # Correction for Mg/Ca salinity effect 
  salcorrco ~ dnorm(0.042, (1/0.004^2))
  
  # pH correction on Mg/Ca
  pHcorrco ~ dnorm(pHpccorr, (1/pHpccorrsd^2))T(0,3.8)
  
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
  A ~ dnorm(A.m, A.p)T(0.067,0.084) # Exponential constant in Mg/Ca-SST calibration (Evans et al., 2016)
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
    
    # Compute pCO2 from DIC and pH 
    hyd[i] <- 10^(-(pH[ai.prox[i]]))
    co2[i] <- (dic[ai.prox[i]]) / (1 + (Ks1[i]/hyd[i]) + ((Ks1[i]*Ks2[i])/(hyd[i]^2)))
    fco2[i] <- co2[i] / K0[i]
    pco2[i] <- fco2[i] / 0.9968
    
    # Compute d11Bforam from pH and d11Bsw
    pKsB[i] <- -(log(KsB[i])/log(10))
    t1[i] <- 10^(pKsB[i]-pH[ai.prox[i]])
    d11Bb[i] <- ((t1[i]*epsilon)-(t1[i]*d11Bsw[ai.prox[i]])-d11Bsw[ai.prox[i]])/(-((t1[i]*alpha)+1))
    
    c.final.1[i] <- c.1 + c.correction1 # adjusts final c value if desired - correction specified in driver
    d11Bf.1[i] <- m.1*d11Bb[i]+ (c.final.1[i])
    
    c.final.2[i] <- c.2 + c.correction2 # adjusts final c value if desired - correction specified in driver
    d11Bf.2[i] <- m.2*d11Bb[i]+ (c.final.2[i])     
    
    # Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
    d18Osw.sc[i] <- d18Osw[ai.prox[i]] + (sw.sens*(sal[ai.prox[i]]-35))
    d18Oswpdb[i] <- d18Osw.sc[i] -0.27
    d18Of.pr[i] <- d18Oswpdb[i] + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC[ai.prox[i]])))^0.5))/(2*0.09))
    d18Of[i] <- d18Of.pr[i]*(1-indexop) + indexop*d18Oseccal

    # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
    mgcasw[i] <- (xmg[ai.prox[i]]/xca[ai.prox[i]])     
    Bcorr[i] <- ((mgcasw[i]^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr[i] <- Bcorr[i]*(exp(A*tempC[ai.prox[i]]))
    mgca_sal[i] <- mgca_corr[i] / (1-(8.05-pH[ai.prox[i]])*pHcorrco)
    mgcaf[i] <- mgca_sal[i] / (1-(sal[ai.prox[i]]-35)*salcorrco)
    
  }
  
############################################################################################
#    ENVIRONMENTAL MODEL
############################################################################################    
 
  # Environmental time-dependent prior initial conditions 
        # .phi = temporal autocorrelation of a parameter
        # .eps =  error term 
        # .tau = error precision for dt = 1
        # .pc = error precision of temporal autocorrelation error term (.eps)

  
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
  
  # [SO4] (mmol kg^-1)
  xso4[1] ~ dnorm(xso4.m, xso4.p)      
  xso4.phi ~ dbeta(5,2) 
  xso4.eps[1] = 0 
  xso4.tau ~ dgamma(1e3, 4e-4)
  
  # d11B of seawater (per mille SRM-951) 
  d11Bsw[1] ~ dnorm(d11Bsw.m, d11Bsw.p)    
  d11Bsw.phi ~ dbeta(5,2)  
  d11Bsw.eps[1] = 0 
  d11Bsw.tau ~ dgamma(4e3, 1e-3)
  
  # d18O of seawater (per mille SMOW) 
  d18Osw[1] ~ dnorm(d18Osw.m, d18Osw.p)    
  d18Osw.phi ~ dbeta(5,2) 
  d18Osw.eps[1] = 0 
  d18Osw.tau ~ dgamma(1e3, 1e-3)
  
  # pH (total scale)
  pH[1] ~ dunif(pH.l, pH.u)   
  pH.phi ~ dbeta(5,2)
  pH.eps[1] = 0 
  pH.tau ~ dgamma(3000, 1e-1) 
  
  # DIC (mmol kg^-1) 
  dic[1] ~ dnorm(dic.sim[1], dic.p)
  
  
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
    
    # [Ca] (mmol kg^-1); linear decline in Cenozoic follows Holland et al. (2020). 20 mmol/kg decline over last 120 Myr
    xca.tdep[i] = -0.00019 * (ages.prox[1]-ages.prox[i])  
    xca.pc[i] = xca.tau*((1-xca.phi^2)/(1-xca.phi^(2*dt[i-1])))
    xca.eps[i] ~ dnorm(xca.eps[i-1]*(xca.phi^dt[i-1]), xca.pc[i])T(-0.3, 0.3)
    xca[i] = xca[1] * (1 + xca.eps[i]) + xca.tdep[i]
    
    # [Mg] (mmol kg^-1); imposing max realistic decline in Mg/Casw; decrease suggested by paired benthic Mg/Ca+d18O from ODP Site 1209 (data suggest even greater decline than what's imposed here; Mg decrease used here is derived from max observed shifts in Mg/Casw over Cenozoic as in Holland et al. (2020)
    xmg.tdep[i] = -0.00274 * (ages.prox[1]-ages.prox[i])    
    xmg.pc[i] = xmg.tau*((1-xmg.phi^2)/(1-xmg.phi^(2*dt[i-1])))
    xmg.eps[i] ~ dnorm(xmg.eps[i-1]*(xmg.phi^dt[i-1]), xmg.pc[i])T(-0.3, 0.3)
    xmg[i] = xmg[1] * (1 + xmg.eps[i]) + xmg.tdep[i]
    
    # [SO4] (mmol kg^-1)
    xso4.pc[i] = xso4.tau*((1-xso4.phi^2)/(1-xso4.phi^(2*dt[i-1])))
    xso4.eps[i] ~ dnorm(xso4.eps[i-1]*(xso4.phi^dt[i-1]), xso4.pc[i])T(-0.3, 0.3)
    xso4[i] = xso4[1] * (1 + xso4.eps[i])
    
    # d11B of seawater (per mille SRM-951) 
    d11Bsw.pc[i] = d11Bsw.tau*((1-d11Bsw.phi^2)/(1-d11Bsw.phi^(2*dt[i-1])))
    d11Bsw.eps[i] ~ dnorm(d11Bsw.eps[i-1]*(d11Bsw.phi^dt[i-1]), d11Bsw.pc[i])T(-1, 1)
    d11Bsw[i] = d11Bsw[1] + d11Bsw.eps[i]
    
    # d18O of seawater (per mille SMOW) 
    d18Osw.pc[i] = d18Osw.tau*((1-d18Osw.phi^2)/(1-d18Osw.phi^(2*dt[i-1])))
    d18Osw.eps[i] ~ dnorm(d18Osw.eps[i-1]*(d18Osw.phi^dt[i-1]), d18Osw.pc[i])T(-2, 2)
    d18Osw[i] = d18Osw[1] + d18Osw.eps[i]
    
    # pH
    pH.pc[i] = pH.tau*((1-pH.phi^2)/(1-pH.phi^(2*dt[i-1])))
    pH.eps[i] ~ dnorm(pH.eps[i-1]*(pH.phi^dt[i-1]), pH.pc[i])T(-1, 1)
    pH[i] = pH[1] + pH.eps[i]
    
    # DIC (mol kg^-1) - change in DIC prior based on LOSCAR output
    dic[i] ~ dnorm(dic.sim[i], dic.p)
    
  }
  
  # Time independent parameters 
  
    # 'm' boron vital effect calibration parm 
    m.1 ~ dnorm(m.Grub, 1/m.Grubu^2)  
    m.2 ~ dnorm(m.Tsac, 1/m.Tsacu^2)  
    
    # 'c' boron vital effect calibration parm
    c.1 ~ dnorm(c.Grub, 1/c.Grubu^2)     
    c.2 ~ dnorm(c.Tsac, 1/c.Tsacu^2) 
    
    # Pressure (bar)
    press ~ dnorm(6, press.p)    
    press.p = 1/1^2 
    
}


  