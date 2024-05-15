
model{
  ############################################################################################
  #    LIKELIHOOD FUNCTION: evaluates the data against modeled values
  ############################################################################################    

    d11Bf.data1 ~ dnorm(d11Bf.1, d11Bf.p1)
    d11Bf.p1 = 1/d11Bfu.data1^2 

    d11Bf.data2 ~ dnorm(d11Bf.2, d11Bf.p2)
    d11Bf.p2 = 1/d11Bfu.data2^2 
    
    mgcaf.data ~ dnorm(mgcaf, mgcaf.p)
    mgcaf.p = 1/(mgcafu.data)^2  
  
    d18Of.data ~ dnorm(d18Of, d18Of.p)
    d18Of.p = 1/0.08^2
  
  
  ############################################################################################
  #    PROXY SYSTEM MODEL
  ############################################################################################    
  
  
  # INPUT VALUES FOR CALCULATIONS 
  # Includes time invariant prior distributions of parameters and fixed parameters 
  
  # Set modern concentrations for Mg, Ca, and SO4
  xcam = 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm = 52.8171 # modern [Mg] (mmol kg^-1)
  xso4m = 28.24  # modern [SO4] (mmol kg^-1)
  mgcaswm <- xmgm/xcam # modern Mg/Ca of seawater
  
  # Set fractionation factor: Klochko et al. (2006)
  alpha ~ dnorm(alpha.m, alpha.p)     
  alpha.m = 1.0272    #gaussian mean
  alpha.p = 1/0.0003^2   #gaussian precision
  epsilon <- (alpha - 1)*1000  # Compute epsilon from alpha
  
  # Index of diagenetic overprint
  indexop ~ dnorm((seccal/100), (1/(seccal.u/100)^2))T(0,1)
  
  # Sensitivity of seawater d18O to salinity 
  # Regression slope and uncertainty (i.e., from standard deviation) of GEOSECS obs. reported in 
  # Charles and Fairbanks (1990) after Duplessey et al. (1991)
  sw.sens ~ dnorm(0.558, 1/0.03^2) 
  
  # Correction for d18O sea surface pH effect 
  d18O_pHcorr ~ dnorm(d18O_pHcorr.avg, 1/d18O_pHcorr.sd^2)T(-2,0) 
  
  # Distributions for coefficients in calcite - water O isotope fractionation (Daeron et al., 2019)
  A_Daeron ~ dnorm(17.57, 1/0.43^2)T(16.71, 18.43)
  B_Daeron ~ dnorm(29.89, 1/0.06^2)T(29.77, 30.01)
  
  # Correction for Mg/Ca salinity effect 
  salcorrco ~ dnorm(0.042, (1/0.004^2))
  
  # Correction for Mg/Ca D[CO3=]bottom water effect (i.e., below 21.3 umol/kg)
  mgca_satdec ~ dnorm(0.054, 1/0.008^2)
  
  # Correction for Mg/Ca surface pH effect (Hollis et al., 2019; Gray et al., 2016; Gray and Evans et al., 2019)
  mgca_pHcorr ~ dnorm((mgca_pHcorr.avg+1e-20)/10, (1/(mgca_pHcorr.sd+1e-20)^2)/10)T(0,1) 
  
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
  
  R = 83.131 # constant (cm^3 bar mol^-1 K^-1)
  

    # CARB CHEM EQUILIBRIUM CONSTANT CALCULATIONS FOLLOWING ZEEBE AND TYRRELL (2019)
    
    # Calculate equil. constants using salinity and temp:
    temp <- tempC+273.15
    Ks1m_st <-exp(2.83655-2307.1266/temp-1.5529413*(log(temp))-((0.20760841+4.0484/temp)*sqrt(sal))+0.0846834*sal-0.00654208*(sal^1.5)+log(1-(0.001005*sal)))
    Ks2m_st <- exp(-9.226508-3351.6106/temp-0.2005743*(log(temp))-((0.106901773+23.9722/temp)*sqrt(sal))+0.1130822*sal-0.00846934*(sal^1.5)+log(1-(0.001005*sal)))
    logKsspcm_st <- ((-171.9065-0.077993*temp+2839.319/temp+71.595*(log(temp)/log(10))+(-0.77712+0.0028426*temp+178.34/temp)*(sal^0.5)-0.07711*sal+0.0041249*(sal^1.5)))
    Ksspcm_st <- 10^(logKsspcm_st)
    lnKsB_st <- ((-8966.9-2890.53*sal^0.5-77.942*sal+1.728*sal^1.5-0.0996*sal^2)/temp)+148.0248+137.1942*sal^0.5+1.62142*sal-(24.4344+25.085*sal^0.5+0.2474*sal)*(log(temp))+(0.053105*sal^0.5*temp)
    KsB_st <- exp(lnKsB_st)
    Ksw_st <- exp(148.96502-13847.26/temp-23.6521*(log(temp))+(118.67/temp-5.977+1.0495*(log(temp)))*(sal^0.5)-0.01615*sal)
    K0 <- exp(9345.17/temp-60.2409+23.3585*(log(temp/100))+sal*(0.023517-0.00023656*temp+0.0047036*((temp/100)^2)))
    
    # Adjust equil. constants for the effect of pressure (Millero 1995):
    delV1 <- (-25.50)+0.1271*tempC
    delV2 <- (-15.82)+(-0.0219*tempC)
    delVspc<- (-48.76)+(0.5304*tempC)
    delVB <- (-29.48)+0.1622*tempC+(2.608/1000)*tempC^2
    delVw <- (-25.60)+0.2324*tempC+(-3.6246/1000)*tempC^2
    
    delk1 <- (-3.08/1000)+(0.0877/1000)*tempC
    delk2 <- (1.13/1000)+(-0.1475/1000)*tempC
    delkspc <- (-11.76/1000)+(0.3692/1000)*tempC
    delkB <- -2.84/1000
    delkw <- (-5.13/1000)+(0.0794/1000)*tempC
    
    Ks1m <- (exp(-((delV1/(R*temp))*press)+((0.5*delk1)/(R*temp))*press^2))*Ks1m_st
    Ks2m <- (exp(-((delV2/(R*temp))*press)+((0.5*delk2)/(R*temp))*press^2))*Ks2m_st
    Ksspcm <- (exp(-((delVspc/(R*temp))*press)+((0.5*delkspc)/(R*temp))*press^2))*Ksspcm_st
    KsB <- (exp(-((delVB/(R*temp))*press)+((0.5*delkB)/(R*temp))*press^2))*KsB_st
    Ksw <- (exp(-((delVw/(R*temp))*press)+((0.5*delkw)/(R*temp))*press^2))*Ksw_st
    
    # K*1, K*2, and K*spc are corrected for past seawater [Ca], [Mg], and [SO4] following ZT19
    Ks1 <- Ks1m*(1+(s1_ca*(xca/xcam-1)+s1_mg*(xmg/xmgm-1)+s1_so4*(xso4/xso4m-1)))
    Ks2 <- Ks2m*(1+(s2_ca*(xca/xcam-1)+s2_mg*(xmg/xmgm-1)+s2_so4*(xso4/xso4m-1)))
    Ksspc <- Ksspcm*(1+(sspc_ca*(xca/xcam-1)+sspc_mg*(xmg/xmgm-1)+sspc_so4*(xso4/xso4m-1)))
    
    
    # DETERMINE FORAMINIFERAL D18O, D11B, AND MG/CA 
    
    # Compute pCO2 from DIC and pH 
    hyd <- 10^(-(pH))
    co2 <- (dic) / (1 + (Ks1/hyd) + ((Ks1*Ks2)/(hyd^2)))
    fco2 <- co2 / K0
    pco2 <- fco2 / 0.9968
    
    # Compute d11Bforam from pH and d11Bsw
    pKsB <- -(log(KsB)/log(10))
    t1 <- 10^(pKsB-pH)
    d11Bb <- ((t1*epsilon)-(t1*d11Bsw)-d11Bsw)/(-((t1*alpha)+1))
    
    c.final.1 <- c.1 + c.correction1 # adjusts final c value if desired - correction specified in driver
    d11Bf.1 <- m.1*d11Bb+ (c.final.1)
    
    c.final.2 <- c.2 + c.correction2 # adjusts final c value if desired - correction specified in driver
    d11Bf.2 <- m.2*d11Bb+ (c.final.2)     
    
    # Compute d18Oforam (Daeron et al., 2019; Hollis et al., 2019)
    d18Osw.sc <- d18Osw + (sw.sens*(sal-35))
    alpha.ccw <- exp((A_Daeron*10^3*(1/temp - 1/297.7) + B_Daeron) / 10^3)
    epsilon.ccw <- (alpha.ccw-1)*10^3 
    d18Ofsmow <- d18Osw.sc + epsilon.ccw 
    d18Ofpdb <- 0.97001 * d18Ofsmow - 29.99
    d18Ofdiag <- d18Ofpdb - (pH - 8)*d18O_pHcorr
    d18Of <- d18Ofdiag*(1-indexop) + d18Oseccal*indexop
    
    # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
    mgcasw <- (xmg/xca)     
    Bcorr <- ((mgcasw^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr <- Bcorr*(exp(A*tempC))
    mgca_corr2 <- mgca_corr / (1-(sal-35)*salcorrco)
    mgca_corr3 <- mgca_corr2 / (1-(pH-8.05)*mgca_pHcorr)
    mgcaf <- mgca_corr3 + DDco3_corr*mgca_satdec

  
  ############################################################################################
  #    ENVIRONMENTAL MODEL
  ############################################################################################    
  
  # Salinity (ppt)  
  sal ~ dnorm(sal.m, sal.p)T(30,40)  
  
  # Temp in C
  tempC ~ dnorm(tempC.m, tempC.p)T(24,36) 
  
  # [Ca] (mmol kg^-1)
  xca ~ dnorm(xca.m, xca.p)

  # [Mg] (mmol kg^-1)
  xmg ~ dnorm(xmg.m, xmg.p)      
 
  # [SO4] (mmol kg^-1)
  xso4 ~ dnorm(xso4.m, xso4.p)      
 
  # d11B of seawater (per mille SRM-951) 
  d11Bsw ~ dnorm(d11Bsw.m, d11Bsw.p)    

  # d18O of seawater (per mille SMOW) 
  d18Osw ~ dnorm(d18Osw.m, d18Osw.p)T(-1,-0.5)    
 
  # pH (total scale)
  pH ~ dunif(pH.l, pH.u)   
  
  # DIC (mmol kg^-1) 
  dic ~ dnorm(dic.m, 1/dic.sd^2)
  
  # D[CO3=] of bottom water (umol kg^-1) difference from dissolution susceptibility limit of 21.3
  DDco3_corr ~ dnorm(DDco3.m, 1/DDco3.sd^2)

  # 'm' boron vital effect calibration parm 
  m.1 ~ dnorm(m.Grub, 1/m.Grubu^2)  
  m.2 ~ dnorm(m.Tsac, 1/m.Tsacu^2)  
  
  # 'c' boron vital effect calibration parm
  c.1 ~ dnorm(c.Grub, 1/c.Grubu^2)     
  c.2 ~ dnorm(c.Tsac, 1/c.Tsacu^2) 
  
  # Pressure (bar)
  press ~ dnorm(6, 1/1^2)    

}
