model{
  
  #GJB - Likelihood function: this block evaluates the data against
  #the modeled values
  mgcaf.data ~ dnorm(mgcaf, mgcaf.p)
  mgcaf.p = 1/0.03^2 
  d18Of.data ~ dnorm(d18Of, d18Of.p)
  d18Of.p = 1/0.1^2
  d11Bf.data ~ dnorm(d11Bf, d11Bf.p)
  d11Bf.p = 1/d11Bfu.data^2

  # Determine 'm' and 'c' d11Bforam-d11Bborate calibration parameters using modern species calibration points. Data is specified in 'driver'

  for(i in 1:length(d11Bcb)){
    d11Bcb[i] ~ dnorm(d11Bcb.m[i], 1/0.2^2) # second term is average measurement precision for calibration points
    d11Bcb.m[i] = d11Bcfo[i]*m + c
  }

  # GJB - Proxy System Model
  
  # INPUT FOR MODERN VALUES USED, CORRECTIONS APPLIED, CALIBRATION APPROACH 
  
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

  
  # Hollis et al., 2019 input: set H, B and A values (i.e., Evans et al., 2012, 2016b)
  # DH - Includes calibration uncertainty in these terms
  Hp ~ dnorm(Hp.m, Hp.p) # nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer is typically used = 0.41)
  Hp.m = 0.6
  Hp.p = 1/0.1^2
  Bmod ~ dnorm(Bmod.m, Bmod.p) # modern pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
  Bmod.m = 0.38  
  Bmod.p = 1/0.02^2
  A ~ dnorm(A.m, A.p) # Exponenital constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans )
  A.m = 0.07    
  A.p = 1/0.01^2
  
  
  #################################################################################################################################
  
  
  #################################################################################################################################
  # CARB CHEM EQUILIBRIUM CONSTANT CALCULATIONS FOLLOWING ZEEBE AND TYRRELL (2019)
  
  # GJB - note, no log10 function in JAGS, using log(x)/log(10)
  # Calculate equil. constants using salinity and temp:
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
  delVspc <- (-48.76)+(0.5304*tempC)
  delVB <- (-29.48)+0.1622*tempC+(2.608/1000)*tempC^2
  delVw <- (-25.60)+0.2324*tempC+(-3.6246/1000)*tempC^2
  
  delk1 <- (-3.08/1000)+(0.0877/1000)*tempC
  delk2 <- (1.13/1000)+(-0.1475/1000)*tempC
  delkspc <- (-11.76/1000)+(0.3692/1000)*tempC
  delkB <- -2.84/1000
  delkw <- (-5.13/1000)+(0.0794/1000)*tempC
  
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
  Ks1 <- Ks1m*(1+(s1_ca*(xca/xcam-1)+s1_mg*(xmg/xmgm-1)+s1_so4*(xso4/xso4m-1)))
  Ks2 <- Ks2m*(1+(s2_ca*(xca/xcam-1)+s2_mg*(xmg/xmgm-1)+s2_so4*(xso4/xso4m-1)))
  Ksspc <- Ksspcm*(1+(sspc_ca*(xca/xcam-1)+sspc_mg*(xmg/xmgm-1)+sspc_so4*(xso4/xso4m-1)))

  #################################################################################################################################
  
  #################################################################################################################################
  # DETERMINE FORAMINIFERAL D18O, D11B, AND MG/CA 
  
  # GJB - In my experiments this is likley to be the biggest source
  # of instability...I think for certain combos of pco2 and dic the 
  # roots h1 and h2 are undefined?
  
  # Compute pH from [co2] and DIC
  fco2 <- pco2*0.9968
  co2 <- fco2*K0
  qa <- co2-dic
  qb <- co2*Ks1
  qc <- co2*Ks1*Ks2
  h1 <- (-qb+sqrt((qb^2)-(4*qa*qc)))/(2*qa)
  h2 <- (-qb-sqrt((qb^2)-(4*qa*qc)))/(2*qa)
  pH <- -(log(max(h1, h2))/log(10))                      
  
  # Compute d11Bforam from pH and d11Bsw
  pKsB <- -(log(KsB)/log(10))
  t1 <- 10^(pKsB-pH)
  d11Bb <- ((t1*epsilon)-(t1*d11Bsw)-d11Bsw)/(-((t1*alpha)+1))
  c.final <- c + c.correction # adjusts final c value if desired - correction specified in driver
  d11Bf <- m*d11Bb + (c.final)
   
  # Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
  sw.sens ~ dnorm(0.558, 1/0.03^2) # Uncertainty represents std dev of regression slope in GEOSECS obs. reported in Charles and Fairbanks (1990)
  d18Osw.sc <- d18Osw + (sw.sens*(sal-35)) # update d18Osw to reflect changes in salinity
  d18Oswpdb <- d18Osw.sc -0.27
  d18Of.pr <- d18Oswpdb + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC)))^0.5))/(2*0.09))
  indexop ~ dnorm(indexop.m, indexop.p)
  indexop.m = seccal/100
  indexop.p = 1/(seccal.u/100)^2
  d18Of <- d18Of.pr + indexop*Dd18Oseccal
  
  # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
    Bcorr <- ((mgcasw^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr <- Bcorr*(exp(A*tempC))
    pHcorrco ~ dnorm(pHcorrco.m, pHcorrco.p)
    pHcorrco.m = 0.70
    pHcorrco.p = 1/0.09^2
    mgca_sal <- mgca_corr / (1-(8.05-pH)*pHcorrco)
    salcorrco ~ dnorm(salcorrco.m, salcorrco.p)
    salcorrco.m = 0.042
    salcorrco.p = 1/0.004^2
    mgcaf <- mgca_sal / (1-(sal-35)*salcorrco)

  
  # GJB - These are the prior distributions for the stochastic
  # model parameters. In this case they're all environmental
  # parms 'external' to the psm, but we could also add some of 
  # the PSM parameters as appropriate. Using gamma, normal, and 
  # uniform distributions for different data types...
    
  # Environmental priors
    
  sal ~ dnorm(sal.m, sal.p)     # temp in C
  sal.m = 35   #gaussian mean
  sal.p = 1/0.5^2   #gaussian precision  
  tempC ~ dnorm(temp.m, temp.p)     # temp in C
  temp.m = 30   #gaussian mean
  temp.p = 1/10^2   #gaussian precision
  temp <- tempC+273.15        # temp in K
  press ~ dnorm(press.m, press.p)     # bar
  press.m = 10
  press.p = 1/2^2
  xca ~ dnorm(xca.m, xca.p)       # [Ca] (mol kg^-1)
  xca.m = 17
  xca.p = 1/0.5^2
  xmg ~ dnorm(xmg.m, xmg.p)       # [Mg] (mol kg^-1)
  xmg.m = 36
  xmg.p = 1/0.5^2
  mgcasw <- (xmg/xca)         # Seawater Mg/Ca 
  xso4 ~ dnorm(xso4.m, xso4.p)      # [SO4] (mmol kg^-1)
  xso4.m = 14
  xso4.p = 1/0.5^2
  d11Bsw ~ dnorm(d11Bsw.m, d11Bsw.p)    # d11B of seawater (per mille SRM-951) 
  d11Bsw.m = 38.45
  d11Bsw.p = 1/0.5^2
  d18Osw ~ dnorm(d18Osw.m, d18Osw.p)    # d18O of seawater (per mille SMOW) 
  d18Osw.m = -1
  d18Osw.p = 1/0.5^2
  pco2 ~ dunif(pco2.l, pco2.u)     # atmospheric pCO2 (atm)
  pco2.l = 300e-6   #uniform lower bound
  pco2.u = 2500e-6   #uniform upper bound
  dic ~ dunif(dic.l, dic.u)      # seawater DIC (mol kg^-1)
  dic.l = 0.0018
  dic.u = 0.0024
  
  # "Vital effect" priors, posterior values are used in PSM calcs
  m ~ dnorm(m.mean, 1/m.sd^2)    # prior distribution of m
  c ~ dnorm(c.mean, 1/c.sd^2)    # prior distribution of c

}


  