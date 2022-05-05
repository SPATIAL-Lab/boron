model{
  
  #GJB - Likelihood function: this block evaluates the data against
  #the modeled values
  mgcaf.data ~ dnorm(mgcaf, mgcaf.p)
  mgcaf.p = 1/0.03^2 #GJB - I've made up measurement uncertainties
  d18Of.data ~ dnorm(d18Of, d18Of.p)
  d18Of.p = 1/0.1^2
  d11Bf.data ~ dnorm(d11Bf, d11Bf.p)
  d11Bf.p = 1/0.2^2
  
  # GJB - Proxy System Model
  
  # INPUT FOR MODERN VALUES USED, CORRECTIONS APPLIED, CALIBRATION APPROACH 
  
  # Set modern concentrations for Mg, Ca, and SO4
  xcam <- 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm <- 52.8171 # modern [Mg] (mmol kg^-1)
  xso4m <- 28.24  # modern [SO4] (mmol kg^-1)
  mgcaswm <- xmgm/xcam # modern Mg/Ca of seawater
  
  # Set the vital effect correction (i.e., d11B of borate to d11B of foram)
  m <- 0.88
  c <- 1.73 #- 0.8  # adjusted to align to d11B PETM data 
  
  # Set fractionation factor 
  alpha <- 1.0272              # Klochko et al. (2006)
  epsilon <- (alpha - 1)*1000  # Compute epsilon from alpha
  
  # GJB - I'm leaving the multiple approaches out for now...
  # Select if apporach follows (0) Holland et al., 2020 (DIC-corrected), 
  # OR (1) Hollis et al., 2019 (pH-corrected) for Mg/Ca-SST
#  sstapp <- 0
  
  # Holland et al., 2020 input: A (Ap), B(Bp), C(Cp), D(Dp), E(Ep) species/seawater specific parameters 
  # Default input follows O. universa calibration (ie Fig 4A) following paleogene application in paper
  # See paper for calibration uncertainties
#  Ap <- 0.785
#  Bp <- 0.338
#  Cp <- 22.762
#  Dp <- 0.087
#  Ep <- -0.27 # E modified following Holland et al., 2020 PETM example calc
  
  # Hollis et al., 2019 input: set H, B and A values (i.e., Evans et al., 2012, 2016b)
  # See paper for calibration uncertainties 
  Hp <- 0.41    # nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer)
  Bmod <- 0.38  # modern pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003) 
  A <- 0.09     # Exponenital constant in Mg/Ca-SST calibration (Anand et al., 2003)
  
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
  d11Bf <- m*d11Bb + c
  
  # Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
  d18Oswpdb <- d18Osw -0.27
  d18Of <- d18Oswpdb + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC)))^0.5))/(2*0.09))
  
  # Compute d18Oforam using salintiy-derived d18Osw (Anand et al., 2003)
#  d18Oswsal <- (-19.264+(0.558*sal))-0.27 # includes conversion factor (0.27) for SMOW to PDB
#  d18Ofsal <- d18Oswsal + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC)))^0.5))/(2*0.09))
  
  # Compute Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
#  if(sstapp > 0){
    Bcorr <- ((mgcasw^Hp)/(mgcaswm^Hp)) * Bmod
    mgca_corr <- Bcorr*(exp(A*tempC))
    mgca_sal <- mgca_corr / (1-(8.05-pH)*0.070)  
    mgcaf <- mgca_sal / (1-(sal-35)*0.042)
#    }
  
  # Compute Mg/Caforam following Holland et al. (2020) Mg/Ca carb chem correction approach
#  if(sstapp < 1){
#    mgcaf <- mgcasw^Ap*dic^Bp*(exp(Cp*xca+Dp*tempC+Ep))}
  
  
  # GJB - These are the prior distributions for the stochastic
  # model parameters. In this case they're all environmental
  # parms 'external' to the psm, but we could also add some of 
  # the PSM parameters as appropriate. Using gamma, normal, and 
  # uniform distributions for different data types...
    
  # Environmental priors
    
  sal ~ dgamma(sal.s, sal.r)       # ppt
  sal.s = 35  #gamma shape parameter
  sal.r = 1   #gamma rate parameter
  tempC ~ dnorm(temp.m, temp.p)     # temp in C
  temp.m = 18   #gaussian mean
  temp.p = 1/5^2   #gaussian precision
  temp <- tempC+273.15        # temp in K
  press ~ dgamma(press.s, press.r)     # bar
  press.s = 1
  press.r = 10
  xca ~ dgamma(xca.s, xca.r)       # [Ca] (mol kg^-1)
  xca.s = 17
  xca.r = 1e3
  xmg ~ dgamma(xmg.s, xmg.r)       # [Mg] (mol kg^-1)
  xmg.s = 36
  xmg.r = 1e3
  mgcasw <- (xmg/xca)         # Seawater Mg/Ca 
  xso4 ~ dgamma(xso4.s, xso4.r)      # [SO4] (mmol kg^-1)
  xso4.s = 28
  xso4.r = 2
  d11Bsw ~ dnorm(d11Bsw.m, d11Bsw.p)    # d11B of seawater (per mille SRM-951) 
  d11Bsw.m = 38.45
  d11Bsw.p = 1/0.5^2
  d18Osw ~ dnorm(d18Osw.m, d18Osw.p)    # d18O of seawater (per mille SMOW) 
  d18Osw.m = -1
  d18Osw.p = 1/0.5^2
  pco2 ~ dunif(pco2.l, pco2.u)     # atmospheric pCO2 (atm)
  pco2.l = 500e-6   #uniform lower bound
  pco2.u = 1500e-6   #uniform upper bound
  dic ~ dunif(dic.l, dic.u)      # seawater [CO3] (mol kg^-1)
  dic.l = 0.0015
  dic.u = 0.003
}