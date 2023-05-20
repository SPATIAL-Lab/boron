
# Driver for time series data inversion using forward foraminifera PSM 
# Reads in time series proxy data and LOSCAR-simulated DIC data for DIC priors at each time step
# This version includes generation of look-up tables for carb chem equil constant calcs
# 
# Dustin T. Harper
# 11 May 2023


############################################################################################
# Load libraries 
library(rjags)
############################################################################################

############################################################################################  
# jags specifications
n.chains = 3
n.iter = 200
n.burnin = 50
n.thin = 1

############################################################################################  
# pH correction on Mg/Ca 

# Input pH correction on Mg/Ca as the percent change per 0.1 pH unit 
# e.g., Hollis et al., 2019 recommend 7% per 0.1 pH unit, plus 1sd unc = 0.9%; Haynes et al., 2023 suggest minimal to zero effect at low Mg/Casw
pHpccorr = 2 
pHpccorrsd = 0.9


############################################################################################  
# Diagenesis correction on d18O

# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 60    # Percent recrystallized; DATA FOR JAGS
seccal.u = 2.5      # 1sd % recrystallized; DATA FOR JAGS
d18Oseccal = 0.85 # Calculated following Edgar et al. (2015); DATA FOR JAGS


############################################################################################  
# Mg/Ca SST proxy vital effects and calibration parameters 

# Nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer) 
Hp.mean = 0.41      
Hp.sd = 0.1         

# Modern (pre-corrected) pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
Bmod.mean = 0.38    
Bmod.sd = 0.02      

# Exponential constant in Mg/Ca-SST calibration (Evans et al., 2016)
A.mean = 0.0757      
A.sd = 0.0045        


############################################################################################  
# d11Bforam-d11Bborate "vital effect" calibration

# G. ruber vitals
m.Grub = 0.62       # mean "m" value for G. ruber distribution 
m.Grubu = 0.055     # s.d. for "m" value for G. ruber distribution
c.Grub = 9.52       # mean "c" value for G. ruber distribution 
c.Grubu = 1.01      # s.d. for "c" value for G. ruber distribution
c.Grub.corr = -3.76 # optional "c" correction for paleo applications

# T. sacculifer vitals
m.Tsac = 0.82       # mean "m" value for T. sacculifer distribution 
m.Tsacu = 0.11      # s.d. for "m" value for T. sacculifer distribution
c.Tsac = 3.94       # mean "c" value for T. sacculifer distribution 
c.Tsacu = 2.01      # s.d. for "c" value for T. sacculifer distribution
c.Tsac.corr = -1.9  # optional "c" correction for paleo applications

# # O. universa vitals
m.Ouni = 0.95       # mean "m" value for O. universa distribution 
m.Ouniu = 0.085     # s.d. for "c" value for O. universa distribution
c.Ouni = -0.42      # mean "c" value for O. universa distribution 
c.Ouniu = 1.43      # s.d. for "c" value for O. universa distribution
c.Ouni.corr = 0     # optional "c" correction for paleo applications

# Custom vitals
m.custom = 0.75     # mean "m" value for custom distribution 
m.customu = 0.1     # s.d. for "c" value for custom distribution
c.custom = 2        # mean "c" value for custom distribution 
c.customu = 1       # s.d. for "c" value for custom distribution

# Borate "vitals" defaults to m=1, c=0 with no uncertainty (i.e., d11Bb = d11Bf exactly)
m.boru = 1e-3       # presribed uncertainty in "borate vitals" m=1 
c.boru = 1e-3       # presribed uncertainty in "borate vitals" c=0 
c.bor.corr = 0      # optional "c" correction for paleo applications (corrects "borate vitals" c=0)

# Vital effect vectors to pass
m.mean <- c(m.Grub, m.Tsac, m.Ouni, m.custom, 1)
m.prec <- c(m.Grubu, m.Tsacu, m.Ouniu, m.customu, m.boru)

c.mean <- c(c.Grub+c.Grub.corr, c.Tsac+c.Tsac.corr, c.Ouni+c.Ouni.corr, c.custom, c.bor.corr)
c.prec <- c(c.Grubu, c.Tsacu, c.Ouniu, c.customu, c.boru)


############################################################################################
#    GENERATE CARBONATE CHEMISTRY EQUIL CONSTANT LOOKUP ARRAYS - follows Zeebe & Wolf-Gladrow 2001
############################################################################################    

# Set upper and lower STP bounds for equil constant array 
tempC.lb = 0
tempC.ub = 65
sal.lb = 15
sal.ub = 60
press.lb = 0
press.ub= 50

# Step increments for sal (ppt) temp (degrees C) and press (bar)
t.inc = 5 #0.25
s.inc = 5 #0.25
p.inc = 5 #1

# Ranges of variables over which to evaluate
tempC.vr = seq(tempC.lb, tempC.ub, by=t.inc)
sal.vr = seq(sal.lb, sal.ub, by=s.inc)
press.vr = seq(press.lb, press.ub, by=p.inc)

# Initiate arrays 
temp.vr = c(1:length(tempC.vr))
delV1 = c(1:length(tempC.vr))
delV2 = c(1:length(tempC.vr))
delVspc = c(1:length(tempC.vr))
delVB = c(1:length(tempC.vr))
delVw = c(1:length(tempC.vr))

delk1 = c(1:length(tempC.vr))
delk2 = c(1:length(tempC.vr))
delkspc = c(1:length(tempC.vr))
delkB = c(1:length(tempC.vr))
delkw = c(1:length(tempC.vr))

base2Darray = c(1:(length(tempC.vr)*length(sal.vr)))
dim(base2Darray) = c((length(tempC.vr)), (length(sal.vr)))

Ks1m_st = base2Darray
Ks2m_st =  base2Darray
logKsspcm_st = base2Darray 
Ksspcm_st = base2Darray
lnKsB_st = base2Darray
KsB_st = base2Darray
Ksw_st = base2Darray
K0a = base2Darray

base3Darray = c(1:(length(tempC.vr)*length(sal.vr)*length(press.vr)))
dim(base3Darray) = c((length(tempC.vr)), (length(sal.vr)), (length(press.vr)))      

Ks1a = base3Darray
Ks2a = base3Darray
Ksspca = base3Darray
KsBa = base3Darray
Kswa = base3Darray

# Constant (cm^3 bar mol^-1 K^-1)
R <- 83.131 

# Nested for loops to calculate 3D array for K*1a (i.e., array of potential K*1 values w/o any major ion effect)
for (i in 1:length(tempC.vr)){
  for (j in 1:length(sal.vr)){
    for (k in 1:length(press.vr)){
      
      temp.vr[i] <- tempC.vr[i]+273.15
      
      Ks1m_st[i,j] <- exp(2.83655-2307.1266/temp.vr[i]-1.5529413*(log(temp.vr[i]))-((0.20760841+4.0484/temp.vr[i])*sqrt(sal.vr[j]))+0.0846834*sal.vr[j]-0.00654208*(sal.vr[j]^1.5)+log(1-(0.001005*sal.vr[j])))
      Ks2m_st[i,j] <- exp(-9.226508-3351.6106/temp.vr[i]-0.2005743*(log(temp.vr[i]))-((0.106901773+23.9722/temp.vr[i])*sqrt(sal.vr[j]))+0.1130822*sal.vr[j]-0.00846934*(sal.vr[j]^1.5)+log(1-(0.001005*sal.vr[j])))
      logKsspcm_st[i,j] <- ((-171.9065-0.077993*temp.vr[i]+2839.319/temp.vr[i]+71.595*(log(temp.vr[i])/log(10))+(-0.77712+0.0028426*temp.vr[i]+178.34/temp.vr[i])*(sal.vr[j]^0.5)-0.07711*sal.vr[j]+0.0041249*(sal.vr[j]^1.5)))
      Ksspcm_st[i,j] <- 10^(logKsspcm_st[i,j])
      lnKsB_st[i,j] <- ((-8966.9-2890.53*sal.vr[j]^0.5-77.942*sal.vr[j]+1.728*sal.vr[j]^1.5-0.0996*sal.vr[j]^2)/temp.vr[i])+148.0248+137.1942*sal.vr[j]^0.5+1.62142*sal.vr[j]-(24.4344+25.085*sal.vr[j]^0.5+0.2474*sal.vr[j])*(log(temp.vr[i]))+(0.053105*sal.vr[j]^0.5*temp.vr[i])
      KsB_st[i,j] <- exp(lnKsB_st[i,j])
      Ksw_st[i,j] <- exp(148.96502-13847.26/temp.vr[i]-23.6521*(log(temp.vr[i]))+(118.67/temp.vr[i]-5.977+1.0495*(log(temp.vr[i])))*(sal.vr[j]^0.5)-0.01615*sal.vr[j])
      K0a[i,j] <- exp(9345.17/temp.vr[i]-60.2409+23.3585*(log(temp.vr[i]/100))+sal.vr[j]*(0.023517-0.00023656*temp.vr[i]+0.0047036*((temp.vr[i]/100)^2)))
      
      delV1[i] <- (-25.50)+0.1271*tempC.vr[i]
      delV2[i] <- (-15.82)+(-0.0219*tempC.vr[i])
      delVspc[i]<- (-48.76)+(0.5304*tempC.vr[i])
      delVB[i] <- (-29.48)+0.1622*tempC.vr[i]+(2.608/1000)*tempC.vr[i]^2
      delVw[i] <- (-25.60)+0.2324*tempC.vr[i]+(-3.6246/1000)*tempC.vr[i]^2
      
      delk1[i] <- (-3.08/1000)+(0.0877/1000)*tempC.vr[i]
      delk2[i] <- (1.13/1000)+(-0.1475/1000)*tempC.vr[i]
      delkspc[i] <- (-11.76/1000)+(0.3692/1000)*tempC.vr[i]
      delkB[i] <- -2.84/1000
      delkw[i] <- (-5.13/1000)+(0.0794/1000)*tempC.vr[i]
      
      Ks1a[i,j,k] <- (exp(-((delV1[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delk1[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ks1m_st[i,j]
      Ks2a[i,j,k] <- (exp(-((delV2[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delk2[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ks2m_st[i,j]
      Ksspca[i,j,k] <- (exp(-((delVspc[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delkspc[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ksspcm_st[i,j]
      KsBa[i,j,k] <- (exp(-((delVB[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delkB[i])/(R*temp.vr[i]))*press.vr[k]^2))*KsB_st[i,j]
      Kswa[i,j,k] <- (exp(-((delVw[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delkw[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ksw_st[i,j]
    }
  }
}


############################################################################################
#    INVERSION DRIVER - Read in and groom data; input environmental model prior distributions 
############################################################################################

# These parameters will be recorded in the output
parms = c("sal", "tempC", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "d18Osw.local",
          "pco2", "dic", "pH", "press", "m.vec", "c.vec")#"m.Grub", "m.Tsac", "m.Ouni", "m.custom", "c.Grub", "c.Tsac", "c.Ouni", "c.custom")

# Read in proxy time series data
prox.in <- read.csv('Harper.et.al./data/ShatskyLPEE.csv')
prox.in <- prox.in[,c(1:6)]
names(prox.in) <- c("age","d11B", "d11Bsd", "d18O", "MgCa", "species")

# Setup age range and bins 

ages.prox = unique(round(prox.in$age))
ages.prox.max = max(ages.prox)
dt = abs(diff(ages.prox, lag=1))
ages.prox.ai = seq(1,length(ages.prox), by=1)

# Age index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1))))

# parse clean.d11B by clean.d11B$species (clean.d11B1, clean.d11B2)
clean.d11B <- prox.in[complete.cases(prox.in$d11B), ]
clean.d11B <- transform(clean.d11B, si=ifelse(clean.d11B$species=="Grub",1,
                                              ifelse(clean.d11B$species=="Tsac",2,
                                                     ifelse(clean.d11B$species=="Ouni",3,
                                                                   ifelse(clean.d11B$species=="custom",4,5)))))
clean.d11B <- clean.d11B[-6]
clean.mgca <- prox.in[complete.cases(prox.in$MgCa), ]
mgcafu <- clean.mgca$MgCa*0.015
clean.d18O <- prox.in[complete.cases(prox.in$d18O), ]

# Vector of age indexes that contain d11B proxy data (with duplicates)
ai.d11B <- c(clean.d11B$ai)    

# Vector of age indexes that contain Mg/Ca proxy data
ai.mgca <- c(clean.mgca$ai)     

# Vector of age indexes that contain d18O proxy data
ai.d18O <- c(clean.d18O$ai)

ai.all <- c(ai.d11B, ai.mgca, ai.d18O)

# Index vector which contains each environmental time step that has one or more proxy data

ai.prox <-  unique(ai.all)     
ai.prox <- sort(ai.prox, decreasing = FALSE) 
n.steps = length(ai.prox)

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.d11B = match(ai.d11B, ai.prox)
ai.mgca = match(ai.mgca, ai.prox)
ai.d18O = match(ai.d18O, ai.prox)

############################################################################################
# Input prior mean and precision estimates for environmental parameters 

sal.m = 35  
sal.p = 1/0.5^2    

tempC.m = 30   
tempC.p = 1/5^2 

xca.m = 21.41 # 21.41 @ 59 Ma, 20.84 @ 56 Ma [Holland et al. (2020); i.e., -0.00019*dt]
xca.p = 1/0.5^2

xmg.m = 68.51 # ~37 throughout LPEE, OR 68.51 @59 and 37 @ 56 (i.e., -0.00274*dt); paired benthic Mg/Ca+d18O from Shatsky suggest large decrease in Mg/Casw during LPEE; max Cenozoic change in Mg is 2.74/Myr
xmg.p = 1/0.5^2

xso4.m = 14
xso4.p = 1/0.5^2

d11Bsw.m = 38.45
d11Bsw.p = 1/0.5^2

d18Osw.m = -1.2
d18Osw.p = 1/0.1^2

pH.l = 7.45
pH.u = 7.75


############################################################################################
# DIC prior(s) 

# Specify precision on DIC prior(s)
dic.p <- 1/0.00015^2

# Priors at each time step, or prior at time step 1 with autoregressive step function?


# Read in DIC time series data and linearly interpolate DIC for ages associated with each time step using input DIC time series 
dic.in <- read.csv('Harper.et.al./data/LOSCAR.DIC.csv')
dic.mod <- lm(dic.in$dic ~ dic.in$age, data = dic.in)
dic.interp <- approx(dic.in$age, dic.in$dic, xout=ages.prox, ties = mean, method = "linear")
dic.sim <- dic.interp$y #mol/kg


############################################################################################
# Data to pass to jags
data <- list("d11Bf.data" = clean.d11B$d11B,
             "d11Bfu.data" = clean.d11B$d11Bsd,
             "si.d11B" = clean.d11B$si,
             "d18Of.data" = clean.d18O$d18O, 
             "mgcaf.data" = clean.mgca$MgCa,
             "mgcafu.data" = mgcafu,
             "n.steps" = n.steps,
             "dt" = dt,
             "ages.prox" = ages.prox,
             "ai.prox" = ai.prox, 
             "ai.d11B" = ai.d11B,
             "ai.d18O" = ai.d18O, 
             "ai.mgca" = ai.mgca, 
             "m.mean" = m.mean,
             "m.prec" = m.prec,
             "c.mean" = c.mean,
             "c.prec" = c.prec,
             "seccal" = seccal, 
             "seccal.u" = seccal.u, 
             "d18Oseccal" = d18Oseccal, 
             "Hp.mean" = Hp.mean, 
             "Hp.sd" = Hp.sd, 
             "Bmod.mean" = Bmod.mean, 
             "Bmod.sd" = Bmod.sd, 
             "A.mean" = A.mean, 
             "A.sd" = A.sd,
             "pHpccorr" = pHpccorr,
             "pHpccorrsd" = pHpccorrsd,
             "sal.m" = sal.m,  
             "sal.p" = sal.p,   
             "tempC.m" = tempC.m, 
             "tempC.p" = tempC.p,
             "xca.m" = xca.m,
             "xca.p" = xca.p,
             "xmg.m" = xmg.m,
             "xmg.p" = xmg.p,
             "xso4.m" = xso4.m,
             "xso4.p" = xso4.p,
             "d11Bsw.m" = d11Bsw.m,
             "d11Bsw.p" = d11Bsw.p,
             "d18Osw.m" = d18Osw.m,
             "d18Osw.p" = d18Osw.p,
             "pH.l" = pH.l,  
             "pH.u" = pH.u, 
             "dic.sim" = dic.sim,
             "dic.p" = dic.p,
             "tempC.lb" = tempC.lb, 
             "sal.lb" = sal.lb,
             "press.lb" = press.lb,
             "press.ub" = press.ub,
             "t.inc" = t.inc,
             "s.inc" = s.inc,
             "p.inc" = p.inc,
             "K0a" = K0a,
             "Ks1a" = Ks1a,
             "Ks2a" = Ks2a,
             "Ksspca" = Ksspca,
             "KsBa" = KsBa,
             "Kswa" = Kswa)


############################################################################################
# Run the inversion

inv.out = R2jags::jags(model.file = "ForamPSM.lutKs.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = n.chains, n.iter = n.iter,
            n.burnin = n.burnin, n.thin = n.thin)
############################################################################################


