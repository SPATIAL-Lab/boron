
# Driver file for single sample data inversion using forward foraminifera PSM


############################################################################################
# Load libraries 
library(rjags)
library(R2jags)
library(tidyverse)
############################################################################################

############################################################################################
#    INPUT TO PASS TO JAGS
############################################################################################    
# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 60         # Percent recrystallized
seccal.u = 5        # 1sd % recrystallized
d18Oseccal = 0.85   # Calculated following Edgar et al. (2015)

# sea surface pH sensitivity of d18O planktic foram (Dd18O/DpH; Spero et al., 1997 for symbiont-bearing O. universa = -0.89)
d18O_pHcorr.avg = 0
d18O_pHcorr.sd = 0

############################################################################################  
# Mg/Ca SST proxy vital effects and calibration parameters 

# Nonlinearity of the relationship b/w shell and Mg/Casw (Haynes et al. 2023, T. sacculifer) 
Hp.mean = 0.74      
Hp.sd = 0.05         

# Modern (pre-corrected) pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
Bmod.mean = 0.38    
Bmod.sd = 0.02      

# Exponential constant in Mg/Ca-SST calibration (all species regression of Gray and Evans 2019), with 
# expanded uncertainty (10 times the sd of their regression) for low (P-E) Mg/Casw
A.mean = 0.061  
A.sd = 0.005     

# sea surface pH sensitivty of Mg/Ca planktic foram (as a percent change per 0.1 pH unit)
mgca_pHcorr.avg = 0
mgca_pHcorr.sd = 0


############################################################################################  
# Modern d11Bforam-d11Bborate "vital effect" calibration

# Adjust the offset 'c' for paleo application
c.correction1 = -3.76   # Correction set to get c = 5.76 for M. vel (i.e., intercept when PETM values are plotted versus borate d11B [calc'd from A. sol as G. ruber])
c.correction2 = -2.8    # Correction set to align PETM d11Borate reconstruction for the two species

# G. ruber 
m.Grub = 0.62     # mean "m" value for G. ruber distribution 
m.Grubu = 0.055   # s.d. for "m" value for G. ruber distribution
c.Grub = 9.52     # mean "c" value for G. ruber distribution 
c.Grubu = 1.01    # s.d. for "c" value for G. ruber distribution

# T. sacculifer
m.Tsac = 0.73     # mean "m" value for T. sacculifer distribution 
m.Tsacu = 0.04    # s.d. for "m" value for T. sacculifer distribution
c.Tsac = 6.42     # mean "c" value for T. sacculifer distribution 
c.Tsacu = 0.82    # s.d. for "c" value for T. sacculifer distribution


############################################################################################
# Input environmental priors

sal.m = 35  
sal.p = 1/0.5^2    

tempC.m = 30   
tempC.p = 1/3^2 

xca.m = 21.41 # 21.41 @ 59 Ma, 20.84 @ 56 Ma [Holland et al. (2020); i.e., -0.00019*dt]
xca.p = 1/0.5^2

xmg.m = 68.51 # ~37 throughout LPEE, OR 68.51 @59 and 37 @ 56 (i.e., -0.00274*dt); paired benthic Mg/Ca+d18O from Shatsky suggest large decrease in Mg/Casw during LPEE; max Cenozoic change in Mg is 2.74/Myr
xmg.p = 1/0.5^2

xso4.m = 14
xso4.p = 1/0.5^2

d11Bsw.m = 38.45
d11Bsw.p = 1/0.5^2

d18Osw.m = -0.7
d18Osw.p = 1/0.1^2

pH.l = 7.45
pH.u = 7.75

dic.m = 0.002
dic.sd = 0.00001
  
DDco3.m = 50
DDco3.sd = 1


############################################################################################
# These parameters will be recorded in the output
parms <- c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "pco2", "dic", 
           "pH", "m.1", "m.2", "c.1", "c.2", "alpha", "d11Bf.1", "d11Bf.2", "d18Of", "mgcaf")


############################################################################################
# Input single sample dataset 
d11Bf.data1 = 15.4
d11Bfu.data1 = 0.2
d11Bf.data2 = 15.7
d11Bfu.data2 = 0.15
d18Of.data = -1
mgcaf.data = 4
mgcafu.data = 0.3
  

############################################################################################
# Data to pass to jags
data <- list("d11Bf.data1" = d11Bf.data1, 
             "d11Bfu.data1" = d11Bfu.data1, 
             "d11Bf.data2" = d11Bf.data2, 
             "d11Bfu.data2" = d11Bfu.data2, 
             "d18Of.data" = d18Of.data, 
             "mgcaf.data" = mgcaf.data,
             "mgcafu.data" = mgcafu.data,
             "m.Grub" = m.Grub,
             "m.Grubu" = m.Grubu,
             "c.Grub" = c.Grub,
             "c.Grubu" = c.Grubu,
             "m.Tsac" = m.Tsac,
             "m.Tsacu" = m.Tsacu,
             "c.Tsac" = c.Tsac,
             "c.Tsacu" = c.Tsacu,
             "seccal" = seccal, 
             "seccal.u" = seccal.u, 
             "d18Oseccal" = d18Oseccal, 
             "c.correction1" = c.correction1,
             "c.correction2" = c.correction2,
             "Hp.mean" = Hp.mean, 
             "Hp.sd" = Hp.sd, 
             "Bmod.mean" = Bmod.mean, 
             "Bmod.sd" = Bmod.sd, 
             "A.mean" = A.mean, 
             "A.sd" = A.sd,
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
             "dic.m" = dic.m,
             "dic.sd" = dic.sd,
             "DDco3.m" = DDco3.m,
             "DDco3.sd" = DDco3.sd,
             "mgca_pHcorr.avg" = mgca_pHcorr.avg,
             "mgca_pHcorr.sd" = mgca_pHcorr.sd,
             "d18O_pHcorr.avg" = d18O_pHcorr.avg,
             "d18O_pHcorr.sd" = d18O_pHcorr.sd)


############################################################################################
# Run the inversion
jout = jags.parallel(model.file = "Harperetal_subm/RevisionApril2024/code/singlesamplePSM.R", 
                     parameters.to.save = parms,
                     data = data, inits = NULL, n.chains = 9, n.iter = 800000,
                     n.burnin = 500000, n.thin = 3)

############################################################################################



