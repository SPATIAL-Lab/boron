#
# Driver for time series data inversion using forward Boron PSM 
# Reads in proxy time series data and prior means, or lower and upper bounds, as time series data
# 
# Dustin T. Harper
# 10 August 2022


# Load libraries 
library(rjags)
library(R2jags)
library(openxlsx)
library(tidyverse)


############################################################################################
#    INPUT DATA
############################################################################################    

# Setup age range and bins 
ages.bin = 0.01   # Myr; DATA FOR JAGS
ages.max = 59     # Myr; DATA FOR JAGS
ages.min = 53     # Myr; DATA FOR JAGS

# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 65        # Percent recrystallized; DATA FOR JAGS
seccal.u = 10      # 1sd % recrystallized; DATA FOR JAGS
Dd18Oseccal = 3.85 # Calculated following Edgar et al. (2015); DATA FOR JAGS


#---------------------------------------------------------------------------------------------------------
# Pressure prior (ie., depth habitat) in bar for each of the modern species' calibration options 

# G. ruber calibrated species' associated pressures at habitat depth 
Grub.press.m = 10        
Grub.press.sd = 2     

# T. sacculifer calibrated species' associated pressures at habitat depth 
Tsac.press.m = 10        
Tsac.press.sd = 2     

# O. universa calibrated species' associated pressures at habitat depth 
Ouni.press.m = 10        
Ouni.press.sd = 2   

# Borate calibrated species' associated pressures at habitat depth 
bor.press.m = 10        
bor.press.sd = 2     


#---------------------------------------------------------------------------------------------------------
# Mg/Ca SST proxy vital effects and calibration parameters 

# Nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer) 
Hp.mean = 0.41      # DATA FOR JAGS
Hp.sd = 0.1         # DATA FOR JAGS

# Modern (pre-corrected) pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
Bmod.mean = 0.38    # DATA FOR JAGS
Bmod.sd = 0.02      # DATA FOR JAGS

# Exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evanset al., 2016)
A.mean = 0.07       # DATA FOR JAGS
A.sd = 0.01         # DATA FOR JAGS


#---------------------------------------------------------------------------------------------------------
# Modern d11Bforam-d11Bborate "vital effect" calibration data 

# You can adjust the offset for the final value of 'c' here if you wish to see the effect
c.correction = 0   # DATA FOR JAGS

# G. ruber 
bor.Grub <- c(14.22, 16.66, 19.76)    # Henehan et al. (2013) G. ruber data (borate)
for.Grub <- c(18.2, 19.63, 21.46)     # Henehan et al. (2013) G. ruber data (calcite)
m.Grub = 0.6   # mean "m" value for G. ruber distribution 
m.Grubu = 0.055   # s.d. for "c" value for G. ruber distribution
c.Grub = 9.52  # mean "c" value for G. ruber distribution 
c.Grubu = 1.01 # s.d. for "c" value for G. ruber distribution

# T. sacculifer
bor.Tsac <- c(19.03, 19.15, 18.8, 19.57, 19.57, 19.57, 18.98, 18.32, 23.86, 14, 18.4)    # T. sacculifer data from Foster (2008), Martinez-Botti et al. core top (2015) [including Sanyal et al (2001) refit data] (borate)
for.Tsac <- c(19.95, 19.82, 19.43, 20.28, 20.13, 20.02, 19.60, 19, 23.66, 15.46, 18.49)     # T. sacculifer data from Foster (2008), Martinez-Botti et al. core top (2015) [including Sanyal et al (2001) refit data] (calcite)
m.Tsac = 0.82   # mean "m" value for T. sacculifer distribution 
m.Tsacu = 0.11   # s.d. for "c" value for T. sacculifer distribution
c.Tsac = 3.94  # mean "c" value for T. sacculifer distribution 
c.Tsacu = 2.01 # s.d. for "c" value for T. sacculifer distribution

# O. universa
bor.Ouni <- c(18.14, 18.34, 16.93, 17.51, 16.28, 16.61, 17.74, 16.93, 16.28, 16.93, 18.59, 16.28, 18.21, 19.57, 18.21, 19.82, 19.82, 19.76, 18.78)  # Henehan et al. (2016) O. universa data (borate)
for.Ouni <- c(16.17, 16.97, 15.98, 15.99, 15.15, 15.65, 16.29, 15.84, 15.19, 16.08, 17.5, 14.98, 16.62, 18.04, 16.82, 18.57, 18.98, 18.82, 19.47)   # Henehan et al. (2016) O. universa data (calcite)
m.Ouni = 0.95   # mean "m" value for O. universa distribution 
m.Ouniu = 0.085   # s.d. for "c" value for O. universa distribution
c.Ouni = -0.42  # mean "c" value for O. universa distribution 
c.Ouniu = 1.43 # s.d. for "c" value for O. universa distribution

# borate
bor.borate <- c(0, 1) # 1:1 line (borate) 
for.borate <- c(0, 1)  # 1:1 line (calcite)
m.borate = 1   # mean "m" value for borate distribution 
m.borateu = 0.1   # s.d. for "c" value for borate distribution
c.borate = 0  # mean "c" value for borate distribution 
c.borateu = 1 # s.d. for "c" value for borate distribution


#---------------------------------------------------------------------------------------------------------
# These parameters will be recorded in the output
parms = c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
          "pco2", "dic", "pH", "m", "c")


############################################################################################
#    INVERSION DRIVER
############################################################################################    

# Read in proxy time series data
input.df = "data/Input_data_TS.xlsx"
prox.in = read.xlsx(input.df, sheet = "proxy4")
prox.in = prox.in[,c(1:6)]
names(prox.in) = c("age","d11B", "d11Bsd", "d18O", "MgCa", "species")

# Age index proxy data
ages = seq(ages.max, ages.min, by = 0 - ages.bin) - ages.bin / 2
n.steps <- (ages.max - ages.min) / ages.bin
ages.len = length(ages)
prox.in$ai = ceiling((ages.max - prox.in$age) / ages.bin) 

clean.d11B <- prox.in[complete.cases(prox.in$d11B), ]
clean.mgca <- prox.in[complete.cases(prox.in$MgCa), ]
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

# Age index vector for prior time bins
ai.env = ceiling((ages.max - ages) / ages.bin)  

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.d11B = match(ai.d11B, ai.prox)
ai.mgca = match(ai.mgca, ai.prox)
ai.d18O = match(ai.d18O, ai.prox)

# Data to pass to jags
data = list("d11Bf.data" = clean.d11B$d11B, 
            "d11Bfu.data" = clean.d11B$d11Bsd, 
            "d18Of.data" = clean.d18O$d18O, 
            "mgcaf.data" = clean.mgca$MgCa,
            "spec.in" = prox.in$species,
            "n.steps" = n.steps, 
            "ai.prox" = ai.prox, 
            "ai.d11B" = ai.d11B, 
            "ai.d18O" = ai.d18O, 
            "ai.mgca" = ai.mgca, 
            "bor.Grub" = bor.Grub,
            "for.Grub" = for.Grub,
            "m.Grub" = m.Grub,
            "m.Grubu" = m.Grubu,
            "c.Grub" = c.Grub,
            "c.Grubu" = c.Grubu,
            "Grub.press.m" = Grub.press.m, 
            "Grub.press.sd" = Grub.press.sd,
            "bor.Tsac" = bor.Tsac,
            "for.Tsac" = for.Tsac,
            "m.Tsac" = m.Tsac,
            "m.Tsacu" = m.Tsacu,
            "c.Tsac" = c.Tsac,
            "c.Tsacu" = c.Tsacu,
            "Tsac.press.m" = Tsac.press.m, 
            "Tsac.press.sd" = Tsac.press.sd,
            "bor.Ouni" = bor.Ouni,
            "for.Ouni" = for.Ouni,
            "m.Ouni" = m.Ouni,
            "m.Ouniu" = m.Ouniu,
            "c.Ouni" = c.Ouni,
            "c.Ouniu" = c.Ouniu,
            "Ouni.press.m" = Ouni.press.m, 
            "Ouni.press.sd" = Ouni.press.sd,
            "bor.borate" = bor.borate,
            "for.borate" = for.borate,
            "m.borate" = m.borate,
            "m.borateu" = m.borateu,
            "c.borate" = c.borate,
            "c.borateu" = c.borateu,
            "borate.press.m" = bor.press.m,
            "borate.press.sd" = bor.press.sd,
            "seccal" = seccal, 
            "seccal.u" = seccal.u, 
            "Dd18Oseccal" = Dd18Oseccal, 
            "c.correction" = c.correction, 
            "Hp.mean" = Hp.mean, 
            "Hp.sd" = Hp.sd, 
            "Bmod.mean" = Bmod.mean, 
            "Bmod.sd" = Bmod.sd, 
            "A.mean" = A.mean, 
            "A.sd" = A.sd)

# Run the inversion

jout = jags(model.file = "boronPSM_TS_inv_pH.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 3, n.iter = 100,
            n.burnin = 1, n.thin = 10)

 # jout = jags(model.file = "boronPSM_TS_inv_pH.R", parameters.to.save = parms,
 #             data = data, inits = NULL, n.chains = 3, n.iter = 1e4,
 #             n.burnin = 1e3, n.thin = 10)

