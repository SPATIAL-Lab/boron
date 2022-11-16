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
dt = 5   # kyr; DATA FOR JAGS
ages.max = 58780     # kyr; DATA FOR JAGS
ages.min = 53080     # kyr; DATA FOR JAGS

# Input pH correction on Mg/Ca as the percent change per 0.1 pH unit 
# (e.g., Hollis et al., 2019 recommend 7% per 0.1 pH unit, plus 1sd unc = 0.9%)
pHpccorr = 0 # 7
pHpccorrsd = 10^-30 # 0.9

pHpccorr = pHpccorr/10
pHpccorrsd = pHpccorrsd/10

# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 65        # Percent recrystallized; DATA FOR JAGS
seccal.u = 10      # 1sd % recrystallized; DATA FOR JAGS
Dd18Oseccal = 3.85 # Calculated following Edgar et al. (2015); DATA FOR JAGS


#---------------------------------------------------------------------------------------------------------
# # Pressure prior (ie., depth habitat) in bar for each of the modern species' calibration options 
# 
# # G. ruber calibrated species' associated pressures at habitat depth 
# Grub.press.m =  5        
# Grub.press.sd = 1     
# 
# # T. sacculifer calibrated species' associated pressures at habitat depth 
# Tsac.press.m = 7        
# Tsac.press.sd = 1     
# 
# # O. universa calibrated species' associated pressures at habitat depth 
# Ouni.press.m = 10        
# Ouni.press.sd = 4   
# 
# # Borate calibrated species' associated pressures at habitat depth 
# bor.press.m = 6       
# bor.press.sd = 2     


#---------------------------------------------------------------------------------------------------------
# Mg/Ca SST proxy vital effects and calibration parameters 

# Nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer) 
Hp.mean = 0.41      # DATA FOR JAGS
Hp.sd = 0.1         # DATA FOR JAGS

# Modern (pre-corrected) pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
Bmod.mean = 0.38    # DATA FOR JAGS
Bmod.sd = 0.02      # DATA FOR JAGS

# Exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evanset al., 2016)
A.mean = 0.09       # DATA FOR JAGS
A.sd = 0.01         # DATA FOR JAGS


#---------------------------------------------------------------------------------------------------------
# Modern d11Bforam-d11Bborate "vital effect" calibration data 

# You can adjust the offset for the final value of 'c' here if you wish to see the effect
c.correction1 = -3.3   # DATA FOR JAGS
c.correction2 = -1.9
# G. ruber 
# bor.Grub <- c(14.22, 16.66, 19.76)    # Henehan et al. (2013) G. ruber data (borate)
# for.Grub <- c(18.2, 19.63, 21.46)     # Henehan et al. (2013) G. ruber data (calcite)
m.Grub = 0.62   # mean "m" value for G. ruber distribution 
m.Grubu = 0.0055   # s.d. for "c" value for G. ruber distribution
c.Grub = 9.52  # mean "c" value for G. ruber distribution 
c.Grubu = 1.01 # s.d. for "c" value for G. ruber distribution

# T. sacculifer
# bor.Tsac <- c(19.03, 19.15, 18.8, 19.57, 19.57, 19.57, 18.98, 18.32, 23.86, 14, 18.4)    # T. sacculifer data from Foster (2008), Martinez-Botti et al. core top (2015) [including Sanyal et al (2001) refit data] (borate)
# for.Tsac <- c(19.95, 19.82, 19.43, 20.28, 20.13, 20.02, 19.60, 19, 23.66, 15.46, 18.49)     # T. sacculifer data from Foster (2008), Martinez-Botti et al. core top (2015) [including Sanyal et al (2001) refit data] (calcite)
m.Tsac = 0.82   # mean "m" value for T. sacculifer distribution 
m.Tsacu = 0.011   # s.d. for "c" value for T. sacculifer distribution
c.Tsac = 3.94  # mean "c" value for T. sacculifer distribution 
c.Tsacu = 2.01 # s.d. for "c" value for T. sacculifer distribution

# # O. universa
# bor.Ouni <- c(18.14, 18.34, 16.93, 17.51, 16.28, 16.61, 17.74, 16.93, 16.28, 16.93, 18.59, 16.28, 18.21, 19.57, 18.21, 19.82, 19.82, 19.76, 18.78)  # Henehan et al. (2016) O. universa data (borate)
# for.Ouni <- c(16.17, 16.97, 15.98, 15.99, 15.15, 15.65, 16.29, 15.84, 15.19, 16.08, 17.5, 14.98, 16.62, 18.04, 16.82, 18.57, 18.98, 18.82, 19.47)   # Henehan et al. (2016) O. universa data (calcite)
# m.Ouni = 0.95   # mean "m" value for O. universa distribution 
# m.Ouniu = 0.085   # s.d. for "c" value for O. universa distribution
# c.Ouni = -0.42  # mean "c" value for O. universa distribution 
# c.Ouniu = 1.43 # s.d. for "c" value for O. universa distribution
# 
# # borate
# bor.borate <- c(0, 1) # 1:1 line (borate) 
# for.borate <- c(0, 1)  # 1:1 line (calcite)
# m.borate = 1   # mean "m" value for borate distribution 
# m.borateu = 0.1   # s.d. for "c" value for borate distribution
# c.borate = 0  # mean "c" value for borate distribution 
# c.borateu = 1 # s.d. for "c" value for borate distribution


#---------------------------------------------------------------------------------------------------------
# These parameters will be recorded in the output
parms = c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
          "pco2", "dic", "pH", "m.1", "m.2", "c.1", "c.2")


############################################################################################
#    INVERSION DRIVER
############################################################################################    

# Read in proxy time series data
input.df = "data/Input_data_TS.xlsx"
prox.in = read.xlsx(input.df, sheet = "ShatskyLPEE")
prox.in = prox.in[,c(1:6)]
names(prox.in) = c("age","d11B", "d11Bsd", "d18O", "MgCa", "species")

# Age index proxy data
ages = head(seq(ages.max, ages.min, by = 0 - dt) - dt / 2, -1)
n.steps <- (ages.max - ages.min) / dt
n.steps <- round(n.steps, digits=0)
ages.len = length(ages)
prox.in$ai = ceiling((ages.max - prox.in$age) / dt) 

clean.d11B <- prox.in[complete.cases(prox.in$d11B), ]
# parse clean.d11B by clean.d11B$species (clean.d11B1, clean.d11B2)
clean.d11Bs <- split(clean.d11B, clean.d11B$species)
clean.d11B1 <- clean.d11Bs$Grub
clean.d11B2 <- clean.d11Bs$Tsac

clean.mgca <- prox.in[complete.cases(prox.in$MgCa), ]
clean.d18O <- prox.in[complete.cases(prox.in$d18O), ]

# Vector of age indexes that contain d11B proxy data (with duplicates)
ai.d11B1 <- c(clean.d11B1$ai)    
ai.d11B2 <- c(clean.d11B2$ai)    

# Vector of age indexes that contain Mg/Ca proxy data
ai.mgca <- c(clean.mgca$ai)     

# Vector of age indexes that contain d18O proxy data
ai.d18O <- c(clean.d18O$ai)     

ai.all <- c(ai.d11B1, ai.d11B2, ai.mgca, ai.d18O)

# Index vector which contains each environmental time step that has one or more proxy data
ai.prox <-  unique(ai.all)     
ai.prox <- sort(ai.prox, decreasing = FALSE) 

# Age index vector for prior time bins
ai.env = ceiling((ages.max - ages) / dt)  

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.d11B1 = match(ai.d11B1, ai.prox)
ai.d11B2 = match(ai.d11B2, ai.prox)
ai.mgca = match(ai.mgca, ai.prox)
ai.d18O = match(ai.d18O, ai.prox)

# Data to pass to jags
data = list("d11Bf.data1" = clean.d11B1$d11B, 
            "d11Bfu.data1" = clean.d11B1$d11Bsd, 
            "d11Bf.data2" = clean.d11B2$d11B, 
            "d11Bfu.data2" = clean.d11B2$d11Bsd, 
            "d18Of.data" = clean.d18O$d18O, 
            "mgcaf.data" = clean.mgca$MgCa,
            "n.steps" = n.steps,
            "dt" = dt,
            "ai.prox" = ai.prox, 
            "ai.d11B1" = ai.d11B1,
            "ai.d11B2" = ai.d11B2,
            "ai.d18O" = ai.d18O, 
            "ai.mgca" = ai.mgca, 
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
            "Dd18Oseccal" = Dd18Oseccal, 
            "c.correction1" = c.correction1,
            "c.correction2" = c.correction2,
            "Hp.mean" = Hp.mean, 
            "Hp.sd" = Hp.sd, 
            "Bmod.mean" = Bmod.mean, 
            "Bmod.sd" = Bmod.sd, 
            "A.mean" = A.mean, 
            "A.sd" = A.sd,
            "pHpccorr" = pHpccorr,
            "pHpccorrsd" = pHpccorrsd)

# Run the inversion

jout = jags.parallel(model.file = "boronPSM_TS_pHdic.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 3, n.iter = 1e3,
            n.burnin = 500, n.thin = 1)

