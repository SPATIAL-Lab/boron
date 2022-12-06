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

# Input pH correction on Mg/Ca as the percent change per 0.1 pH unit 
# (e.g., Hollis et al., 2019 recommend 7% per 0.1 pH unit, plus 1sd unc = 0.9%)
pHpccorr = 2
pHpccorrsd = 0.9

pHpccorr = pHpccorr/10
pHpccorrsd = pHpccorrsd/10

# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 60     # Percent recrystallized; DATA FOR JAGS
seccal.u = 5      # 1sd % recrystallized; DATA FOR JAGS
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
c.correction1 = -3.76   # Correction set to get c = 5.76 for M. vel (i.e., intercept when PETM values are plotted versus borate d11B [calc'd from A. sol as G. ruber])
c.correction2 = -1.9   # correction set to align PETM d11Borate reconstruction for the two species
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
#    INVERSION DRIVER - 1st data set 
############################################################################################
# Read in proxy time series data
input.df = "data/Input_data_TS.xlsx"
prox.in = read.xlsx(input.df, sheet = "ShatskyLPEE")
prox.in = prox.in[,c(1:6)]
names(prox.in) = c("age","d11B", "d11Bsd", "d18O", "MgCa", "species")

# Setup age range and bins 

ages.prox = unique(round(prox.in$age))
ages.prox.max = max(ages.prox)
dt = abs(diff(ages.prox, lag=1))
ages.prox.ai = seq(1,length(ages.prox), by=1)

# Age index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1))))

# parse clean.d11B by clean.d11B$species (clean.d11B1, clean.d11B2)
clean.d11B <- prox.in[complete.cases(prox.in$d11B), ]
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
n.steps = length(ai.prox)

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.d11B1 = match(ai.d11B1, ai.prox)
ai.d11B2 = match(ai.d11B2, ai.prox)
ai.mgca = match(ai.mgca, ai.prox)
ai.d18O = match(ai.d18O, ai.prox)


# Input prior mean and precision estimates 

sal.m = 35  
sal.p = 1/0.5^2    

tempC.m = 30   
tempC.p = 1/5^2 

xca.m = 21.5842 # 21.5842 @ 59 Ma, 21.0028 @ 56 Ma [Holland et al. (2020); i.e., -0.0001938*dt]
xca.p = 1/0.5^2

xmg.m = 45.21 # ~37 throughout LPEE, OR 45.21 @59 and 37 @ 56 (i.e., -0.00274*dt); paired benthic Mg/Ca+d18O from Shatsky suggest large decrease in Mg/Casw during LPEE; max Cenozoic change in Mg is 2.74/Myr
xmg.p = 1/0.5^2

xso4.m = 14
xso4.p = 1/0.5^2

d11Bsw.m = 38.45
d11Bsw.p = 1/0.5^2

d18Osw.m = -1
d18Osw.p = 1/0.1^2

pH.l = 7.5  
pH.u = 7.9

# DIC priors read in from LOSCAR output 
dic.p <- 1/0.0001^2
# Read in DIC time series data
input.dic = "data/Input_data_dic.xlsx"
dic.in = read.xlsx(input.dic, sheet = "LPEE")
dic.in.x <- dic.in[,1]
dic.in.y <- dic.in[,2]     #mol/kg
dic.in.df <- data.frame(dic.in.x, dic.in.y)
# Linearly interpolate DIC for ages associated with each time step using input DIC time series 
dic.mod <- lm(dic.in.y ~ dic.in.x, data = dic.in.df)
dic.interp <- approx(dic.in.df$dic.in.x, dic.in.df$dic.in.y, xout=ages.prox, method="linear") 
dic.sim <- dic.interp[["y"]]


# Data to pass to jags
data = list("d11Bf.data1" = clean.d11B1$d11B, 
            "d11Bfu.data1" = clean.d11B1$d11Bsd, 
            "d11Bf.data2" = clean.d11B2$d11B, 
            "d11Bfu.data2" = clean.d11B2$d11Bsd, 
            "d18Of.data" = clean.d18O$d18O, 
            "mgcaf.data" = clean.mgca$MgCa,
            "n.steps" = n.steps,
            "dt" = dt,
            "ages.prox" = ages.prox,
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
            "dic.p" = dic.p)

# Run the inversion

jout = jags.parallel(model.file = "boronPSM_TS_pHdicv6.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 3, n.iter = 2000,
            n.burnin = 1000, n.thin = 1)


