#
# Driver file for time series data inversion using forward foraminifera PSM 
# Reads in time series proxy data from rds objects
# Reads in Haynes & Hönisch (2020) DIC and LOSCAR-simulated DIC data, determines DIC priors at each time step
# Runs JAGS inversion using 'ForamPSMLPEE.R' model file
# Designed for late Paleocene early Eocene Shatsky Rise Data Inversion (Harper et al., submission to PNAS)



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
seccal = 70         # Percent recrystallized
seccal.u = 5        # 1sd % recrystallized
d18Oseccal = 0.85   # Calculated following Edgar et al. (2015)


############################################################################################  
# Mg/Ca SST proxy vital effects and calibration parameters 

# Nonlinearity of the relationship b/w shell and Mg/Casw (Haynes et al. subm., T. sacculifer) 
Hp.mean = 0.74      
Hp.sd = 0.05         

# Modern (pre-corrected) pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
Bmod.mean = 0.38    
Bmod.sd = 0.02      

# Exponential constant in Mg/Ca-SST calibration (Evans et al., 2016) for low (P-E) Mg/Casw
A.mean = 0.0757       
A.sd = 0.0045         


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
#    INVERSION DRIVER - Read in and groom data; input environmental model prior distributions 
############################################################################################

# These parameters will be recorded in the output
parms <- c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
           "pco2", "dic", "pH", "m.1", "m.2", "c.1", "c.2", "alpha")

# Read in proxy time series data
prox.in <- readRDS(file = "Harperetal_subm/data/ShatskyLPEE_data.rds")

# Setup age range and bins 
ages.prox <- unique(round(prox.in$age))
ages.prox.max <- max(ages.prox)
dt <- abs(diff(ages.prox, lag=1))
ages.prox.ai <- seq(1,length(ages.prox), by=1)

# Age index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1))))

# parse clean.d11B by clean.d11B$species (clean.d11B1, clean.d11B2)
clean.d11B <- prox.in[complete.cases(prox.in$d11B), ]
clean.d11Bs <- split(clean.d11B, clean.d11B$species)
clean.d11B1 <- clean.d11Bs$Grub
clean.d11B2 <- clean.d11Bs$Tsac

clean.mgca <- prox.in[complete.cases(prox.in$MgCa), ]
mgcafu <- clean.mgca$MgCa*0.015
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
n.steps <- length(ai.prox)

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.d11B1 <- match(ai.d11B1, ai.prox)
ai.d11B2 <- match(ai.d11B2, ai.prox)
ai.mgca <- match(ai.mgca, ai.prox)
ai.d18O <- match(ai.d18O, ai.prox)


############################################################################################
# Input prior mean and precision estimates for environmental parameters 

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

d18Osw.m = -1.2
d18Osw.p = 1/0.1^2

pH.l = 7.45
pH.u = 7.75


############################################################################################
# Read in DIC from LOSCAR simulation output: mean of 2 sims for PETM and 2 sims for ETM-2 with ZT19 long-term carb chem

dic_LOSCAR <- readRDS(file = "Harperetal_subm/data/dic_LOSCAR.rds")
# Linearly interpolate DIC for ages associated with each time step using input DIC time series 
dic.interp.LOSCAR <- approx(dic_LOSCAR$dic.LOSCAR.x, dic_LOSCAR$dic.LOSCAR.y, xout=ages.prox, method="linear") 
dic.interp.LOSCAR.err <- approx(dic_LOSCAR$dic.LOSCAR.x, dic_LOSCAR$dic.LOSCAR.2s, xout=ages.prox, method="linear") 
dic.LOSCAR <- dic.interp.LOSCAR[["y"]]
dic.LOSCAR.err <- dic.interp.LOSCAR.err[["y"]] / 2

# Read in DIC time series from Haynes and Hönisch (2020)
dic_HH <- readRDS("Harperetal_subm/data/dic_HH.rds")
# Linearly interpolate DIC for ages associated with each time step using input DIC time series 
dic.interp.HH <- approx(dic_HH$dic.HH.x, dic_HH$dic.HH.y, xout=ages.prox, method="linear") 
dic.HH.meanerr <- rowMeans(cbind(dic_HH$dic.HH.2sp, dic_HH$dic.HH.2sn))
dic.interp.HH.err <- approx(dic_HH$dic.HH.x, dic.HH.meanerr, xout=ages.prox, method="linear") 
dic.HH <- dic.interp.HH[["y"]]
dic.HH.err <- dic.interp.HH.err[["y"]] / 2 

# Compute inverse variance weighted average + variance of DIC from LOSCAR output and HH20 B/Ca-based reconstruction 
dic.avg.pri <- (dic.LOSCAR*(1/(dic.LOSCAR.err^2))+dic.HH*(1/(dic.HH.err^2))) / (1/dic.LOSCAR.err^2 + 1/dic.HH.err^2)
dic.var.pri <- 1 / (1/dic.LOSCAR.err^2 + 1/dic.HH.err^2)
dic.var.pri <- replace(dic.var.pri, dic.var.pri<0.0001^2, 0.0001^2)

# Truncation values for DIC prior normal distribution; min and max of HH20 + LOSCAR sim approaches, plus 
# added error representative of difference in PETM cGENIE (Gutjahr et al 2017) and LOSCAR (this study) DIC
minmaxerr <- 0.0003 # i.e., from LOSCAR and cGENIE PETM DIC differences 
dic.max <- pmax(dic.LOSCAR, dic.HH) + minmaxerr
dic.min <- pmin(dic.LOSCAR, dic.HH) - minmaxerr

# Plot the DIC prior with inverse variance weighted average, 95% CI and distribution truncation - Figure S3 in Harper et al., in prep. 
dic.pri <- data.frame((ages.prox/10^3), (dic.avg.pri*10^3), (2*sqrt(dic.var.pri)*10^3), (dic.min*10^3), (dic.max*10^3))
names(dic.pri) <- c("age", "wavg","twosd", "min","max")

ggplot() +
  geom_ribbon(data = dic.pri, aes(x=age, ymin=(wavg+twosd), ymax=wavg-twosd), fill = "gray") +
  geom_line(data = dic.pri, aes(x=age, y=wavg), color = "black") +
  geom_line(data = dic.pri, aes(x=age, y=min), color = "black", linetype=3) +
  geom_line(data = dic.pri, aes(x=age, y=max), color = "black", linetype=3) +
  scale_x_reverse() +
  labs(x = "Age (Ma)", y = expression("DIC (mmol/kg)")) +
  theme_bw() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"),
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"))

############################################################################################
# Data to pass to jags
data <- list("d11Bf.data1" = clean.d11B1$d11B, 
            "d11Bfu.data1" = clean.d11B1$d11Bse, 
            "d11Bf.data2" = clean.d11B2$d11B, 
            "d11Bfu.data2" = clean.d11B2$d11Bse, 
            "d18Of.data" = clean.d18O$d18O, 
            "mgcaf.data" = clean.mgca$MgCa,
            "mgcafu.data" = mgcafu,
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
            "dic.avg.pri" = dic.avg.pri,
            "dic.var.pri" = dic.var.pri,
            "dic.min" = dic.min,
            "dic.max" = dic.max)


############################################################################################
# Run the inversion

jout = jags.parallel(model.file = "Harperetal_subm/code/_ForamPSMLPEE.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 9, n.iter = 800000,
            n.burnin = 500000, n.thin = 100)
# 500k burn in, 9 chains, 800k iterations takes 10 hours


############################################################################################
# Display summary statistics and save summary as .csv

View(jout$BUGSoutput$summary)
write.csv(jout$BUGSoutput$summary, "Harperetal_subm/out/inversion_sum.csv")

############################################################################################


