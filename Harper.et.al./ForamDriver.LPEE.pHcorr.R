#
# Driver for time series data inversion using forward foraminifera PSM 
# Reads in time series proxy data and LOSCAR-simulated DIC data for DIC priors at each time step
# Runs JAGS inversion using 'Foram.PSM.LPEE.pHcorr.R'
# Designed for late Paleocene early Eocene Shatsky Rise Data Inversion (Harper et al., in prep.)
# 
# Dustin T. Harper
# 27 April 2023


############################################################################################
# Load libraries 
library(rjags)
library(R2jags)
############################################################################################

############################################################################################
#    INPUT TO PASS TO JAGS
############################################################################################    

# Input pH correction on Mg/Ca as the percent change per 0.1 pH unit 
  # Hollis et al., 2019 recommend 7% per 0.1 pH unit, plus 1sd unc = 0.9%; 
  # Haynes et al., submitted show Mg/Casw dependency, ~0% to 2% per 0.1 pH unit for Paleocene - Eocene Mg/Casw
pHpccorr = 2
pHpccorrsd = 0.9

pHpccorr = pHpccorr/10
pHpccorrsd = pHpccorrsd/10

# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 60         # Percent recrystallized
seccal.u = 2.5      # 1sd % recrystallized
d18Oseccal = 0.85   # Calculated following Edgar et al. (2015)


############################################################################################  
# Mg/Ca SST proxy vital effects and calibration parameters 

# Nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer) 
Hp.mean = 0.41      
Hp.sd = 0.1         

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
c.correction2 = -1.9   # correction set to align PETM d11Borate reconstruction for the two species

# G. ruber 
m.Grub = 0.62   # mean "m" value for G. ruber distribution 
m.Grubu = 0.055   # s.d. for "m" value for G. ruber distribution
c.Grub = 9.52  # mean "c" value for G. ruber distribution 
c.Grubu = 1.01 # s.d. for "c" value for G. ruber distribution

# T. sacculifer
m.Tsac = 0.82   # mean "m" value for T. sacculifer distribution 
m.Tsacu = 0.11   # s.d. for "m" value for T. sacculifer distribution
c.Tsac = 3.94  # mean "c" value for T. sacculifer distribution 
c.Tsacu = 2.01 # s.d. for "c" value for T. sacculifer distribution


############################################################################################
#    INVERSION DRIVER - Read in and groom data; input environmental model prior distributions 
############################################################################################

# These parameters will be recorded in the output
parms <- c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
           "pco2", "dic", "pH", "m.1", "m.2", "c.1", "c.2", "pH.phi")

# Read in proxy time series data
prox.in <- read.csv('Harper.et.al./data/ShatskyLPEE.csv')
prox.in <- prox.in[,c(1:6)]
names(prox.in) <- c("age","d11B", "d11Bsd", "d18O", "MgCa", "species")

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
# DIC priors read in from LOSCAR output 
dic.p <- 1/0.00015^2

# Read in DIC time series data
dic.in <- read.csv('Harper.et.al./data/LOSCAR.DIC.csv')
dic.in.x <- dic.in[,1]
dic.in.y <- dic.in[,2]     #mol/kg
dic.in.df <- data.frame(dic.in.x, dic.in.y)

# Linearly interpolate DIC for ages associated with each time step using input DIC time series 
dic.mod <- lm(dic.in.y ~ dic.in.x, data = dic.in.df)
dic.interp <- approx(dic.in.df$dic.in.x, dic.in.df$dic.in.y, xout=ages.prox, method="linear") 
dic.sim <- dic.interp[["y"]]


############################################################################################
# Data to pass to jags
data <- list("d11Bf.data1" = clean.d11B1$d11B, 
            "d11Bfu.data1" = clean.d11B1$d11Bsd, 
            "d11Bf.data2" = clean.d11B2$d11B, 
            "d11Bfu.data2" = clean.d11B2$d11Bsd, 
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


############################################################################################
# Run the inversion

jout = jags.parallel(model.file = "Harper.et.al./ForamPSM.LPEE.pHcorr.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 9, n.iter = 800000,
            n.burnin = 500000, n.thin = 400)
# 500k burn in, 9 chains, 800k iterations takes 10 hours
############################################################################################

write.csv(jout$BUGSoutput$summary, "Harper.et.al./inversion.sum.csv")


