#
# Driver file for time series data inversion using forward foraminifera PSM 
# Reads in time series proxy data from rds objects
# Reads in Haynes & Hönisch (2020) DIC and LOSCAR-simulated DIC data, determines DIC priors at each time step
# Runs JAGS inversion using 'ForamPSMLPEE.R' model file
# Designed for late Paleocene early Eocene Shatsky Rise Data Inversion (Harper et al., submission to PNAS)



############################################################################################
# Load libraries 
#library(rjags)
library(R2jags)

############################################################################################

############################################################################################
#    INPUT TO PASS TO JAGS
############################################################################################    
# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 60         # Percent recrystallized
seccal.u = 5        # 1sd % recrystallized
d18Oseccal = 0.85   # Calculated following Edgar et al. (2015)

# sea surface pH sensitivity of d18O planktic foram (Dd18O/DpH; Spero et al., 1997 for symbiont-bearing O. universa = -0.89, for pH > 8.0)
d18O_pHcorr.avg = -0.89
d18O_pHcorr.sd = 0.2


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
#    INVERSION DRIVER - Read in and groom data; input environmental model prior distributions 
############################################################################################

# These parameters will be recorded in the output
parms <- c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "pco2", "alk", "pH", 
           "m.1", "m.2", "c.1", "c.2", "alpha", "d11Bf.1", "d11Bf.2", "d18Of", "mgcaf")

# Read in proxy time series data
prox.in <- readRDS(file = "Harperetal_resubm/data/ShatskyLPEE_data.rds")

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

d18Osw.m = -0.7
d18Osw.p = 1/0.1^2

pH.l = 7.45
pH.u = 7.75


############################################################################################
# Read in D[CO3=] from file and linearly interpolate for ages associated with each time step 

Dco3_in <- as.data.frame(readRDS(file = "Harperetal_resubm/data/Dco3_bw.rds"))
Dco3.interp <- approx(Dco3_in$age, Dco3_in$Dco3, xout=ages.prox, method="linear") 
Dco3_2s.interp <- approx(Dco3_in$age, Dco3_in$Dco3_2s, xout=ages.prox, method="linear") 
Dco3.avg.pri <- Dco3.interp[["y"]]
Dco3.sd.pri <- Dco3_2s.interp[["y"]] / 2
DDco3 <- vector("numeric")
DDco3.sd <- vector("numeric")
for (i in 1:length(Dco3.avg.pri)){
  if (Dco3.avg.pri[i] < 21.3){
    DDco3[i] <- Dco3.avg.pri[i]-21.3
    DDco3.sd[i] <- Dco3.sd.pri[i]
  } else {
    DDco3[i] <- 0
    DDco3.sd[i] <- 1e-20
    }
}

############################################################################################
# # Read in alk from LOSCAR simulation output: 5800 PETM and 3800 ETM-2 with 2.2 mmol/kg LPEE background
# library(readxl)
# alk_PETM_hi <- read_xlsx("Harperetal_subm/RevisionApril2024/data/alk_LOSCAR.xlsx", sheet = "alk_PETM_hi")
# alk_PETM_low <-read_xlsx("Harperetal_subm/RevisionApril2024/data/alk_LOSCAR.xlsx", sheet = "alk_PETM_low")
# alk_ETM2_hi <- read_xlsx("Harperetal_subm/RevisionApril2024/data/alk_LOSCAR.xlsx", sheet = "alk_ETM2_hi")
# alk_ETM2_low <- read_xlsx("Harperetal_subm/RevisionApril2024/data/alk_LOSCAR.xlsx", sheet = "alk_ETM2_low")
# 
# alk_PETM_hi_interp <- approx(alk_PETM_hi$age, alk_PETM_hi$alk, xout=alk_PETM_low$age, method="linear") 
# alk_PETM <- cbind(alk_PETM_low$age, (alk_PETM_hi_interp$y + alk_PETM_low$alk)/2)
# colnames(alk_PETM) <- c("age", "alk")
# 
# alk_ETM2_hi_interp <- approx(alk_ETM2_hi$age, alk_ETM2_hi$alk, xout=alk_ETM2_low$age, method="linear") 
# alk_ETM2 <- cbind(alk_ETM2_low$age, (alk_ETM2_hi_interp$y + alk_ETM2_low$alk)/2)
# colnames(alk_ETM2) <- c("age", "alk")
# 
# alk_LOSCAR <- data.frame(rbind(alk_PETM, alk_ETM2))
# alk.interp.LOSCAR <- approx(alk_LOSCAR$age, alk_LOSCAR$alk, xout=ages.prox, method="linear")
# save(alk.interp.LOSCAR, file = "Harperetal_subm/RevisionApril2024/data/alk.interp.LOSCAR.rda")
load(file= "Harperetal_resubm/data/alk.interp.LOSCAR.rda")

alk.avg.pri <- (alk.interp.LOSCAR[["y"]]) / 1e3
alk.sd.pri <- rep(0.0001, times= length(alk.avg.pri))

alk.min <- alk.avg.pri - 0.0003
alk.max <- alk.avg.pri + 0.0003

# Plot the sea surface alkalinity prior mean 
plot(x = ages.prox, y = alk.avg.pri)


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
            "alk.avg.pri" = alk.avg.pri,
            "alk.sd.pri" = alk.sd.pri,
            "alk.min" = alk.min,
            "alk.max" = alk.max,
            "DDco3" = DDco3,
            "DDco3.sd" = DDco3.sd,
            "mgca_pHcorr.avg" = mgca_pHcorr.avg,
            "mgca_pHcorr.sd" = mgca_pHcorr.sd,
            "d18O_pHcorr.avg" = d18O_pHcorr.avg,
            "d18O_pHcorr.sd" = d18O_pHcorr.sd)


############################################################################################
# Run the inversion

system.time({jout = jags.parallel(model.file = "Harperetal_resubm/code/_ForamPSMLPEE_CCvar.R", 
                     parameters.to.save = parms, data = data, inits = NULL, 
                     n.chains = 9, n.iter = 800000, n.burnin = 500000, n.thin = 200)})
# 500k burn in, 9 chains, 800k iterations takes 10 hours


############################################################################################
# Display summary statistics and save summary as .csv

#View(jout$BUGSoutput$summary)
#write.csv(jout$BUGSoutput$summary, "Harperetal_resubm/out_senstest/inversion_sum.csv")

############################################################################################
save(jout, file = "Harperetal_resubm/out_senstest/LPEE_CCvar.rda")
