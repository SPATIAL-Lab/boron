#
# Driver for time series data inversion using forward MagOxyPSM 
# Reads in proxy time series data and prior means, or lower and upper bounds
# 
# Dustin T. Harper
# 27 April 2023


#########################################################################################
# Load libraries 
library(rjags)
library(R2jags)
############################################################################################

############################################################################################
#    INPUT TO PASS TO JAGS
############################################################################################    
# Input approximate % of calcite non-primary (recrystallized), with 1sd uncertainty,
# and estimated Dd18O of primary-inorganic (secondary) calcite
seccal = 28    # Percent recrystallized
seccal.u = 2.5     # 1sd % recrystallized
d18Oseccal = 0.85 # Calculated following Edgar et al. (2015)


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
#    INVERSION DRIVER - 1st data set 
############################################################################################
# These parameters will be recorded in the output
parms <- c("sal", "tempC", "xca", "xmg", "d18Osw", "d18Osw.sc", "indexop")

# Read in proxy time series data
prox.in <- read.csv('1258PETM/data/1258PETM.csv')
prox.in <- prox.in[,c(1:5)]
names(prox.in) <- c("depth", "age", "d13C", "d18O", "MgCa")

# Setup age range and bins 

ages.prox <- unique(round(prox.in$age))
ages.prox.max <- max(ages.prox)
dt <- abs(diff(ages.prox, lag=1))
ages.prox.ai <- seq(1,length(ages.prox), by=1)

# Age index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1))))

clean.mgca <- prox.in[complete.cases(prox.in$MgCa), ]
mgcafu <- clean.mgca$MgCa*0.015
clean.d18O <- prox.in[complete.cases(prox.in$d18O), ]

# Vector of age indexes that contain Mg/Ca proxy data
ai.mgca <- c(clean.mgca$ai)     

# Vector of age indexes that contain d18O proxy data
ai.d18O <- c(clean.d18O$ai)

ai.all <- c(ai.mgca, ai.d18O)

# Index vector which contains each environmental time step that has one or more proxy data

ai.prox <-  unique(ai.all)     
ai.prox <- sort(ai.prox, decreasing = FALSE) 
n.steps <- length(ai.prox)

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.mgca <- match(ai.mgca, ai.prox)
ai.d18O <- match(ai.d18O, ai.prox)


############################################################################################
# Input prior mean and precision estimates for environmental parameters 

sal.m = 35  
sal.p = 1/0.5^2    

tempC.m = 33  
tempC.p = 1/2^2 

xca.m = 21.0028 
xca.p = 1/0.5^2

xmg.m = 37
xmg.p = 1/0.5^2

d18Osw.m = -1.2
d18Osw.p = 1/0.1^2


############################################################################################
# Data to pass to jags
data <- list("d18Of.data" = clean.d18O$d18O, 
            "mgcaf.data" = clean.mgca$MgCa,
            "mgcafu.data" = mgcafu, 
            "n.steps" = n.steps,
            "dt" = dt,
            "ages.prox" = ages.prox,
            "ai.prox" = ai.prox, 
            "ai.d18O" = ai.d18O, 
            "ai.mgca" = ai.mgca, 
            "seccal" = seccal, 
            "seccal.u" = seccal.u, 
            "d18Oseccal" = d18Oseccal, 
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
            "d18Osw.m" = d18Osw.m,
            "d18Osw.p" = d18Osw.p)


############################################################################################
# Run the inversion
jout = jags.parallel(model.file = "1258PETM/MagOxyPSM.R", parameters.to.save = parms,
                     data = data, inits = NULL, n.chains = 9, n.iter = 10000,
                     n.burnin = 3000, n.thin = 10)
############################################################################################

write.csv(jout$BUGSoutput$summary, "1258PETM/inversion.sum.csv")

