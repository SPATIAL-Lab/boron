#
# Driver for Boron PSM inversion (model file = boronPSM_inv.R) using multiple, complete 
# foram d11B, d18O, and Mg/Ca datums from different species read from from a csv
# 18 July 2022
# Dustin Harper
#

# Load libraries and make empty lists 
library(R2jags)
library(tidyverse)
library(openxlsx)
multi.out <- list()
parms.out <- list()

# Input approximate % (seccal) of sample foram calcite non-primary (recrystalized), with 1sd uncertainty (ssecal.u),
# and estimated Dd18O of primary vs. inorganic (secondary) calcite
seccal <- 65        # Percent recrystallized (i.e., diagenetic effect)
seccal.u <- 10      # 1sd uncertainty in % seccal
Dd18Oseccal <- 3.85 # calculated following Edgar et al. (2015)

# you can adjust the offset for the final value of 'c' here if you wish to see the effect
c.correction <- 0


bor.Grub <- c(14.22, 16.66, 19.76)    # Henehan et al. (2013) G. ruber data (borate)
for.Grub <- c(18.2, 19.63, 21.46)     # Henehan et al. (2013) G. ruber data (calcite)
m.Grub <- 0.6   # mean "m" value for G. ruber distribution 
m.Grubu <- 0.055   # s.d. for "c" value for G. ruber distribution
c.Grub <- 9.52  # mean "c" value for G. ruber distribution 
c.Grubu <- 1.01 # s.d. for "c" value for G. ruber distribution

bor.Tsac <- c(19.03, 19.15, 18.8, 19.57, 19.57, 19.57, 18.98, 18.32, 23.86, 14, 18.4)    # T. sacculifer data from Foster (2008), Martinez-Botti et al. core top (2015) [including Sanyal et al (2001) refit data] (borate)
for.Tsac <- c(19.95, 19.82, 19.43, 20.28, 20.13, 20.02, 19.60, 19, 23.66, 15.46, 18.49)     # T. sacculifer data from Foster (2008), Martinez-Botti et al. core top (2015) [including Sanyal et al (2001) refit data] (calcite)
m.Tsac <- 0.82   # mean "m" value for T. sacculifer distribution 
m.Tsacu <- 0.11   # s.d. for "c" value for T. sacculifer distribution
c.Tsac <- 3.94  # mean "c" value for T. sacculifer distribution 
c.Tsacu <- 2.01 # s.d. for "c" value for T. sacculifer distribution

bor.Ouni <- c(18.14, 18.34, 16.93, 17.51, 16.28, 16.61, 17.74, 16.93, 16.28, 16.93, 18.59, 16.28, 18.21, 19.57, 18.21, 19.82, 19.82, 19.76, 18.78)  # Henehan et al. (2016) O. universa data (borate)
for.Ouni <- c(16.17, 16.97, 15.98, 15.99, 15.15, 15.65, 16.29, 15.84, 15.19, 16.08, 17.5, 14.98, 16.62, 18.04, 16.82, 18.57, 18.98, 18.82, 19.47)   # Henehan et al. (2016) O. universa data (calcite)
m.Ouni <- 0.95   # mean "m" value for O. universa distribution 
m.Ouniu <- 0.085   # s.d. for "c" value for O. universa distribution
c.Ouni <- -0.42  # mean "c" value for O. universa distribution 
c.Ouniu <- 1.43 # s.d. for "c" value for O. universa distribution

bor.borate <- c(0, 1) # 1:1 line (borate) 
for.borate <- c(0, 1)  # 1:1 line (calcite)
m.borate <- 1   # mean "m" value for borate distribution 
m.borateu <- 0.1   # s.d. for "c" value for borate distribution
c.borate <- 0  # mean "c" value for borate distribution 
c.borateu <- 1 # s.d. for "c" value for borate distribution


# These parameters will be recorded in the output
parms = c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
          "pco2", "dic","pH", "m.mean", "m.sd", "c.mean", "c.sd")


# Read in datums from file d11Bf, d18Of, Mg/Caf
input.df = "data/Input_data_multisample.xlsx"
input_data = read.xlsx(input.df, sheet = "data")
d11Bf.id <- input_data %>% pull(d11B)
d11Bfu.id <- input_data %>% pull(d11Bsd)
d18Of.id <- input_data %>% pull(d18O)
MgCaf.id <- input_data %>% pull(MgCa)
age <- input_data %>% pull(age)
species <- input_data %>% pull(species)
rows.id <- nrow(input_data)


# Specify input data for each datum and run the inversion 
for(i in 1:rows.id) 
    { 
  species.dat <- species[i]
  
  if (species.dat == "Grub")
  {bor.dat <- bor.Grub
    for.dat <-for.Grub
    m.mean <- m.Grub
    m.sd <- m.Grubu
    c.mean <- c.Grub
    c.sd <- c.Grubu} 
  else if (species.dat == "Tsac")
  {bor.dat <- bor.Tsac
    for.dat <-for.Tsac
    m.mean <- m.Tsac
    m.sd <- m.Tsacu
    c.mean <- c.Tsac
    c.sd <- c.Tsacu} 
  else if (species.dat == "Ouni")
  {bor.dat <- bor.Ouni
    for.dat <-for.Ouni
    m.mean <- m.Ouni
    m.sd <- m.Ouniu
    c.mean <- c.Ouni
    c.sd <- c.Ouniu} 
  else if (species.dat == "borate")
  {bor.dat <- bor.borate
    for.dat <-for.borate
    m.mean <- m.borate
    m.sd <- m.borateu
    c.mean <- c.borate
    c.sd <- c.borateu} 
  
  d11Bf.ld <- d11Bf.id[i]
  d11Bfu.ld <- d11Bfu.id[i]
  MgCaf.ld <- MgCaf.id[i]
  d18Of.ld <- d18Of.id[i]
  
  data = list("d11Bf.data" = d11Bf.ld, "d11Bfu.data" = d11Bfu.ld, "d18Of.data" = d18Of.ld, "mgcaf.data" = MgCaf.ld, 
              "d11Bcb" = bor.dat, "d11Bcfo" = for.dat, "m.mean" = m.mean, "m.sd" = m.sd, "c.mean" = c.mean, 
              "c.sd" = c.sd, "seccal" = seccal, "seccal.u" = seccal.u, "Dd18Oseccal" = Dd18Oseccal, 
              "c.correction" = c.correction)
  
  # Run the inversion
  jout = jags(model.file = "boronPSM_inv.R", parameters.to.save = parms,
                 data = data, inits = NULL, n.chains = 3, n.iter = 1e4,
                 n.burnin = 1e3, n.thin = 10)
  
  # Store datum output in compiled data structure
  multi.out[[i]] <- jout
  parms.out[[i]] <- jout$BUGSoutput$summary
}


