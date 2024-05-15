# Environment
sal = 35
tempC = 30
temp = tempC + 273.15
d18Osw = -0.75
pH = 8

# Index of diagenetic overprint
indexop = 60/100
d18Oseccal = 0.85

# Sensitivity of seawater d18O to salinity 
sw.sens = 0.558

# Correction for d18O sea surface pH effect (Spero et al., 1997 for symbiont-bearing O. universa Dd18O/DpH = -0.89)
d18O_pHcorr = -0.89 

# Compute d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
d18Osw.sc <- d18Osw + (sw.sens*(sal-35))
alpha.ccw <- exp(((17.57*10^3)/temp - 29.13) / 10^3)
epsilon.ccw <- (alpha.ccw-1)*10^3 
d18Ofsmow <- d18Osw.sc + epsilon.ccw 
d18Ofpdb <- 0.97001 * d18Ofsmow - 29.99
d18Ofdiag <- d18Ofpdb - (pH - 8)*d18O_pHcorr
d18Of <- d18Ofdiag*(1-indexop) + d18Oseccal*indexop

print(d18Of)
