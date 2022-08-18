#
# This script takes output from the boron PSM driver, extracts and plots the data
# Dustin T. Harper
# August 9, 2022
#
#

# Load libraries
library(tidyverse)

# Extract data for plotting
step.vector <- seq(1, n.steps, by=1)
parms.out <- jout$BUGSoutput$summary
parms.out <- data.frame(parms.out)

pco2.v <- paste("pco2[", step.vector, "]", sep="")
pco2.out <- parms.out[c(pco2.v),]
pco2.out <- data.frame(ages,pco2.out)

dic.v <- paste("dic[", step.vector, "]", sep="")
dic.out <- parms.out[c(dic.v),]
dic.out <- data.frame(ages,dic.out)

pH.v <- paste("pH[", step.vector, "]", sep="")
pH.out <- parms.out[c(pH.v),]
pH.out <- data.frame(ages,pH.out)

sal.v <- paste("sal[", step.vector, "]", sep="")
sal.out <- parms.out[c(sal.v),]
sal.out <- data.frame(ages,sal.out)

tempC.v <- paste("tempC[", step.vector, "]", sep="")
tempC.out <- parms.out[c(tempC.v),]
tempC.out <- data.frame(ages,tempC.out)

xca.v <- paste("xca[", step.vector, "]", sep="")
xca.out <- parms.out[c(xca.v),]
xca.out <- data.frame(ages,xca.out)

xmg.v <- paste("xmg[", step.vector, "]", sep="")
xmg.out <- parms.out[c(xmg.v),]
xmg.out <- data.frame(ages,xmg.out)

xso4.v <- paste("xso4[", step.vector, "]", sep="")
xso4.out <- parms.out[c(xso4.v),]
xso4.out <- data.frame(ages,xso4.out)

d18Osw.v <- paste("d18Osw[", step.vector, "]", sep="")
d18Osw.out <- parms.out[c(d18Osw.v),]
d18Osw.out <- data.frame(ages,d18Osw.out)

d11Bsw.v <- paste("d11Bsw[", step.vector, "]", sep="")
d11Bsw.out <- parms.out[c(d11Bsw.v),]
d11Bsw.out <- data.frame(ages,d11Bsw.out)

c1.out <- parms.out[c("c.1"),]
c2.out <- parms.out[c("c.2"),]
m1.out <- parms.out[c("m.1"),]
m2.out <- parms.out[c("m.2"),]
press.out <- parms.out[c("press"),]


# Plot parms of interest

ggplot() + 
  geom_pointrange(data = d11Bsw.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(34,44) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "d11Bsw") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = d18Osw.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(-3,2) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "d18Osw") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = xca.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(15,25) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "[Ca]") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = xmg.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(30,70) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "[Mg]") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = xso4.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(10,22) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "[SO4]") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = tempC.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(20,40) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "temp (C)") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = sal.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(33,37) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "salinity (ppt)") +  
  theme_bw()

ggplot() + 
  geom_pointrange(data = pco2.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(0,0.002) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "pCO2 (atm)") +   
  theme_bw()

ggplot() + 
  geom_pointrange(data = dic.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(0.0015, 0.003) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "DIC (mol/kg)") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = pH.out, mapping = aes(x=ages, y=mean, ymin=X25., ymax=X75.)) +
  ylim(7.4, 8.4) +
  xlim(53,59) +
  labs(x= "age (Ma)", y = "pH") +
  theme_bw()


