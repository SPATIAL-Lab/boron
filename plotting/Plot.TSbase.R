
# This script takes output from the boron PSM driver, extracts and plots the output for  
# each time step w/ 95% CI error bars
# Dustin T. Harper
# October 31, 2022


###########################################################################################
# Load libraries
###########################################################################################

library(tidyverse)
library(ggpubr)


###########################################################################################
# Extract data for plotting
###########################################################################################

step.vector <- seq(1, n.steps, by=1)
parms.out <- inv.out$BUGSoutput$summary
parms.out <- data.frame(parms.out)

pco2.v <- paste("pco2[", step.vector, "]", sep="")
pco2.out <- parms.out[c(pco2.v),]
pco2.out <- na.omit(pco2.out)
pco2.out <- data.frame(ages.prox,pco2.out)
pco2.out <- na.omit(pco2.out)

dic.v <- paste("dic[", step.vector, "]", sep="")
dic.out <- parms.out[c(dic.v),]
dic.out <- data.frame(ages.prox,dic.out)

pH.v <- paste("pH[", step.vector, "]", sep="")
pH.out <- parms.out[c(pH.v),]
pH.out <- data.frame(ages.prox,pH.out)

sal.v <- paste("sal[", step.vector, "]", sep="")
sal.out <- parms.out[c(sal.v),]
sal.out <- data.frame(ages.prox,sal.out)

tempC.v <- paste("tempC[", step.vector, "]", sep="")
tempC.out <- parms.out[c(tempC.v),]
tempC.out <- data.frame(ages.prox,tempC.out)

xca.v <- paste("xca[", step.vector, "]", sep="")
xca.out <- parms.out[c(xca.v),]
xca.out <- data.frame(ages.prox,xca.out)

xmg.v <- paste("xmg[", step.vector, "]", sep="")
xmg.out <- parms.out[c(xmg.v),]
xmg.out <- data.frame(ages.prox,xmg.out)

xso4.v <- paste("xso4[", step.vector, "]", sep="")
xso4.out <- parms.out[c(xso4.v),]
xso4.out <- data.frame(ages.prox,xso4.out)

d18Osw.v <- paste("d18Osw[", step.vector, "]", sep="")
d18Osw.out <- parms.out[c(d18Osw.v),]
d18Osw.out <- data.frame(ages.prox,d18Osw.out)

d11Bsw.v <- paste("d11Bsw[", step.vector, "]", sep="")
d11Bsw.out <- parms.out[c(d11Bsw.v),]
d11Bsw.out <- data.frame(ages.prox,d11Bsw.out)

c1.out <- parms.out[c("c.1"),]
c2.out <- parms.out[c("c.2"),]
m1.out <- parms.out[c("m.1"),]
m2.out <- parms.out[c("m.2"),]
press.out <- parms.out[c("press"),]


###########################################################################################
# 95% CI plots for each time step 
###########################################################################################


xmax <- 59000 #56200 
xmin <- 53000 #55400 

# Plot parms of interest

ggplot() + 
  geom_errorbar(data = d11Bsw.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = d11Bsw.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(34,44) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "d11Bsw") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = d18Osw.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = d18Osw.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(-4.5,2.5) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "d18Osw") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = xca.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = xca.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(15,25) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "[Ca]") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = xmg.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = xmg.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(30,70) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "[Mg]") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = xso4.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = xso4.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(10,22) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "[SO4]") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = tempC.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = tempC.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(10,40) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "temp (C)") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = sal.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = sal.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(27,43) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "salinity (ppt)") +  
  theme_bw()

ggplot() + 
  geom_errorbar(data = pco2.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = pco2.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(0,0.006) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "pCO2 (atm)") +   
  theme_bw()

ggplot() + 
  geom_errorbar(data = dic.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = dic.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(0.0015, 0.003) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "DIC (mol/kg)") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = pH.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = pH.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(6.8, 8.4) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "pH") +
  theme_bw()



