
# This script takes output from the MagOxy PSM driver, extracts and plots the output for  
# each time step w/ 95% CI error bars
# Dustin T. Harper
# February 3, 2023


###########################################################################################
# Load libraries
###########################################################################################

library(tidyverse)
library(ggpubr)


###########################################################################################
# Extract data for plotting
###########################################################################################

step.vector <- seq(1, n.steps, by=1)
parms.out <- jout$BUGSoutput$summary
parms.out <- data.frame(parms.out)

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

d18Osw.v <- paste("d18Osw[", step.vector, "]", sep="")
d18Osw.out <- parms.out[c(d18Osw.v),]
d18Osw.out <- data.frame(ages.prox,d18Osw.out)

d18Osw.sc.v <- paste("d18Osw.sc[", step.vector, "]", sep="")
d18Osw.sc.out <- parms.out[c(d18Osw.sc.v),]
d18Osw.sc.out <- data.frame(ages.prox,d18Osw.sc.out)

###########################################################################################
# 95% CI plots for each time step 
###########################################################################################


xmax <- 56200 
xmin <- 55400 

# Plot parms of interest

ggplot() + 
  geom_errorbar(data = d18Osw.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = d18Osw.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(-4.5,2.5) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "d18Osw.global") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = d18Osw.sc.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray") +
  geom_point(data = d18Osw.sc.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(-4.5,2.5) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = "d18Osw.local") +
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




