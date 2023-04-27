
# This script takes output from the boron PSM driver, extracts and plots the data in 
# 50% & 95% confidence interval envelopes
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
parms.out <- jout$BUGSoutput$summary
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
# Envelope plots - three
###########################################################################################

fig.font <- "Arial"
fontsize.axislabels <- 10
fontsize.scalelabels <- 8

xmax <- 59
xmin <- 53.2
age.in <- ages.prox/1000

#pco2 in uatm
pco2min <- 0
pco2max <- 3000
# pH 
pHmin <- 7.3
pHmax <- 8.0
# temp in C 
tempCmin <- 27
tempCmax <- 36


pco2.plot <- ggplot(data = pco2.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=(X2.5.*10^6), ymax=(X97.5.*10^6)), fill = "#5F9EA0")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=(X25.*10^6), ymax=(X75.*10^6)), fill = "#008080")  +
  geom_line(mapping = aes(x=age.in, y=(mean*10^6)), color = "#000000") +
  ylim(pco2min, pco2max) +
  xlim(xmax,xmin) +
  labs(x = "age (Ma)", y = expression(paste("paleo-CO"[2]," (",mu,"atm)"))) + 
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25))

pH.plot <- ggplot(data = pH.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#8FBC8F")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#008000")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(pHmin, pHmax) +
  xlim(xmax,xmin) +
  labs(x = "age (Ma)", y = expression(pH)) +
  theme_classic() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25))

tempC.plot <- ggplot(data = tempC.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#FFA07A")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#CD5C5C")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(tempCmin, tempCmax) +
  xlim(xmax,xmin) +
  labs(x = "age (Ma)", y = expression(paste("SST (", degree,"C)"))) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25))


Figure <- ggarrange(pco2.plot, tempC.plot, pH.plot, 
                    labels = c("a.", "b.", "c."),
                    ncol = 1, nrow = 3, align = "hv")
print(Figure)




###########################################################################################
# Envelope plots - six
###########################################################################################

fig.font <- "Helvetica"
fontsize.axislabels <- 11
fontsize.scalelabels <- 9

xmax <- 59
xmin <- 53
age.in <- ages.prox/1000

#pco2 in uatm
pco2min <- 0
pco2max <- 3000
# pH 
pHmin <- 7.2
pHmax <- 8.1
# dic in umol/kgsw
dicmin <- 1500  
dicmax <- 3000
# temp in C 
tempCmin <- 25
tempCmax <- 42
# sal in ppt
salmin <- 32
salmax <- 40
# seawater d11B in per mille
d11Bswmin <- 38
d11Bswmax <- 40

pco2.plot <- ggplot(data = pco2.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=(X2.5.*10^6), ymax=(X97.5.*10^6)), fill = "#5F9EA0")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=(X25.*10^6), ymax=(X75.*10^6)), fill = "#008080")  +
  geom_line(mapping = aes(x=age.in, y=(mean*10^6)), color = "#000000") +
  ylim(pco2min, pco2max) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression('paleo-CO'[2])) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0.1,0.1,-1.5,0.1), "cm")) 

pH.plot <- ggplot(data = pH.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#8FBC8F")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#008000")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(pHmin, pHmax) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression(pH)) +
  theme_classic() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) 

dic.plot <- ggplot(data = dic.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=(X2.5.*10^6), ymax=(X97.5.*10^6)), fill = "#D2B48C")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=(X25.*10^6), ymax=(X75.*10^6)), fill = "#CD853F")  +
  geom_line(mapping = aes(x=age.in, y=(mean*10^6)), color = "#000000") +
  ylim(dicmin, dicmax) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression(paste("DIC (", mu, "mol/kg)"))) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0,0.1,-1.5,0.1), "cm")) 

tempC.plot <- ggplot(data = tempC.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#FFA07A")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#CD5C5C")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(tempCmin, tempCmax) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression(paste("SST (", degree,"C)"))) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0.1,0.1,-1.5,0.1), "cm")) 

sal.plot <- ggplot(data = sal.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#B0C4DE")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#778899")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(salmin, salmax) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression("Salinity (ppt)")) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0,0.1,-1.5,0.1), "cm")) 

d11Bsw.plot <- ggplot(data = d11Bsw.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#D3D3D3")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#999999")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(d11Bswmin, d11Bswmax) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression(paste(delta^11, "B sw (‰)"))) +
  theme_classic() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) 


Figure <- ggarrange(pco2.plot, tempC.plot, dic.plot, sal.plot, pH.plot, d11Bsw.plot,
                   labels = c("a.", "b.", "c.", "d.", "e.", "f."),
                   ncol = 2, nrow = 3, align = "hv")
print(Figure)






