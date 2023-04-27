
# This script takes output from the MagOxy PSM driver, extracts and plots the data in 
# 50% & 95% confidence interval envelopes
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

d18Osw.sc.v <- paste("d18Osw.sc[", step.vector, "]", sep="")
d18Osw.sc.out <- parms.out[c(d18Osw.sc.v),]
d18Osw.sc.out <- data.frame(ages.prox,d18Osw.sc.out)

d18Osw.v <- paste("d18Osw[", step.vector, "]", sep="")
d18Osw.out <- parms.out[c(d18Osw.v),]
d18Osw.out <- data.frame(ages.prox,d18Osw.out)


###########################################################################################
# Envelope plots - three
###########################################################################################

fig.font <- "Arial"
fontsize.axislabels <- 10
fontsize.scalelabels <- 8

xmax <- 56
xmin <- 55.7
age.in <- ages.prox/1000

# temp in C 
tempCmin <- 26
tempCmax <- 40
# sal in ppt
salmin <- 32
salmax <- 38
# seawater d18O in per mille
d18Oswmin <- -2
d18Oswmax <- 1

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

d18Osw.sc.plot <- ggplot(data = d18Osw.sc.out) + 
  geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#D3D3D3")  +
  geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#999999")  +
  geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
  ylim(d18Oswmin, d18Oswmax) +
  xlim(xmax,xmin) +
  labs(x= "age (Ma)", y = expression(paste(delta^18, "O sw (‰)"))) +
  theme_classic() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
        plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) 

# d18Osw.plot <- ggplot(data = d18Osw.out) + 
#   geom_ribbon(mapping = aes(x=age.in, ymin=X2.5., ymax=X97.5.), fill = "#D3D3D3")  +
#   geom_ribbon(mapping = aes(x=age.in, ymin=X25., ymax=X75.), fill = "#999999")  +
#   geom_line(mapping = aes(x=age.in, y=mean), color = "#000000") +
#   ylim(d18Oswmin, d18Oswmax) +
#   xlim(xmax,xmin) +
#   labs(x= "age (Ma)", y = expression(paste(delta^18, "O glob.sw (‰)"))) +
#   theme_classic() +
#   theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
#         axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
#         axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
#         panel.grid.major.y = element_line(color = "#B0C4DE", size = 0.25),
#         plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) 


Figure <- ggarrange(tempC.plot, sal.plot, d18Osw.sc.plot,
                   labels = c("a.", "b.", "c."),
                   ncol = 1, nrow = 3, align = "hv")
print(Figure)






