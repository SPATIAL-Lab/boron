#
# This contains the R code to plot Figures 2, 3, and S5 in Harper et al., submission to PNAS 
# 


###########################################################################################
# Temp, pCO2, and pH time series draw plots
###########################################################################################
plot(ages.prox, jout$BUGSoutput$sims.list$tempC[1,], type="l", axes = FALSE, xlab = "", ylab = expression(paste("SST ",degree,"C")), xlim = rev(range(ages.prox)), ylim = c(30,43), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages.prox, jout$BUGSoutput$sims.list$tempC[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
  lines(ages.prox, jout$BUGSoutput$median$tempC, col="red", lwd=1.5)
  axis(4)

plot(ages.prox, jout$BUGSoutput$sims.list$pco2[1,], type="l", axes = FALSE, xlab = "", ylab = "", xlim = rev(range(ages.prox)), ylim = c(0.0004,0.0033), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages.prox, jout$BUGSoutput$sims.list$pco2[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
  lines(ages.prox, jout$BUGSoutput$median$pco2, col="deepskyblue2",lwd=2)
  axis(2)
  axis(3)

plot(ages.prox, jout$BUGSoutput$sims.list$pH[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("pH")), xlim = rev(range(ages.prox)), ylim = c(7.3,8.1), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages.prox, jout$BUGSoutput$sims.list$pH[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
  lines(ages.prox, jout$BUGSoutput$median$pH, col="darkgoldenrod2", lwd=2)
  axis(1)
  axis(2)
  
  
###########################################################################################  
# Climate sensitivity density plot
###########################################################################################  
LPEE.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,1:25]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,431:435])
LPEE.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,431:435])
LPEE.pc.pco2 <- LPEE.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25]) * 100
LPEE.CS <- LPEE.DSST / LPEE.pc.pco2 * 100
LPEE.CS <- subset(LPEE.CS, LPEE.CS<=10 & LPEE.CS >= 0)

PCIM.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,7:12]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,75:79])
PCIM.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,7:12]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,75:79])

PETM.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,202:209]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,230:233])
PETM.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233])
PETM.pc.pco2 <- PETM.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]) * 100
PETM.CS <- PETM.DSST / PETM.pc.pco2 * 100
PETM.CS <- subset(PETM.CS, PETM.CS<=10 & PETM.CS >= 0)

ETM2.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,359:360]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,383:384])
ETM2.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,368:384])
ETM2.pc.pco2 <- ETM2.Dpco2 / min(jout$BUGSoutput$sims.list$pco2[,359:360]) * 100
ETM2.CS <- ETM2.DSST / ETM2.pc.pco2 * 100
ETM2.CS <- subset(ETM2.CS, ETM2.CS<=10 & ETM2.CS >= 0)

plot(density(LPEE.CS), xlab = expression(paste("ESCS (",degree,"C / CO"[2], " doubling)")), col = rgb(0,0,0, alpha = 0), ylim=c(0,0.55), xlim=c(1,10), lwd=1.5)
polygon(density(LPEE.CS), col = rgb(0,0,0, alpha = 0.4))
polygon(density(PETM.CS), col = rgb(1,0,0, alpha = 0.4))
polygon(density(ETM2.CS), col = rgb(0,0,1, alpha = 0.4))
abline(v=median(PETM.CS), col="red")
abline(v=median(ETM2.CS), col="blue")
abline(v=median(LPEE.CS), col="black")


###########################################################################################
# Print averages and confidence intervals of interest
###########################################################################################

# PETM
print(c(median(PETM.CS), quantile(PETM.CS, 0.025), quantile(PETM.CS, 0.975)))
print(c(mean(-PETM.DSST), quantile(-PETM.DSST, 0.025), quantile(-PETM.DSST, 0.975)))
print(c(mean(-PETM.Dpco2), quantile(-PETM.Dpco2, 0.025), quantile(-PETM.Dpco2, 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]), 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233]), 0.975)))

# ETM-2
print(c(median(ETM2.CS), quantile(ETM2.CS, 0.025), quantile(ETM2.CS, 0.975)))
print(c(mean(-ETM2.DSST), quantile(-ETM2.DSST, 0.025), quantile(-ETM2.DSST, 0.975)))
print(c(mean(-ETM2.Dpco2), quantile(-ETM2.Dpco2, 0.025), quantile(-ETM2.Dpco2, 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360]), 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,368:384])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,368:384]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,368:384]), 0.975)))

# PCIM
print(c(mean(-PCIM.DSST), quantile(-PCIM.DSST, 0.025), quantile(-PCIM.DSST, 0.975)))
print(c(mean(-PCIM.Dpco2), quantile(-PCIM.Dpco2, 0.025), quantile(-PCIM.Dpco2, 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,7:12])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,7:12]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,7:12]), 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,75:79])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,75:79]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,75:79]), 0.975)))

# LPEE
print(c(median(LPEE.CS), quantile(LPEE.CS, 0.025), quantile(LPEE.CS, 0.975)))
print(c(mean(-LPEE.DSST), quantile(-LPEE.DSST, 0.025), quantile(-LPEE.DSST, 0.975)))
print(c(mean(-LPEE.Dpco2), quantile(-LPEE.Dpco2, 0.025), quantile(-LPEE.Dpco2, 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25]), 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,431:435])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,431:435]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,431:435]), 0.975)))

# Latest Paleocene
LP.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,170:171]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,182:183])
LP.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,170:171]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,182:183])
print(c(mean(-LP.DSST), quantile(-LP.DSST, 0.025), quantile(-LP.DSST, 0.975)))
print(c(mean(-LP.Dpco2), quantile(-LP.Dpco2, 0.025), quantile(-LP.Dpco2, 0.975)))


###########################################################################################
# LOSCAR Supplemental Figure - 5800 Gt release of -35 per mille C
###########################################################################################
LOSCAR <- readRDS("Harperetal/data/LOSCAR_PETM5800.rds")

# read in and format PSM output 
jout_summ <- jout$BUGSoutput$summary
step.vector <- seq(1, n.steps, by=1)
jout_summ <- data.frame(jout$BUGSoutput$summary)

pco2.v <- paste("pco2[", step.vector, "]", sep="")
pco2.out <- jout_summ[c(pco2.v),]
pco2.out <- na.omit(pco2.out)
pco2.out <- data.frame(ages.prox, pco2.out$mean*10^6, pco2.out$X2.5.*10^6, pco2.out$X97.5.*10^6)
names(pco2.out) <- c("age","mean", "low", "high") 
pco2.out$median <- jout$BUGSoutput$median$pco2*10^6
pco2.out <- pco2.out[pco2.out$age < 56000 & pco2.out$age > 55700, ]

pH.v <- paste("pH[", step.vector, "]", sep="")
pH.out <- jout_summ[c(pH.v),]
pH.out <- na.omit(pH.out)
pH.out <- data.frame(ages.prox, pH.out$mean, pH.out$X2.5., pH.out$X97.5.)
names(pH.out) <- c("age","mean", "low", "high") 
pH.out$median <- jout$BUGSoutput$median$pH
pH.out <- pH.out[pH.out$age < 56000 & pH.out$age > 55700, ]

# define LOSCAR C emissions for release of 4100 GtC over 6kyr followed by 1700 GtC over 140kyr simulation
emist <- c(-25,0,6,140,250)
emisr <- c(0,0.683,0.012,0,0)
emis <- data.frame(emist, emisr)

# set figure aesthetics 
fig.font <- "Arial"
fontsize.axislabels <- 12
fontsize.scalelabels <- 12

# make plots
library(tidyverse)
library(ggpubr)

emis.plot <- ggplot() + 
  geom_step(data = emis, aes(x=emist, y=emisr), color = "black") +
  ylim(0, 0.75) +
  xlim(-50, 250) +
  labs(x = "time (kyr)", y = expression("C emissions (Gt / year)")) + 
  theme_bw() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"))

co2.plot <- ggplot() + 
  geom_ribbon(data = pco2.out, aes(x=((age-55950)*-1), ymin=low, ymax=high), fill = "grey") +
  geom_line(data = pco2.out, aes(x=((age-55950)*-1), y=pco2.out$median), color ='deepskyblue2', size = 1) +
  geom_line(data = LOSCAR, aes(x=time_kyr, y=co2a), color = "black", size = 1) +
  ylim(800, 3000) +
  xlim(-50, 250) +
  labs(x = "time (kyr)", y = expression(paste("paleo-CO"[2]," (ppm)"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"))

pH.plot <- ggplot() + 
  geom_ribbon(data = pH.out, aes(x=((age-55950)*-1), ymin=low, ymax=high), fill = "grey") +
  geom_line(data = pH.out, aes(x=((age-55950)*-1), y=pH.out$median), color ='darkgoldenrod2', size = 1) +
  geom_line(data = LOSCAR, aes(x=time_kyr, y=pH_surfpac), color = "black", size = 1) +
  ylim(7.33,7.75) +
  xlim(-50, 250) +
  labs(x = "time (kyr)", y = expression(pH)) +
  theme_bw() +
  theme(axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels, color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"),
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000"))

Figure <- ggarrange(emis.plot, co2.plot, pH.plot,
                    labels = c("a.", "b.", "c."),
                    ncol = 1, nrow = 3, align = "hv")
print(Figure)





