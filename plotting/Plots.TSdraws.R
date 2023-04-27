# Temp, pCO2, and pH time series draw plots

plot(ages.prox, jout$BUGSoutput$sims.list$tempC[1,], type="l", axes = FALSE, xlab = "", ylab = expression(paste("SST ",degree,"C")), xlim = rev(range(ages.prox)), ylim = c(30,40), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages.prox, jout$BUGSoutput$sims.list$tempC[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
  lines(ages.prox, jout$BUGSoutput$median$tempC, col=rgb(red=1, green=0, blue=0), lwd=1.5)
  axis(4)

  # plot(ages.prox, jout$BUGSoutput$median$tempC, type="l", axes = FALSE, xlab = "", ylab = expression(paste("SST ",degree,"C")), xlim = rev(range(ages.prox)), ylim = c(25,40), col=rgb(red=1, green=0, blue=0), lwd=1.5)
  # temp.quant <- apply(jout$BUGSoutput$sims.list$tempC, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE)
  # lines(ages.prox, temp.quant[1,], col=rgb(red=1, green=0, blue=0), lwd=0.5, lty="dotted")
  # lines(ages.prox, temp.quant[2,], col=rgb(red=1, green=0, blue=0), lwd=0.5, lty="dotted")
  # axis(4)
  # 
  
  
plot(ages.prox, jout$BUGSoutput$sims.list$pco2[1,], type="l", axes = FALSE, xlab = "", ylab = "", xlim = rev(range(ages.prox)), ylim = c(0.0002,0.003), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages.prox, jout$BUGSoutput$sims.list$pco2[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
  lines(ages.prox, jout$BUGSoutput$median$pco2, col=rgb(red=0.3, green=0.4, blue=1),lwd=2)
  axis(2)
  axis(3)

plot(ages.prox, jout$BUGSoutput$sims.list$pH[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("pH")), xlim = rev(range(ages.prox)), ylim = c(7.3,8.1), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages.prox, jout$BUGSoutput$sims.list$pH[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
  lines(ages.prox, jout$BUGSoutput$median$pH, col=rgb(red=0.3, green=0.8, blue=0.4), lwd=2)
  axis(1)
  axis(2)
  
# Climate sensitivity density plot
  
LPEE.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,1:19]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,430:433])
LPEE.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,1:19]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,430:433])
LPEE.pc.pco2 <- LPEE.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,1:19]) * 100
LPEE.CS <- LPEE.DSST / LPEE.pc.pco2 * 100
LPEE.CS <- subset(LPEE.CS, LPEE.CS<=10 & LPEE.CS >= 0)

# PCIM.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,5:15]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,75:77])
# PCIM.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,5:15]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,75:77])
# PCIM.pc.pco2 <- PCIM.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,5:15]) * 100
# PCIM.CS <- PCIM.DSST / PCIM.pc.pco2 * 100
# PCIM.CS <- subset(PCIM.CS, PCIM.CS<=12 & PCIM.CS >= 0)

PETM.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,202:209]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,230:233])
PETM.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233])
PETM.pc.pco2 <- PETM.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]) * 100
PETM.CS <- PETM.DSST / PETM.pc.pco2 * 100
PETM.CS <- subset(PETM.CS, PETM.CS<=10 & PETM.CS >= 0)

ETM2.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,359:360]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,370:371])
ETM2.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,370:371])
ETM2.pc.pco2 <- ETM2.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360]) * 100
ETM2.CS <- ETM2.DSST / ETM2.pc.pco2 * 100
ETM2.CS <- subset(ETM2.CS, ETM2.CS<=10 & ETM2.CS >= 0)

plot(density(LPEE.CS), xlab = expression(paste("ESCS (",degree,"C / CO"[2], " doubling)")), col="black", ylim=c(0,0.55), xlim=c(1,10), lwd=1.5)
lines(density(PETM.CS), col="red", lwd=1.5)
lines(density(ETM2.CS), col="blue", lwd=1.5)
abline(v=median(PETM.CS), col="red")
abline(v=median(ETM2.CS), col="blue")
abline(v=median(LPEE.CS), col="black")
#lines(density(PCIM.CS), col="blue")
#abline(v=median(PCIM.CS), col="blue")