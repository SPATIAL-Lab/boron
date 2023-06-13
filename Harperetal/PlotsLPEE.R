# Temp, pCO2, and pH time series draw plots

plot(ages.prox, jout$BUGSoutput$sims.list$tempC[1,], type="l", axes = FALSE, xlab = "", ylab = expression(paste("SST ",degree,"C")), xlim = rev(range(ages.prox)), ylim = c(28,41), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
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
  
# Climate sensitivity density plot
  
LPEE.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,1:25]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,430:435])
LPEE.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,430:435])
LPEE.pc.pco2 <- LPEE.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,1:25]) * 100
LPEE.CS <- LPEE.DSST / LPEE.pc.pco2 * 100
LPEE.CS <- subset(LPEE.CS, LPEE.CS<=10 & LPEE.CS >= 0)

PETM.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,202:209]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,230:233])
PETM.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233])
PETM.pc.pco2 <- PETM.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]) * 100
PETM.CS <- PETM.DSST / PETM.pc.pco2 * 100
PETM.CS <- subset(PETM.CS, PETM.CS<=10 & PETM.CS >= 0)

ETM2.DSST <- rowMeans(jout$BUGSoutput$sims.list$tempC[,359:360]) - rowMeans(jout$BUGSoutput$sims.list$tempC[,369:373])
ETM2.Dpco2 <- rowMeans(jout$BUGSoutput$sims.list$pco2[,359:360]) - rowMeans(jout$BUGSoutput$sims.list$pco2[,369:373])
ETM2.pc.pco2 <- ETM2.Dpco2 / rowMeans(jout$BUGSoutput$sims.list$pco2[,358:360]) * 100
ETM2.CS <- ETM2.DSST / ETM2.pc.pco2 * 100
ETM2.CS <- subset(ETM2.CS, ETM2.CS<=10 & ETM2.CS >= 0)

plot(density(LPEE.CS), xlab = expression(paste("ESCS (",degree,"C / CO"[2], " doubling)")), col="black", ylim=c(0,0.55), xlim=c(1,10), lwd=1.5)
lines(density(PETM.CS), col="red", lwd=1.5)
lines(density(ETM2.CS), col="blue", lwd=1.5)
abline(v=median(PETM.CS), col="red")
abline(v=median(ETM2.CS), col="blue")
abline(v=median(LPEE.CS), col="black")

# Print averages and confidence intervals of interest

print(c(median(PETM.CS), quantile(PETM.CS, 0.025), quantile(PETM.CS, 0.975)))
print(c(mean(-PETM.DSST), quantile(-PETM.DSST, 0.025), quantile(-PETM.DSST, 0.975)))
print(c(mean(-PETM.Dpco2), quantile(-PETM.Dpco2, 0.025), quantile(-PETM.Dpco2, 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,202:209]), 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,230:233]), 0.975)))

print(c(median(ETM2.CS), quantile(ETM2.CS, 0.025), quantile(ETM2.CS, 0.975)))
print(c(mean(-ETM2.DSST), quantile(-ETM2.DSST, 0.025), quantile(-ETM2.DSST, 0.975)))
print(c(mean(-ETM2.Dpco2), quantile(-ETM2.Dpco2, 0.025), quantile(-ETM2.Dpco2, 0.975)))

print(c(median(LPEE.CS), quantile(LPEE.CS, 0.025), quantile(LPEE.CS, 0.975)))
print(c(mean(-LPEE.DSST), quantile(-LPEE.DSST, 0.025), quantile(-LPEE.DSST, 0.975)))
print(c(mean(-LPEE.Dpco2), quantile(-LPEE.Dpco2, 0.025), quantile(-LPEE.Dpco2, 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,1:19])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,1:19]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,1:19]), 0.975)))
print(c(mean(rowMeans(jout$BUGSoutput$sims.list$pco2[,430:433])), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,430:433]), 0.025), quantile(rowMeans(jout$BUGSoutput$sims.list$pco2[,430:433]), 0.975)))



