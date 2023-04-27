# Temp, pCO2, and pH time series draw plots
ages.ma = ages.prox/1000

plot(ages.ma, jout$BUGSoutput$sims.list$tempC[1,], type="l", axes = FALSE, xlab = "", ylab = expression(paste("SST ",degree,"C")), xlim = c(56.000,55.700), ylim = c(29,40), col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
for (i in 2:500) {
  lines(ages.ma, jout$BUGSoutput$sims.list$tempC[i,], col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
}
  lines(ages.ma, jout$BUGSoutput$median$tempC, col=rgb(red=1, green=0, blue=0), lwd=2)
  axis(2)
  axis(3)

plot(ages.ma, jout$BUGSoutput$sims.list$sal[1,], type="l", axes = FALSE, xlab = "", ylab = "", xlim = c(56.000,55.700),  ylim = c(32,37), col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
for (i in 2:500) {
  lines(ages.ma, jout$BUGSoutput$sims.list$sal[i,], col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
}
  lines(ages.ma, jout$BUGSoutput$median$sal, col=rgb(red=0, green=0, blue=1), lwd=2)
  axis(4)

plot(ages.ma, jout$BUGSoutput$sims.list$d18Osw.sc[1,], type="l", axes = FALSE, xlab = "", ylab = expression(paste(delta^18, "O sw (‰)")), xlim = c(56.000,55.700), ylim = c(-2.5,0.2), col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
for (i in 2:500) {
  lines(ages.ma, jout$BUGSoutput$sims.list$d18Osw.sc[i,], col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
}
  lines(ages.ma, jout$BUGSoutput$median$d18Osw.sc, col=rgb(red=0, green=0.8, blue=0.8), lwd=2)
  axis(2)
  
# Data plots  

plot(ages.ma, prox.in$MgCa, axes =FALSE, ylab = "Mg/Capf (mmol/mol)", xlab = "", xlim = c(56.000,55.700), ylim = c(3.5,6.5), col="navy")
lines(ages.ma, prox.in$MgCa, col="navy")
axis(4)
axis(1)

plot(ages.ma, prox.in$d18O, axes =FALSE, ylab = expression(paste(delta^18, "Opf (‰)")), xlab = "", xlim = c(56.000,55.700), ylim = c(-4.2,-2), col="purple")
lines(ages.ma, prox.in$d18O, col="purple")
axis(2)

plot(ages.ma, prox.in$d13C, xlab = "Age (Ma)", axes =FALSE, ylab = expression(paste(delta^13, "Cpf (‰)")), xlim = c(56.000,55.700), ylim = c(0,3), col="brown")
lines(ages.ma, prox.in$d13C, col="brown")
axis(4)
axis(1)

