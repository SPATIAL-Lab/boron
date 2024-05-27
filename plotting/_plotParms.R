
# Plots time series, or individual sample posteriors from MCMC output object

load(file = "Harperetal_subm/RevisionApril2024/out/ages.prox.rda")
load(file = "Harperetal_subm/RevisionApril2024/out/LPEE_pHMgCa7.rda")

### INPUT ###
############################################################################################
### MCMC object from MCMC inversion output 
mcmcout <- jout
  
### List of parameters to plot. Parameter names must match saved parameters in MCMC object
parms2plot <- c("sal", "tempC", "d11Bsw", "d18Osw", "pco2", "pH", "xca", "xmg", "xso4")

### Time series ages vector (for time series data) 
ages <- ages.prox

############################################################################################



### PULL DATA FROM MCMC OBJECT AND PLOT SELECTED PARAMETERS ###
############################################################################################
for (i in 1:length(parms2plot)){
  parm <- parms2plot[i]
  parm.sl <- mcmcout[["BUGSoutput"]][["sims.list"]][[parm]]
  parm.med <- mcmcout[["BUGSoutput"]][["median"]][[parm]]
  
  if (length(parm.med) < 2){
    plot(density(parm.sl), xlab = parm, col = rgb(0,0,0, alpha = 0))
    polygon(density(parm.sl), col = adjustcolor("dodgerblue", alpha.f = 0.40))
    abline(v=parm.med, col="black")
    
  } else{
    parm.q025 <- vector(mode = "numeric", length = length(parm.sl[1,]))
    parm.q250 <- vector(mode = "numeric", length = length(parm.sl[1,]))
    parm.q750 <- vector(mode = "numeric", length = length(parm.sl[1,]))
    parm.q975 <- vector(mode = "numeric", length = length(parm.sl[1,]))
    
    for (i in 1:length(parm.sl[1,])){
      parm.q025[i] <- quantile(parm.sl[,i], 0.025)
      parm.q250[i] <- quantile(parm.sl[,i], 0.250)
      parm.q750[i] <- quantile(parm.sl[,i], 0.750)
      parm.q975[i] <- quantile(parm.sl[,i], 0.975)
    }
    
    plot(ages, parm.med, type="l", xlab = "age (ka)", 
         ylab = parm, xlim = rev(range(ages)), ylim = c(min(parm.q025),max(parm.q975)))
    polygon(x = c(ages, rev(ages)), y = c(parm.q025, rev(parm.q975)), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
    polygon(x = c(ages, rev(ages)), y = c(parm.q250, rev(parm.q750)), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.80), border = NA)
    lines(ages, parm.med, col=rgb(red=0, green=0, blue=0, alpha=1), lwd=1)
  }
}  



