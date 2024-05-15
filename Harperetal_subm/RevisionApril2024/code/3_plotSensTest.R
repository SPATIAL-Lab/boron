
# Plots time series posteriors from MCMC output object for senstivity tests

### INPUT ###
############################################################################################
### Load MCMC objects from MCMC inversion output 
load(file = "Harperetal_subm/RevisionApril2024/out_primary/LPEE_prim.rda")
jout_prim <- jout
rm(jout)
  
load(file = "Harperetal_subm/RevisionApril2024/out_senstest/LPEE_MgCasw.rda")
jout_st <- jout
rm(jout)


### List of time series parameters to plot. Parameter names must match saved parameters in MCMC object
parms2plot <- c("sal", "tempC", "pco2", "pH")

### Time series ages vector (for time series data) 
load(file = "Harperetal_subm/RevisionApril2024/out_primary/ages.prox.rda")
ages <- ages.prox

############################################################################################



### PULL DATA FROM MCMC OBJECT AND PLOT SELECTED PARAMETERS ###
############################################################################################
for (i in 1:length(parms2plot)){
  parm <- parms2plot[i]
  
  parm.sl <- jout_prim[["BUGSoutput"]][["sims.list"]][[parm]]
  parm.med <- jout_prim[["BUGSoutput"]][["median"]][[parm]]
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
  
  parm.sl.st <- jout_st[["BUGSoutput"]][["sims.list"]][[parm]]
  parm.med.st <- jout_st[["BUGSoutput"]][["median"]][[parm]]
  parm.q025.st <- vector(mode = "numeric", length = length(parm.sl.st[1,]))
  parm.q250.st <- vector(mode = "numeric", length = length(parm.sl.st[1,]))
  parm.q750.st <- vector(mode = "numeric", length = length(parm.sl.st[1,]))
  parm.q975.st <- vector(mode = "numeric", length = length(parm.sl.st[1,]))
    
  for (i in 1:length(parm.sl.st[1,])){
    parm.q025.st[i] <- quantile(parm.sl.st[,i], 0.025)
    parm.q250.st[i] <- quantile(parm.sl.st[,i], 0.250)
    parm.q750.st[i] <- quantile(parm.sl.st[,i], 0.750)
    parm.q975.st[i] <- quantile(parm.sl.st[,i], 0.975)
  }
    
    plot(ages, parm.med, type="l", xlab = "age (ka)", 
         ylab = parm, xlim = rev(range(ages)), ylim = c(min(parm.q025),max(parm.q975)))
    polygon(x = c(ages, rev(ages)), y = c(parm.q025, rev(parm.q975)), 
            col =  adjustcolor("black", alpha.f = 0.40), border = NA)
    polygon(x = c(ages, rev(ages)), y = c(parm.q250, rev(parm.q750)), 
            col =  adjustcolor("black", alpha.f = 0.60), border = NA)
    lines(ages, parm.med, col=rgb(red=0, green=0, blue=0, alpha=1), lwd=1)
    polygon(x = c(ages, rev(ages)), y = c(parm.q025.st, rev(parm.q975.st)), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    polygon(x = c(ages, rev(ages)), y = c(parm.q250.st, rev(parm.q750.st)), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
    lines(ages, parm.med.st, col=adjustcolor("dodgerblue", alpha.f = 1), lwd=1.5)
}  


