
# Plots time series posteriors from MCMC output object for sensitivity tests
# Note that full inversion output files used here ('jout') are too big to include
# in the repository. Thus, individual sensitivity test inversions must be re-run 
# and inversion output saved in "out_..." folders for this script to function.

### INPUT ###
############################################################################################
### Load MCMC objects from MCMC inversion output 
load(file = "Harperetal_resubm/out_primary/LPEE_prim.rda")
jout_prim <- jout
rm(jout)
  
load(file = "Harperetal_resubm/out_senstest/LPEE_MgCasw.rda")
jout_MgCasw <- jout
rm(jout)

load(file = "Harperetal_resubm/out_senstest/LPEE_pHMgCa3.rda")
jout_pHMgCa3 <- jout
rm(jout)

load(file = "Harperetal_resubm/out_senstest/LPEE_pHMgCa7.rda")
jout_pHMgCa7 <- jout
rm(jout)

load(file = "Harperetal_resubm/out_senstest/LPEE_nopHd18O.rda")
jout_nopHd18O <- jout
rm(jout)

load(file = "Harperetal_resubm/out_senstest/LPEE_CO3.rda")
jout_CO3 <- jout
rm(jout)

load(file = "Harperetal_resubm/out_senstest/LPEE_CCvar.rda")
jout_CCvar <- jout
rm(jout)

### List of time series parameters to plot. Parameter names must match saved parameters in MCMC object
parms2plot <- c("sal", "tempC", "pco2", "pH")

### Time series ages vector (for time series data) 
load(file = "Harperetal_resubm/out_primary/ages.prox.rda")
ages <- ages.prox/1e3

############################################################################################



### PULL DATA FROM MCMC OBJECT AND PLOT SELECTED PARAMETERS ###
############################################################################################
for (i in 1:length(parms2plot)){
  parm <- parms2plot[i]
  
  parm.sl <- jout_prim[["BUGSoutput"]][["sims.list"]][[parm]]
  parm.med <- jout_prim[["BUGSoutput"]][["median"]][[parm]]
  parm.med.pHMgCa3 <- jout_pHMgCa3[["BUGSoutput"]][["median"]][[parm]]
  parm.med.pHMgCa7 <- jout_pHMgCa7[["BUGSoutput"]][["median"]][[parm]]
  parm.med.nopHd18O <- jout_nopHd18O[["BUGSoutput"]][["median"]][[parm]]
  parm.med.MgCasw <- jout_MgCasw[["BUGSoutput"]][["median"]][[parm]]
  parm.med.CO3 <- jout_CO3[["BUGSoutput"]][["median"]][[parm]]
  parm.med.CCvar <- jout_CCvar[["BUGSoutput"]][["median"]][[parm]]
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
  
  
    plot(ages, parm.med, type="l", xlab = "age (Ma)", 
         ylab = parm, xlim = rev(range(ages)), ylim = c(min(parm.q025),max(parm.q975)))
    polygon(x = c(ages, rev(ages)), y = c(parm.q025, rev(parm.q975)), 
            col =  adjustcolor("black", alpha.f = 0.40), border = NA)
    polygon(x = c(ages, rev(ages)), y = c(parm.q250, rev(parm.q750)), 
            col =  adjustcolor("black", alpha.f = 0.60), border = NA)
    lines(ages, parm.med, col=rgb(red=0, green=0, blue=0, alpha=1), lwd=1)
    lines(ages, parm.med.pHMgCa3, col=adjustcolor(rgb(230,159,0, maxColorValue=255), alpha.f = 1), lwd=1)
    lines(ages, parm.med.pHMgCa7, col=adjustcolor(rgb(86,180,233, maxColorValue=255), alpha.f = 1), lwd=1)
    lines(ages, parm.med.nopHd18O, col=adjustcolor(rgb(0,158,115, maxColorValue=255), alpha.f = 1), lwd=1)
    lines(ages, parm.med.MgCasw, col=adjustcolor(rgb(0,114,178, maxColorValue=255), alpha.f = 1), lwd=1)
    lines(ages, parm.med.CO3, col=adjustcolor(rgb(213,94,0, maxColorValue=255), alpha.f = 1), lwd=1)
    lines(ages, parm.med.CCvar, col=adjustcolor(rgb(204,121,167, maxColorValue=255), alpha.f = 1), lwd=1)
    legend("topleft", legend=c("3% Mg/Ca pH sens.", "7% Mg/Ca pH sens.", expression(No ~ δ^18*O ~ pH ~ sens.), 
                          expression("Constant Mg/Ca"[sw]), expression("Mg/Ca" ~ CO[3]^"2-" ~ adjusted), "TA as 2nd variable"),  
           fill = c(rgb(230,159,0, maxColorValue=255),rgb(86,180,233, maxColorValue=255),
                    rgb(0,158,115, maxColorValue=255), rgb(0,114,178, maxColorValue=255),
                    rgb(213,94,0, maxColorValue=255),rgb(204,121,167, maxColorValue=255)))
}  


