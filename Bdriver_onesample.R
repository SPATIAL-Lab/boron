library(R2jags)

data = list("mgcaf.data" = 3.5, "d18Of.data" = 0.6, "d11Bf.data" = 15.8)

parms = c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
          "pco2", "dic")

jout = jags(model.file = "boronPSM_inv.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 3, n.iter = 1000,
            n.burnin = 100, n.thin = 1)
