library(R2jags)

# These are the observed data, for now just one sample at a time
# First your non-PETM
data = list("mgcaf.data" = 3.5, "d18Of.data" = -0.8, 
            "d11Bf.data" = 15.7)

# This is what parameters you want recorded in the output
parms = c("sal", "tempC", "press", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", 
          "pco2", "dic")

# Run the inversion
jout = jags(model.file = "boronPSM_inv.R", parameters.to.save = parms,
            data = data, inits = NULL, n.chains = 3, n.iter = 1e4,
            n.burnin = 1e3, n.thin = 10)

# Look at the results, first check visually for burn-in effect and 
# convergence
traceplot(jout) #pretty good for most of parms

# Look at the summary stats for the posterior
# looking here at the Rhat and sample size, which we'd like to be 
# < about 1.03 and > several hundred, ideally. Suggests we should
# run more iterations before publishing!
View(jout$BUGSoutput$summary)

# Now let's plot some posterior distributions and compare them with
# the priors to see the influence of the data
plot_g = function(s, r, d){
  dens.pri = density(rgamma(1e6, s, r))
  dens.post = density(d)
  xrange = range(c(dens.pri$x, dens.post$x))
  yrange = range(c(dens.pri$y, dens.post$y))
  plot(dens.pri, type = "l", xlim = xrange, ylim = yrange, lty = 3)
  lines(dens.post)
} #some functions to make plotting easier
plot_n = function(m, s, d){
  dens.pri = density(rnorm(1e6, m, s))
  dens.post = density(d)
  xrange = range(c(dens.pri$x, dens.post$x))
  yrange = range(c(dens.pri$y, dens.post$y))
  plot(dens.pri, type = "l", xlim = xrange, ylim = yrange, lty = 3)
  lines(dens.post)
}
plot_u = function(l, u, d){
  dens.pri = density(runif(1e6, l, u))
  dens.post = density(d)
  xrange = range(c(dens.pri$x, dens.post$x))
  yrange = range(c(dens.pri$y, dens.post$y))
  plot(c(l, l, u, u), c(0, rep(1 / (u-l), 2), 0), type = "l", 
       xlim = xrange, ylim = yrange, lty = 3)
  lines(dens.post)
}

# shorthand so we don't have to type this every time
sl = jout$BUGSoutput$sims.list

# pco2
# posterior distribution looks very different from prior, which means
# the data are informative.
# this sample is providing pretty strong evidence for pCO2 < 1k ppm
plot_u(500e-6, 1500e-6, sl$pco2)

# temperature
# moderate effect of the data, wants to be colder
# than the prior which I set to have a mean of 20 degrees 
# (arbitrarily).
plot_n(20, 5, sl$tempC)

# DIC
# this one shows little impact. There's nothing in the data that 
# 'requires' DIC to be much different than the prior, maybe a slight
# preference for the higher end of the range
plot_u(0.0015, 0.003, sl$dic)

#you can experiment with some more...

# Now run your PETM sample, will use the same prior so differences
# are only a result of the data
data2 = list("mgcaf.data" = 5.5, "d18Of.data" = -1.5, 
            "d11Bf.data" = 14.5)
jout2 = jags(model.file = "boronPSM_inv.R", parameters.to.save = parms,
            data = data2, inits = NULL, n.chains = 3, n.iter = 1e4,
            n.burnin = 1e3, n.thin = 10)

# This one's less well converged on my machine
View(jout2$BUGSoutput$summary)

# So let's run some more iterations
jout2 = update(jout2, n.iter = 4e4, n.thin = 40)
# That's a little bit better...
View(jout2$BUGSoutput$summary)

sl2 = jout2$BUGSoutput$sims.list

# PETM pCO2 does indeed look higher!
plot(density(sl$pco2))
lines(density(sl2$pco2), col = "red")

# And temperature
# Warmer PETM but not that different, temperature is not well
# constrained overall in this version of the model.
plot(density(sl$tempC))
lines(density(sl2$tempC), col = "red")
# Interesting and worth looking into, I notice many of the temp-
# related parameters are not being sampled as effectively as 
# the carbonate system ones...some of this may relate to my 
# pretty general and broad priors