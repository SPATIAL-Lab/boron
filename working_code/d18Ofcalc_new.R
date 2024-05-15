
d18Osw = -0.7
temp = 310

alpha.ccw[i] <-  exp(((17.57*10^3)/tempC[ai.prox[i]] - 29.13) / 10^3)
epsilon.ccw[i] <- (alpha.ccw[i]-1)*10^3 
d18Ofsmow[i] <- d18Osw.sc[i] + epsilon.ccw[i] 
d18Of[i] = 0.97001 * d18Ofsmow[i] - 29.99



##########################

MgCa_satdec ~ dnorm(0.054, 1/0.008^2)

ifelse(Dco3 < 21.3, MgCa_satcorr = Dco3-21.3*MgCa_satdec)

MgCaf <- MgCa00 + MgCa_satcorr