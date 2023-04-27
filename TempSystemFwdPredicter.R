

# Load libraries 
library(truncnorm)

# Monte Carlo uncertainty propagation for Mg/Caf and d18Of pseudo-measurement predictions
n=1000

tempC = 32.7
sal = 35
d18Osw=-1.2
mgcasw=2.7

d18Oseccal=0.85 # diagenetic end member d18O (per mille)
mgcaswm=5.2 #modern Mg/Casw

n=1000

# Initialize vectors
init.vec <- vector("numeric", n)
Bmod <- init.vec
A <- init.vec
Hp <- init.vec
salcorrco <- init.vec
Bcorr <- init.vec
mgca_corr <- init.vec
mgcaf <- init.vec
sw.sens <- init.vec
indexop <- init.vec
d18Osw.sc <- init.vec
d18Oswpdb <- init.vec
d18Of.pr <- init.vec
d18Of <- init.vec

for (i in 1:n){
Bmod[i] <- rnorm(1,mean=0.38,sd=0.02)
A[i] <- rtruncnorm(n=1,a=0.067,b=0.084,mean=0.0757,sd=0.0045)
Hp[i] <- rnorm(1,mean=0.41,sd=0.1)
salcorrco[i] <- rnorm(1,mean=0.042,sd=0.004)

Bcorr[i] <- ((mgcasw^Hp[i])/(mgcaswm^Hp[i])) * Bmod[i]
mgca_corr[i] <- Bcorr[i]*(exp(A[i]*tempC))
mgcaf[i] <- mgca_corr[i] / (1-(sal-35)*salcorrco[i])

sw.sens[i] <- rnorm(1,mean=0.558,sd=0.03)
indexop[i] <- rnorm(1,mean=0.60,sd=0.025)

d18Osw.sc[i] <- d18Osw + (sw.sens[i]*(sal-35))
d18Oswpdb[i] <- d18Osw.sc[i] - 0.27
d18Of.pr[i] <- d18Oswpdb[i] + ((4.64-(((4.64^2)-(4*0.09*(16.1-tempC)))^0.5))/(2*0.09))
d18Of[i] <- d18Of.pr[i]*(1-indexop[i]) + indexop[i]*d18Oseccal
}

# index of diagenetic overprint
hist(indexop, probability = TRUE, xlim=c(0,1), ylim=c(0,0.02), breaks=seq(0,1,0.1))
lines(density(indexop), col="blue")
abline(v=median(indexop), col="red")
# d18Of
hist(d18Of, probability = TRUE, xlim=c(-5,0), ylim=c(0,2), breaks=seq(-5,0,0.2))
lines(density(d18Of), col="blue")
abline(v=median(d18Of), col="red")
print(median(d18Of))
# Mg/Caf
hist(mgcaf, probability = TRUE, xlim=c(0,20), ylim=c(0,1), breaks=seq(0,20,0.2))
lines(density(mgcaf), col="blue")
abline(v=median(mgcaf), col="red")
print(median(mgcaf))

