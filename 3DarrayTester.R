
library(rjags)
library(R2jags)
# Parameters to save in inversion output
parms = c("A", "tempC") 

# Data to pass
multiarray = c(0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095)
dim(multiarray) = c(2,2,2)      
mgca.data <- c(3.2,3.0,3.3,3.4)
ai.mgca <- c(1,2,3,4)
data.pass = list("multiarray", "mgca.data", "ai.mgca")


#############################################################
# Model 1 RUNS; 4 nodes, graph size = 42
model.string = "model {

  for (i in 1:length(ai.mgca)){
  mgca.data[i] ~ dnorm(mgcaf[ai.mgca[i]], 0.08)
}

  x.index <- 1
  y.index <- 1
  z.index <- 2
  
  for (i in 1:length(ai.mgca)){
  A[i] <- multiarray[x.index,  y.index, z.index]
  mgcaf[i] <- 0.41*(exp(A[i]*tempC[ai.mgca[i]]))
  }
  
  tempC[1] ~ dnorm(25,0.25)
  for (i in 2:length(ai.mgca)){
    tempC.sig[i] ~ dnorm(0, 1)
    tempC[i] = tempC[i-1] + tempC.sig[i]
}
}"
#############################################################

#############################################################
# Model 2 RUNS; 4 nodes, graph size = 42
model.string = "model {

  for (i in 1:length(ai.mgca)){
  mgca.data[i] ~ dnorm(mgcaf[ai.mgca[i]], 0.08)
}
  
  for (i in 1:length(ai.mgca)){
  x.index[i] <- 1
  y.index[i] <- 1
  z.index[i] <- 2
  
  A[i] <- multiarray[x.index[i],  y.index[i], z.index[i]]
  mgcaf[i] <- 0.41*(exp(A[i]*tempC[ai.mgca[i]]))
  }
  
  tempC[1] ~ dnorm(25,0.25)
  for (i in 2:length(ai.mgca)){
    tempC.sig[i] ~ dnorm(0, 1)
    tempC[i] = tempC[i-1] + tempC.sig[i]
}
}"
#############################################################

#############################################################
# Model 3 RUNS; 4 nodes, graph size = 57
model.string = "model {

  for (i in 1:length(ai.mgca)){
  mgca.data[i] ~ dnorm(mgcaf[ai.mgca[i]], 0.08)
}
  
  for (i in 1:length(ai.mgca)){
  x.index[i] <- round(tempC[ai.mgca[i]]/tempC[ai.mgca[i]])
  y.index[i] <- round(tempC[ai.mgca[i]]/tempC[ai.mgca[i]])+1
  z.index[i] <- round(tempC[ai.mgca[i]]/tempC[ai.mgca[i]])
  
  A[i] <- multiarray[x.index[i],  y.index[i], z.index[i]]
  mgcaf[i] <- 0.41*(exp(A[i]*tempC[ai.mgca[i]]))
  }
  
  tempC[1] ~ dnorm(25,0.25)
  for (i in 2:length(ai.mgca)){
    tempC.sig[i] ~ dnorm(0, 1)
    tempC[i] = tempC[i-1] + tempC.sig[i]
}
}"
#############################################################

writeLines(model.string, con = "3DarrayTester.model.txt")

jout = jags(data=data.pass, model.file = "3DarrayTester.model.txt", parameters.to.save = parms,
                          inits = NULL, n.chains = 3, n.iter = 2000,
                          n.burnin = 500, n.thin = 1)


