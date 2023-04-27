

# Step increments for sal (ppt) temp (degrees C) and press (bar)
rnd.dig = 0.1

# Ranges of variables over which to evaluate
tempC = seq(10,45,by=rnd.dig)
sal = seq(30,40,by=rnd.dig)
press = seq(0,20,by=rnd.dig)

# Initiate arrays 
temp = c(1:length(tempC))
delV1 = c(1:length(tempC))             
delk1 = c(1:length(tempC))

Ks1m_st = c(1:(length(tempC)*length(sal)))
dim(Ks1m_st) = c((length(tempC)), (length(sal)))
               
Ks1a = c(1:(length(tempC)*length(sal)*length(press)))
dim(Ks1a) = c((length(tempC)), (length(sal)), (length(press)))      

# Constant (cm^3 bar mol^-1 K^-1)
R <- 83.131 

# Nested for loops to calculate 3D array for K*1a (i.e., array of potential K*1 values w/o any major ion effect)
for (i in 1:length(tempC)){
  for (j in 1:length(sal)){
    for (k in 1:length(press)){
      temp[i] <- tempC[i]+273.15
      Ks1m_st[i,j] <-exp(2.83655-2307.1266/temp[i]-1.5529413*(log(temp[i]))-((0.20760841+4.0484/temp[i])*sqrt(sal[j]))+0.0846834*sal[j]-0.00654208*(sal[j]^1.5)+log(1-(0.001005*sal[j])))
      delV1[i] <- (-25.50)+0.1271*tempC[i]
      delk1[i] <- (-3.08/1000)+(0.0877/1000)*tempC[i]
      Ks1a[i,j,k] <- (exp(-((delV1[i]/(R*temp[i]))*press[k])+((0.5*delk1[i])/(R*temp[i]))*press[k]^2))*Ks1m_st[i,j]
}
}
}

# To look up a value in the PSM, pass Ks1a array and rnd.dig to jags then calculate Ks1m[i] as: 
Ks1m[i] <- Ks1a[((round(tempC[ai.prox[i]], digits=rnd.dig)-10)*(1/rnd.dig)), 
                ((round(sal[ai.prox[i]], digits=rnd.dig)-30)*(1/rnd.dig)), 
                (round(press, digits=rnd.dig)*(1/rnd.dig))]
# This currently assumes the time series structure  (ie PSM); it won't work without PSM time series for loop:

