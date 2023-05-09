

############################################################################################
#    GENERATE CARBONATE CHEMISTRY EQUIL CONSTANT LOOKUP ARRAYS - follows Zeebe & Wolf-Gladrow 2001
############################################################################################    

# Set upper and lower STP bounds for equil constant array 
tempC.lb = 10
tempC.ub = 45
sal.lb = 30
sal.ub = 40
press.lb = 0
press.ub= 20

# Step increments for sal (ppt) temp (degrees C) and press (bar) and determine rounding digits
inc = 0.1
rnd.dig <- as.numeric(match(TRUE, round(inc, 1:20) == inc))

# Ranges of variables over which to evaluate
tempC = seq(tempC.lb, tempC.ub, by=inc)
sal = seq(sal.lb, sal.ub, by=inc)
press = seq(press.lb, press.ub, by=inc)

# Initiate arrays 
temp = c(1:length(tempC))
delV1 = c(1:length(tempC))
delV2 = c(1:length(tempC))
delVspc = c(1:length(tempC))
delVB = c(1:length(tempC))
delVw = c(1:length(tempC))

delk1 = c(1:length(tempC))
delk2 = c(1:length(tempC))
delkspc = c(1:length(tempC))
delkB = c(1:length(tempC))
delkw = c(1:length(tempC))

base2Darray = c(1:(length(tempC)*length(sal)))
dim(base2Darray) = c((length(tempC)), (length(sal)))

Ks1m_st = base2Darray
Ks2m_st =  base2Darray
logKsspcm_st = base2Darray 
Ksspcm_st = base2Darray
lnKsB_st = base2Darray
KsB_st = base2Darray
Ksw_st = base2Darray
K0a = base2Darray

base3Darray = c(1:(length(tempC)*length(sal)*length(press)))
dim(base3Darray) = c((length(tempC)), (length(sal)), (length(press)))      
               
Ks1a = base3Darray
Ks2a = base3Darray
Ksspca = base3Darray
KsBa = base3Darray
Kswa = base3Darray

# Constant (cm^3 bar mol^-1 K^-1)
R <- 83.131 

# Nested for loops to calculate 3D array for K*1a (i.e., array of potential K*1 values w/o any major ion effect)
for (i in 1:length(tempC)){
  for (j in 1:length(sal)){
    for (k in 1:length(press)){
      
      temp[i] <- tempC[i]+273.15
      
      Ks1m_st[i,j] <- exp(2.83655-2307.1266/temp[i]-1.5529413*(log(temp[i]))-((0.20760841+4.0484/temp[i])*sqrt(sal[j]))+0.0846834*sal[j]-0.00654208*(sal[j]^1.5)+log(1-(0.001005*sal[j])))
      Ks2m_st[i,j] <- exp(-9.226508-3351.6106/temp[i]-0.2005743*(log(temp[i]))-((0.106901773+23.9722/temp[i])*sqrt(sal[j]))+0.1130822*sal[j]-0.00846934*(sal[j]^1.5)+log(1-(0.001005*sal[j])))
      logKsspcm_st[i,j] <- ((-171.9065-0.077993*temp[i]+2839.319/temp[i]+71.595*(log(temp[i])/log(10))+(-0.77712+0.0028426*temp[i]+178.34/temp[i])*(sal[j]^0.5)-0.07711*sal[j]+0.0041249*(sal[j]^1.5)))
      Ksspcm_st[i,j] <- 10^(logKsspcm_st[i,j])
      lnKsB_st[i,j] <- ((-8966.9-2890.53*sal[j]^0.5-77.942*sal[j]+1.728*sal[j]^1.5-0.0996*sal[j]^2)/temp[i])+148.0248+137.1942*sal[j]^0.5+1.62142*sal[j]-(24.4344+25.085*sal[j]^0.5+0.2474*sal[j])*(log(temp[i]))+(0.053105*sal[j]^0.5*temp[i])
      KsB_st[i,j] <- exp(lnKsB_st[i,j])
      Ksw_st[i,j] <- exp(148.96502-13847.26/temp[i]-23.6521*(log(temp[i]))+(118.67/temp[i]-5.977+1.0495*(log(temp[i])))*(sal[j]^0.5)-0.01615*sal[j])
      K0a[i,j] <- exp(9345.17/temp[i]-60.2409+23.3585*(log(temp[i]/100))+sal[j]*(0.023517-0.00023656*temp[i]+0.0047036*((temp[i]/100)^2)))
      
      delV1[i] <- (-25.50)+0.1271*tempC[i]
      delV2[i] <- (-15.82)+(-0.0219*tempC[i])
      delVspc[i]<- (-48.76)+(0.5304*tempC[i])
      delVB[i] <- (-29.48)+0.1622*tempC[i]+(2.608/1000)*tempC[i]^2
      delVw[i] <- (-25.60)+0.2324*tempC[i]+(-3.6246/1000)*tempC[i]^2
      
      delk1[i] <- (-3.08/1000)+(0.0877/1000)*tempC[i]
      delk2[i] <- (1.13/1000)+(-0.1475/1000)*tempC[i]
      delkspc[i] <- (-11.76/1000)+(0.3692/1000)*tempC[i]
      delkB[i] <- -2.84/1000
      delkw[i] <- (-5.13/1000)+(0.0794/1000)*tempC[i]
      
      Ks1a[i,j,k] <- (exp(-((delV1[i]/(R*temp[i]))*press[k])+((0.5*delk1[i])/(R*temp[i]))*press[k]^2))*Ks1m_st[i,j]
      Ks2a[i,j,k] <- (exp(-((delV2[i]/(R*temp[i]))*press[k])+((0.5*delk2[i])/(R*temp[i]))*press[k]^2))*Ks2m_st[i,j]
      Ksspca[i,j,k] <- (exp(-((delVspc[i]/(R*temp[i]))*press[k])+((0.5*delkspc[i])/(R*temp[i]))*press[k]^2))*Ksspcm_st[i,j]
      KsBa[i,j,k] <- (exp(-((delVB[i]/(R*temp[i]))*press[k])+((0.5*delkB[i])/(R*temp[i]))*press[k]^2))*KsB_st[i,j]
      Kswa[i,j,k] <- (exp(-((delVw[i]/(R*temp[i]))*press[k])+((0.5*delkw[i])/(R*temp[i]))*press[k]^2))*Ksw_st[i,j]
}
}
}




# To look up a value in the PSM, pass K arrays, rnd.dig, and STP bounds to jags:
"tempC.lb" = tempC.lb, 
"tempC.ub" = tempC.ub,
"sal.lb" = sal.lb,
"sal.ub" = sal.ub,
"press.lb" = press.lb,
"press.ub" = press.ub,
"rnd.dig" = rnd.dig,
"K0a" = K0a,
"Ks1a" = Ks1a,
"Ks2a" = Ks2a,
"Ksspca" = Ksspca,
"KsBa" = KsBa,
"Kswa" = Kswa,

# then replace equil constant calcs in PSM time series for loop to: 
K0[i] <- K0a[(round((tempC[ai.prox[i]]-tempC.lb)/inc)+1), 
                (round((sal[ai.prox[i]]-sal.lb)/inc)+1)]

Ks1m[i] <- Ks1a[(round((tempC[ai.prox[i]]-tempC.lb)/inc)+1), 
                (round((sal[ai.prox[i]]-sal.lb)/inc)+1), 
                (round((press-press.lb)/inc)+1)]

Ks2m[i] <- Ks2a[(round((tempC[ai.prox[i]]-tempC.lb)/inc)+1), 
                (round((sal[ai.prox[i]]-sal.lb)/inc)+1), 
                (round((press-press.lb)/inc)+1)]

Ksspcm[i] <- Ksspca[(round((tempC[ai.prox[i]]-tempC.lb)/inc)+1), 
                    (round((sal[ai.prox[i]]-sal.lb)/inc)+1), 
                    (round((press-press.lb)/inc)+1)]

KsB[i] <- KsBa[(round((tempC[ai.prox[i]]-tempC.lb)/inc)+1), 
               (round((sal[ai.prox[i]]-sal.lb)/inc)+1), 
               (round((press-press.lb)/inc)+1)]

Ksw[i] <- Kswa[(round((tempC[ai.prox[i]]-tempC.lb)/inc)+1), 
               (round((sal[ai.prox[i]]-sal.lb)/inc)+1), 
               (round((press-press.lb)/inc)+1)]


# Scratch test calcs
tempC= 24.96521
sal = 35.032456

K0 <- K0a[(round((tempC-tempC.lb)/inc)+1), 
          (round((sal-sal.lb)/inc)+1)]
pK0 = -log(K0)/log(10)
print(pK0)




my_array <- array(rep(1:50,2), dim=c(10,5,5))
my_array2 <- c(1:(10*5*5))
dim(my_array2) <- c(10*5,5)
for (i in 1:5){
  my_array2 <- rbind(my_array[1:10,1:5,i])
  print(my_array2)
}
