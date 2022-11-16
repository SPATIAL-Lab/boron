#
# This script takes output from the boron PSM driver, extracts and plots the data
# Dustin T. Harper
# August 9, 2022
#
#

# Extract data for plotting; construct parm data frames
c.multiout <- list()
m.multiout <- list()
pco2.multiout <- list()
dic.multiout <- list()
pH.multiout <- list()
sal.multiout <- list()
tempC.multiout <- list()
xca.multiout <- list()
xmg.multiout <- list()
xso4.multiout <- list()

for (i in 1:rows.id)
{
  pl.out <- parms.out[i]
  p.out <- unlist(pl.out)
  pm.out <- matrix(p.out, ncol = 9, byrow = FALSE)
  c.out <- pm.out[1,]
  m.out <- pm.out[7,]
  pco2.out <- pm.out[10,]
  dic.out <- pm.out[6,]
  pH.out <- pm.out[9,]
  sal.out <- pm.out[12,]
  tempC.out <- pm.out[13,]
  xca.out <- pm.out[14,]
  xmg.out <- pm.out[15,]
  xso4.out <- pm.out[16,]
  
  c.multiout[[i]] <-c.out
  m.multiout[[i]] <-m.out
  pco2.multiout[[i]] <-pco2.out
  dic.multiout[[i]] <-dic.out
  pH.multiout[[i]] <-pH.out
  sal.multiout[[i]] <-sal.out
  tempC.multiout[[i]] <-tempC.out
  xca.multiout[[i]] <-xca.out
  xmg.multiout[[i]] <-xmg.out
  xso4.multiout[[i]] <-xso4.out
}
c.df <- data.frame(matrix(unlist(c.multiout), nrow = rows.id, byrow = TRUE))
colnames <- c('mean', 'sd', 'dist2.5', 'dist25', 'dist50', 'dist75', 'dist97.5', 'Rhat','n.eff', 'age')
c.df$age <- age
colnames(c.df) <- colnames
m.df <- data.frame(matrix(unlist(m.multiout), nrow = rows.id, byrow = TRUE))
m.df$age <- age
colnames(m.df) <- colnames
pco2.df <- data.frame(matrix(unlist(pco2.multiout), nrow = rows.id, byrow = TRUE))
pco2.df$age <- age
colnames(pco2.df) <- colnames
dic.df <- data.frame(matrix(unlist(dic.multiout), nrow = rows.id, byrow = TRUE))
dic.df$age <- age
colnames(dic.df) <- colnames
pH.df <- data.frame(matrix(unlist(pH.multiout), nrow = rows.id, byrow = TRUE))
pH.df$age <- age
colnames(pH.df) <- colnames
sal.df <- data.frame(matrix(unlist(sal.multiout), nrow = rows.id, byrow = TRUE))
sal.df$age <- age
colnames(sal.df) <- colnames
tempC.df <- data.frame(matrix(unlist(tempC.multiout), nrow = rows.id, byrow = TRUE))
tempC.df$age <- age
colnames(tempC.df) <- colnames
xca.df <- data.frame(matrix(unlist(xca.multiout), nrow = rows.id, byrow = TRUE))
xca.df$age <- age
colnames(xca.df) <- colnames
xmg.df <- data.frame(matrix(unlist(xmg.multiout), nrow = rows.id, byrow = TRUE))
xmg.df$age <- age
colnames(xmg.df) <- colnames
xso4.df <- data.frame(matrix(unlist(xso4.multiout), nrow = rows.id, byrow = TRUE))
xso4.df$age <- age
colnames(xso4.df) <- colnames

# Plot parms of interest

ggplot() + 
  geom_pointrange(data = m.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(0,1.2) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "m") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = c.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(0,10) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "c") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = xca.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(10,20) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "[Ca]") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = xmg.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(30,40) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "[Mg]") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = xso4.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(10,22) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "[SO4]") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = tempC.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(20,40) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "temp (C)") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = sal.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(33,37) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "salinity (ppt)") +  
  theme_bw()

ggplot() + 
  geom_pointrange(data = pco2.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(0,0.002) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "pCO2 (atm)") +   
  theme_bw()

ggplot() + 
  geom_pointrange(data = dic.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(0.0015, 0.003) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "DIC (mol/kg)") +
  theme_bw()

ggplot() + 
  geom_pointrange(data = pH.df, mapping = aes(x=age, y=mean, ymin=dist25, ymax=dist75)) +
  ylim(7.4, 8.4) +
  xlim(55,59) +
  labs(x= "age (Ma)", y = "pH") +
  theme_bw()
