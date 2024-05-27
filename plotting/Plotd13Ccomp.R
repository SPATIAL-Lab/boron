d13Cbulk_1209 <- ShatskyLPEE_data[complete.cases(ShatskyLPEE_data$d18O), ]

d13Cbenth_1209
d13Cplank_1209

d13Cbenth_1262
d13Cbulk_1262





d18O <- d18O[1:length(d18O$age),1:6]
rownames(d18O) = NULL

Mvel <- d18O[d18O$species == "Grub", ] 
Asol <- d18O[d18O$species == "Tsac", ] 

plot(Mvel$age/1e3, Mvel$d18O, pch = 16, col = "coral", xlim = c(59,53), ylim = c(-0.5, -2.5), cex =0.4)
points(Asol$age/1e3, Asol$d18O, pch = 16, col = "cornflowerblue", cex =0.4)

read.csv(file = "Harperetal_resubm/out_primrary/inversion_sum.csv")