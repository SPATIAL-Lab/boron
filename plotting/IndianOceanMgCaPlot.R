ind_mg <- read.csv(file = "Harperetal_resubm/data/IndianOceanMgCa.csv")


plot(ind_mg$age, ind_mg$Ms_MgCa, xlim = c(59000, 53000), ylim = c(3,6), pch =15, cex =0.5)
points(ind_mg$age, ind_mg$Mv_MgCa, pch=15, cex = 0.5)
arrows(ind_mg$age, ind_mg$Ms_MgCa*0.97, ind_mg$age, ind_mg$Ms_MgCa*1.03, length=0.015, angle=90, code=3)
arrows(ind_mg$age, ind_mg$Mv_MgCa*0.97, ind_mg$age, ind_mg$Mv_MgCa*1.03, length=0.015, angle=90, code=3)


