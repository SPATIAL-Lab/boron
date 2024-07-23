
d13Cbulk_1209 <- read_xlsx(path = "plotting/data/LPEE1209_1262.xlsx", sheet = "1209bulk")
d13Cbenth_1209 <- read_xlsx(path = "plotting/data/LPEE1209_1262.xlsx", sheet = "1209benth")
d13Cplank_1209 <- read_xlsx(path = "plotting/data/LPEE1209_1262.xlsx", sheet = "1209plank")

d13Cbenth_1262 <- read_xlsx(path = "plotting/data/LPEE1209_1262.xlsx", sheet = "1262benth")
d13Cbulk_1262 <- read_xlsx(path = "plotting/data/LPEE1209_1262.xlsx", sheet = "1262bulk")

plot(d13Cplank_1209$age, d13Cplank_1209$d13C,  col = rgb(0,0,0, maxColorValue = 255), xlim = c(59,53), 
     ylim = c(-2.5, 6), xlab = "Age (Ma)", type = "p", cex =0.4)
lines(d13Cbulk_1262$age, d13Cbulk_1262$d13C,  col = rgb(230,159,0, maxColorValue = 255))
lines(d13Cbulk_1209$age, d13Cbulk_1209$d13C,  col = rgb(86,180,233, maxColorValue = 255))
lines(d13Cbenth_1262$age, d13Cbenth_1262$d13C,  col = rgb(213,94,0, maxColorValue = 255))
lines(d13Cbenth_1209$age, d13Cbenth_1209$d13C,  col = rgb(0,114,178, maxColorValue = 255))
lines(d13Cbenth_1262$age, d13Cbenth_1262$d18O-1,  col = rgb(204,121,167, maxColorValue = 255))
lines(d13Cbenth_1209$age, d13Cbenth_1209$d18O-1,  col = rgb(0,158,115, maxColorValue = 255))
legend("topright", legend=c(expression(paste("Site 1209 planktic ", delta^{13}, "C")),
                            expression(paste("Site 1262 bulk ", delta^{13}, "C")), 
                            expression(paste("Site 1209 bulk ", delta^{13}, "C")), 
                            expression(paste("Site 1262 benthic ", delta^{13}, "C")), 
                            expression(paste("Site 1209 benthic ", delta^{13}, "C")), 
                            expression(paste("Site 1262 benthic ", delta^{18}, "O")),
                            expression(paste("Site 1209 benthic ", delta^{18}, "O"))),  
       fill = c(rgb(0,0,0, maxColorValue = 255),rgb(230,159,0, maxColorValue = 255),
                rgb(86,180,233, maxColorValue = 255), rgb(213,94,0, maxColorValue = 255),
                rgb(0,114,178, maxColorValue = 255),rgb(204,121,167, maxColorValue = 255),
                rgb(0,158,115, maxColorValue = 255)))

