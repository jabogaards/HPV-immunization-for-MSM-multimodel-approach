library(interval)
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

clear <- read.csv("clear.penile.csv",header=T)
clear$T_UPPER[is.na(clear$T_UPPER)] <- Inf

fit <- icfit(Surv(clear$T_LOWER,clear$T_UPPER,type="interval2")~1, conf.int=T, data=clear)
summary(fit)

plot(fit, XLAB="", YLAB="")
title("Penile HPV16 clearance", cex.main=1.5)

x=0:1000

gamma1 <- exp(-6.034)
y1 <- exp(-gamma1*x)
lines(y1~x, lwd=2, col=cols[9])

shape <- 0.6796
gamma2 <- exp(-6.0922/1.4715)
y2 <- exp(-gamma2*x^shape)
lines(y2~x, lwd=2, col=cols[11])

gamma <- 3.4009980
zeta <- 5.4783869
xi <- 0.4250297
y3 <- (gamma-xi)/(gamma+zeta-xi)*exp(-(gamma+zeta)*x/365)+zeta/(gamma+zeta-xi)*exp(-xi*x/365)
lines(y3~x, lwd=2, col=cols[4])

frac <- 0.1
psi <- 0.5070302
latent <- exp(0.7009352)/(1+exp(0.7009352))
y4 <- (1-latent)*exp(-psi/frac*x/365) + latent*exp(-psi*x/365)
lines(y4~x, lwd=2, col=cols[3])

frac <- 0.33
psi <- 0.6118052
latent <- exp(0.5627178)/(1+exp(0.5627178))
y5 <- (1-latent)*exp(-psi/frac*x/365) + latent*exp(-psi*x/365)
lines(y5~x, lwd=2, col=cols[1])

legend("topright", legend=c("Exponential","Weibull","Biphasic", "10% Latent mixture", "33% Latent mixture"),
 lty=1, lwd=2, col=cols[c(9,11,4,3,1)], title="Survival function", cex=1.25, bty="n")
box()
