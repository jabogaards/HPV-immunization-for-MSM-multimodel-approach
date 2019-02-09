library(interval)
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

clear <- read.csv("clear.anal.sens.csv",header=T)
clear$T_UPPER[is.na(clear$T_UPPER)] <- Inf

fit <- icfit(Surv(clear$T_LOWER,clear$T_UPPER,type="interval2")~1, conf.int=T, data=clear)
summary(fit)

plot(fit)
title("Anal HPV16 clearance")

x=0:1000

gamma1 <- 0.99/365
y1 <- exp(-gamma1*x)
lines(y1~x, lwd=2, col=cols[9])

shape <- 0.83
gamma2 <- 1.03/365
y2 <- exp(-(gamma2*x)^shape)
lines(y2~x, lwd=2, col=cols[11])

gamma <- 1.7854680
zeta <- 3.1086502
xi <- 0.7450531
y3 <- (gamma-xi)/(gamma+zeta-xi)*exp(-(gamma+zeta)*x/365)+zeta/(gamma+zeta-xi)*exp(-xi*x/365)
lines(y3~x, lwd=2, col=cols[4])

frac <- 0.1
psi <- 0.7633826
latent <- exp(1.2619900)/(1+exp(1.2619900))
y4 <- (1-latent)*exp(-psi/frac*x/365) + latent*exp(-psi*x/365)
lines(y4~x, lwd=2, col=cols[3])

frac <- 0.33
psi <- 0.7353683
latent <- exp(0.6530062)/(1+exp(0.6530062))
y5 <- (1-latent)*exp(-psi/frac*x/365) + latent*exp(-psi*x/365)
lines(y5~x, lwd=2, col=cols[1])

legend("topright", legend=c("Exponential","Weibull","Biphasic", "10% Latent mixture", "33% Latent mixture"),
 lty=1, lwd=2, col=cols[c(9,11,4,3,1)], title="Survival function", bty="n")
box()
