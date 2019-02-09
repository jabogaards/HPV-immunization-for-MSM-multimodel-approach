vac.rate <- array(0,100)

vac.rate[15:19] <- 0.00905
vac.rate[20:24] <- 0.01980
vac.rate[25:29] <- 0.01899
vac.rate[30:34] <- 0.01708
vac.rate[35:39] <- 0.01643
vac.rate[40:44] <- 0.00711
vac.rate[45:49] <- 0.00711
vac.rate[50:54] <- 0.00471
vac.rate[55:59] <- 0.00471
vac.rate[60:64] <- 0.00300
vac.rate[65:69] <- 0.00300

plot(vac.rate, type="s", xlim=c(15,90), ylim=c(0,0.04),
 xlab="Age", ylab="Annual vaccination rate", main="MSM vaccine uptake")

vac.rate[70:100] <- 0.00609
lines(vac.rate, type="s", lty=2)

x <- 1:100
smooth <- dgamma((x-15)/5, shape=2, scale=2)

lines(0.02*smooth/max(smooth), col="firebrick1", lwd=2)
lines(0.04*smooth/max(smooth), col="firebrick1", lwd=1)
legend("topright", legend=c("HepB historic","HPV base-case","HPV alternative"),
 lty=1, lwd=c(1,2,1), col=c(1,2,2), bty="n")
