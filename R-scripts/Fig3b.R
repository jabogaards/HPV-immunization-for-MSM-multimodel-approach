library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### IMPROVED UPTAKE: UNPACK uptake.zip

timeOut <- 90
endTime <- 180
vac.start <- 100
time <- timeOut:endTime-vac.start

x11()
load("uptake.double.RData")
plot(dyn.uptake[[1]]$vacc ~ time,
 type="l", col=cols[2], lwd=2, xlim=c(0,50), ylim=c(0, 0.75),
 xlab="Years since start vaccination program", ylab="Proportion MSM vaccinated")

cross.age <- c(11.2,19.7,22.1)
for (i in 2:4) {
  lines(dyn.uptake[[i]]$vacc ~ time, lwd=2, col=cols[7+i])
  ref <- ceiling(cross.age[i-1])+10
  segments(x0=-10, x1=cross.age[i-1], y0=dyn.uptake[[i]]$vacc[ref], y1=dyn.uptake[[i]]$vacc[ref], lty=3, col=cols[7+i])
  segments(x0=cross.age[i-1], x1=cross.age[i-1], y0=-1, y1=dyn.uptake[[i]]$vacc[ref], lty=3, col=cols[7+i])
}

load("uptake.base.RData")
lines(dyn.uptake[[5]]$vacc ~ c(timeOut:endTime-vac.start), lwd=2, col=cols[4])

legendexp = c(expression("40% boys at 12y and all MSM"),
              expression("80% boys at 12y"),
              expression(paste("MSM", phantom(0)<=phantom(0), "26y (double uptake)", sep="")),
              expression(paste("MSM", phantom(0)<=phantom(0), "40y (double uptake)", sep="")),
              expression("all MSM (double uptake)"))
legend("topleft", legend=legendexp, lty=1, lwd=2,
 col=cols[c(4,2,9,10,11)], title="Vaccine eligibility:", title.adj=0.1, bty="n")
title("Projected HPV vaccine uptake, improved")
