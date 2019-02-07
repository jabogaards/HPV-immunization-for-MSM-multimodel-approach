library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### BASE-CASE UPTAKE
timeOut <- 90
endTime <- 180
vac.start <- 100

x11()
load("uptake.base.RData")
plot(dyn.uptake[[1]]$vacc ~ c(timeOut:endTime-(vac.start)),
 type="l", col=cols[2], lwd=2, xlim=c(0,50), ylim=c(0, 0.5),
 xlab="Years since start vaccination program", ylab="Proportion MSM vaccinated")

cross.age <- c(11.5,21.5,24)
for (i in 2:4) {
  lines(dyn.uptake[[i]]$vacc ~ c(timeOut:endTime-vac.start), lwd=2, col=cols[7+i])
  ref <- sum(dyn.uptake[[1]]$vacc<=dyn.uptake[[i]]$vacc)
  segments(x0=-10, x1=cross.age[i-1], y0=dyn.uptake[[i]]$vacc[ref+2], y1=dyn.uptake[[i]]$vacc[ref+2], lty=3, col=cols[7+i])
  segments(x0=cross.age[i-1], x1=cross.age[i-1], y0=-1, y1=dyn.uptake[[i]]$vacc[ref+2], lty=3, col=cols[7+i])
}

lines(dyn.uptake[[5]]$vacc ~ c(timeOut:endTime-vac.start), lwd=2, col=cols[4])

legendexp = c(expression("40% boys at 12y and all MSM"),
              expression("40% boys at 12y"),
              expression(paste("MSM", phantom(0)<=phantom(0), "26y", sep="")),
              expression(paste("MSM", phantom(0)<=phantom(0), "40y", sep="")),
              expression("all MSM"))
legend("topleft", legend=legendexp, lty=1, lwd=2,
 col=cols[c(4,2,9,10,11)], title="Vaccine eligibility:", title.adj=0.1, bty="n")
title("Projected HPV vaccine uptake, base-case")
