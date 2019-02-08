library("foreign")

sexpartners <- read.csv("sexpartners.csv",header=T)

hist(sexpartners$age)
summary(sexpartners$age)
agecat <- ifelse(sexpartners$age < 25, 1, NA)
agecat[sexpartners$age >= 25 & sexpartners$age < 35 ] <- 2
agecat[sexpartners$age >= 35 & sexpartners$age < 45 ] <- 3
agecat[sexpartners$age >= 45 & sexpartners$age < 55 ] <- 4
agecat[sexpartners$age >= 55 & sexpartners$age < 65 ] <- 5
agecat[sexpartners$age >= 65] <- 6
table(sexpartners$hivstat, agecat)
table(agecat)

# activity all participants
boxplot(log(sexpartners$nranpa6m+1) ~ agecat, 
 names=c("<25","25-35","35-45","45-55","55-65","65+"),
 xlab="Age", ylab="Anal sex partners last 6 months, logarithmic")
title("All H2M participants")

# activity HIV-negative participants
boxplot(log(sexpartners$nranpa6m[sexpartners$hivstat==0]+1) ~ agecat[sexpartners$hivstat==0], 
 names=c("<25","25-35","35-45","45-55","55-65","65+"),
 xlab="Age", ylab="Anal sex partners last 6 months, logarithmic")
title("HIV-negative participants")

plot(nranpa6m ~ age, data=sexpartners[!sexpartners$hivstat,])
title("HIV-negative participants")

mean_all <- array(NA, length(table(agecat)))
mean_hivneg <- array(NA, length(table(agecat)))
sd_hivneg <- array(NA, length(table(agecat)))
cv_hivneg <- array(NA, length(table(agecat)))

for (i in 1:length(table(agecat))) {
  mean_all[i] <- mean(sexpartners$nranpa6m[agecat==i], na.rm=T)
  mean_hivneg[i] <- mean(sexpartners$nranpa6m[sexpartners$hivstat==0 & agecat==i], na.rm=T)
  sd_hivneg[i] <- sd(sexpartners$nranpa6m[sexpartners$hivstat==0 & agecat==i], na.rm=T)
  cv_hivneg[i] <- sd_hivneg[i]/mean_hivneg[i]
}

plot(mean_hivneg, type="l", ylim=c(0,15), xaxt="n")
lines(sd_hivneg, col="red")
lines(cv_hivneg, col="red", lty=2)
axis(1, at=1:6, c("<25","25-35","35-45","45-55","55-65","65+"))
title("HIV-negative participants")

# estimate mean contact rate plus cv assuming simple parabolic function...
select <- !is.na(sexpartners$nranpa6m) & !sexpartners$hivstat
maxAge <- 80
myfun <- function(c_pars) {
  peak <- c_pars[1]
  cv <- c_pars[2]
  expect <- 4 * peak * sexpartners$age[select]/maxAge * (1 - sexpartners$age[select]/maxAge)
  logLiki <- log(dnorm(sexpartners$nranpa6m[select], expect, cv*expect))
  return(-sum(logLiki))
}
optim(c(3,1),myfun)

### Ergo: E(y|age40) = 4.636 with cv = 1.862
c_peak <- 4.64
x <- c(15,20,25,30,35,40,45,50,55,60,65,70,75)
simanpa6m <- 4 * c_peak * x/maxAge * (1 - x/maxAge)
lines(simanpa6m ~ c(0.5*(1:length(x))), col="magenta")
abline(h=1.86, col="magenta", lty=2)

plot(mean_hivneg, type="b", pch=19, ylim=c(0,11), xaxt="n",
 ylab="Mean no. anal sex partners/6m", xlab="Age")
for (i in 1:length(table(agecat))) {
  segments(x0=i, y0=mean_hivneg[i]-1.96*sd_hivneg[i]/sqrt(table(sexpartners$hivstat, agecat)[1,i]),
           x1=i, y1=mean_hivneg[i]+1.96*sd_hivneg[i]/sqrt(table(sexpartners$hivstat, agecat)[1,i]),
  lty=2)
}
axis(1, at=1:6, c("<25","25-35","35-45","45-55","55-65","65+"))
lines(simanpa6m ~ c(0.5*(1:length(x))), col="firebrick1", lwd=2)
lines(mean_all, type="b", lty=3, pch=1)

# ... what if we'd leave 1% extremes out
discard <- floor(sum(select)/100)
indices <- order(sexpartners$nranpa6m[select])[1:(sum(select)-discard)]
keep <- sexpartners$dummyid[select][indices]

for (i in 1:length(table(agecat))) {
  mean_hivneg[i] <- mean(sexpartners$nranpa6m[sexpartners$dummyid %in% keep & agecat==i])
  sd_hivneg[i] <- sd(sexpartners$nranpa6m[sexpartners$dummyid %in% keep & agecat==i])
}
lines(mean_hivneg, type="b")

legend("topleft", lty=c(3,1,1), pch=c(1,19,1), bty="n",
 legend=c("All H2M participants","HIV-negatives w/ 95% CI","HIV-negatives w/o 1% extremes"))
title("Partner acquisition")
