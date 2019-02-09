library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### H2M Data
sexage <- read.csv("sexage.csv",header=T)

##### POPULATION AGE DISTRIBUTION: USE sisprs.output (in sisprs.zip)
load("sisprs.output.RData")
N_age <- rowSums(output$N[dim(output$N)[1],,])

minAge <- 9
maxAge <- 80
gridAge <- 1/4
axisAge <- (minAge/gridAge):(maxAge/gridAge)*gridAge

sum(N_age*axisAge) # mean age in model
sum(N_age[axisAge<=36.75]) # median in model

age.cdf <- array(NA, dim=maxAge)
age.ecdf <- array(NA, dim=maxAge)
for (i in 1:maxAge) {
  age.cdf[i] <- sum(N_age[axisAge<=i])
  age.ecdf[i] <- sum(sexage$age<=i)/dim(sexage)[1]
}
plot(age.ecdf, type="s")
lines(age.cdf, col="red")

# pdf more informative than cdf...
age.pdf <- diff(age.cdf, 1)
age.epdf <- diff(age.ecdf, 1)

plot(age.epdf, type="s", xlab="Age", ylab="Density")
lines(age.pdf, col="red")
legend("topleft", lty=1, col=c(1,2), bty="n",
 legend=c("H2M participants","Modelled population"))
title("Age distribution")
