library("foreign")
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### DATA
baseline <- read.csv("baseline.csv",header=T)
select <- !is.na(baseline$hpv16_anal) &  !is.na(baseline$hpv16_penile)

x_a <- round(baseline$age)[select]
agecat <- list()
Age_breaks <- matrix(c(21:30,31:40,41:50,51:60,61:70), byrow=T, nrow=5)
Age_midPoints <- list()

for (j in 1:dim(Age_breaks)[2]) {
  agecat[[j]] <- cut(floor(baseline$age[select]), breaks=c(1,Age_breaks[,j],99))
  print(table(agecat[[j]]))
  Age_midPoints[[j]] <- array(NA, length(table(agecat[[j]])))
  for (i in 1:length(table(agecat[[j]])))
    Age_midPoints[[j]][i] <- mean(baseline$age[select][agecat[[j]]==levels(agecat[[j]])[i]])
}

# penile HPV16 infection amongst HIV-negative MSM
prev16_penile <- array(NA, dim=c(dim(Age_breaks)[2],length(levels(agecat[[j]]))))
prev16_penile_lcl <- array(NA, dim(prev16_penile))
prev16_penile_ucl <- array(NA, dim(prev16_penile))
for (j in 1:dim(Age_breaks)[2]) {
  for (i in 1:length(levels(agecat[[j]]))) {
    mask <- agecat[[j]] == levels(agecat[[j]])[i]
    print(sum(mask))
    success <- sum(baseline$hpv16_penile[select][mask])
    prev16_penile[j,i] <- success / sum(mask)
    prev16_penile_lcl[j,i] <- binom.test(success, sum(mask), conf.level=0.95)$conf.int[1]
    prev16_penile_ucl[j,i] <- binom.test(success, sum(mask), conf.level=0.95)$conf.int[2]
  }
}

# anal HPV16 infection amongst HIV-negative MSM
prev16_anal <- array(NA, dim=c(dim(Age_breaks)[2],length(levels(agecat[[j]]))))
prev16_anal_lcl <- array(NA, dim(prev16_anal))
prev16_anal_ucl <- array(NA, dim(prev16_anal))
for (j in 1:dim(Age_breaks)[2]) {
  for (i in 1:length(levels(agecat[[j]]))) {
    mask <- agecat[[j]] == levels(agecat[[j]])[i]
    print(sum(mask))
    success <-  sum(baseline$hpv16_anal[select][mask])
    prev16_anal[j,i] <- success / sum(mask)
    prev16_anal_lcl[j,i] <- binom.test(success, sum(mask), conf.level=0.95)$conf.int[1]
    prev16_anal_ucl[j,i] <- binom.test(success, sum(mask), conf.level=0.95)$conf.int[2]
  }
}

##### MODELS
minAge <- 9
maxAge <- 80
gridAge <- 1/4
axisAge <- (minAge/gridAge):(maxAge/gridAge)*gridAge

##### PENILE PREVALENCE
x11()
load("base.penile.RData")
plot(penile.prev[1,] ~ axisAge, type="n", xlim=c(18,70), ylim=c(0,0.2),
  xlab="Age", ylab="Proportion HPV16 DNA positive")
for (j in 1:dim(Age_breaks)[2]) {
  for (i in 1:length(levels(agecat[[j]]))) {
    segments(x0=Age_midPoints[[j]][i], y0=prev16_penile_lcl[j,i],
             x1=Age_midPoints[[j]][i], y1=prev16_penile_ucl[j,i], col="grey")
  }
}

for (i in 1:360) lines(penile.prev[i,] ~ axisAge, col=cols[8])
# all models with latency
for (i in 289:360) lines(penile.prev[i,] ~ axisAge, col=cols[9])

for (j in 1:dim(Age_breaks)[2]) {
  lines(prev16_penile[j,] ~ Age_midPoints[[j]], type="p", pch=16, cex=1.2)
  lines(prev16_penile[j,] ~ Age_midPoints[[j]], type="p", pch=20, cex=1.4, col="white")
  lines(prev16_penile[j,] ~ Age_midPoints[[j]], type="p", pch=16, cex=0.7, col="grey")
}
box()
title("Penile HPV16 prevalence")

##### ANAL PREVALENCE
x11()
load("base.anal.RData")
plot(anal.prev[1,] ~ axisAge, type="n", xlim=c(18,70), ylim=c(0,0.3),
  xlab="Age", ylab="Proportion HPV16 DNA positive")
#lines(base.SIL[[7]]$mod16_x01 ~ axisAge)
for (j in 1:dim(Age_breaks)[2]) {
  for (i in 1:length(levels(agecat[[j]]))) {
    segments(x0=Age_midPoints[[j]][i], y0=prev16_anal_lcl[j,i],
             x1=Age_midPoints[[j]][i], y1=prev16_anal_ucl[j,i], col="grey")
  }
}

for (i in 1:360) lines(anal.prev[i,] ~ axisAge, col=cols[8])
# all models with latency
for (i in 289:360) lines(anal.prev[i,] ~ axisAge, col=cols[9])

for (j in 1:dim(Age_breaks)[2]) {
  lines(prev16_anal[j,] ~ Age_midPoints[[j]], type="p", pch=16, cex=1.2)
  lines(prev16_anal[j,] ~ Age_midPoints[[j]], type="p", pch=20, cex=1.4, col="white")
  lines(prev16_anal[j,] ~ Age_midPoints[[j]], type="p", pch=16, cex=0.7, col="grey")
}
box()
title("Anal HPV16 prevalence")
