library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### CONSERVATIVE SCENARIOS: UNPACK sens.leakyifss.zip

timeOut <- 90
endTime <- 180
vac.start <- 100
time <- timeOut:endTime-vac.start

load("delta.AIC.RData")
aic.weights <- exp(-0.5*delta.aic)/sum(exp(-0.5*delta.aic))

##### MODELS
setwd("msm.leakyifss")
models <- c("SIS","SISPS","SISPminRS","SISPmaxRS","SIRS","SIS10RS","SIS33RS","SIRlocalS","SIS10RlocalS","SIS33RlocalS",
		"SIRanalS","SIS10RanalS","SIS33RanalS","SIRpenileS","SIS10RpenileS","SIS33RpenileS","SIL","SIS10L","SIS33L","SIR33L")


windows(width = 12, height = 8)
layout(matrix(1:8, 2, 4, byrow = TRUE), respect = F)

### PENILE
fill.y1_ <- function() {
  y.model <- array(NA,dim=c(360,length(time)))
  for (listindex in 1:18) {
    y.model[0*18+listindex,] <- dyn.SIS[[listindex]]$mod16_x10+dyn.SIS[[listindex]]$mod16_x11
    y.model[1*18+listindex,] <- dyn.SISPS[[listindex]]$mod16_x10+dyn.SISPS[[listindex]]$mod16_x11
    y.model[2*18+listindex,] <- dyn.SISPminRS[[listindex]]$mod16_x10+dyn.SISPminRS[[listindex]]$mod16_x11
    y.model[3*18+listindex,] <- dyn.SISPmaxRS[[listindex]]$mod16_x10+dyn.SISPmaxRS[[listindex]]$mod16_x11
    y.model[4*18+listindex,] <- dyn.SIRS[[listindex]]$mod16_x10+dyn.SIRS[[listindex]]$mod16_x11
    y.model[5*18+listindex,] <- dyn.SIS10RS[[listindex]]$mod16_x10+dyn.SIS10RS[[listindex]]$mod16_x11
    y.model[6*18+listindex,] <- dyn.SIS33RS[[listindex]]$mod16_x10+dyn.SIS33RS[[listindex]]$mod16_x11
    y.model[7*18+listindex,] <- dyn.SIRlocalS[[listindex]]$mod16_x10+dyn.SIRlocalS[[listindex]]$mod16_x11
    y.model[8*18+listindex,] <- dyn.SIS10RlocalS[[listindex]]$mod16_x10+dyn.SIS10RlocalS[[listindex]]$mod16_x11
    y.model[9*18+listindex,] <- dyn.SIS33RlocalS[[listindex]]$mod16_x10+dyn.SIS33RlocalS[[listindex]]$mod16_x11
    y.model[10*18+listindex,] <- dyn.SIRanalS[[listindex]]$mod16_x10+dyn.SIRanalS[[listindex]]$mod16_x11
    y.model[11*18+listindex,] <- dyn.SIS10RanalS[[listindex]]$mod16_x10+dyn.SIS10RanalS[[listindex]]$mod16_x11
    y.model[12*18+listindex,] <- dyn.SIS33RanalS[[listindex]]$mod16_x10+dyn.SIS33RanalS[[listindex]]$mod16_x11
    y.model[13*18+listindex,] <- dyn.SIRpenileS[[listindex]]$mod16_x10+dyn.SIRpenileS[[listindex]]$mod16_x11
    y.model[14*18+listindex,] <- dyn.SIS10RpenileS[[listindex]]$mod16_x10+dyn.SIS10RpenileS[[listindex]]$mod16_x11
    y.model[15*18+listindex,] <- dyn.SIS33RpenileS[[listindex]]$mod16_x10+dyn.SIS33RpenileS[[listindex]]$mod16_x11
    y.model[16*18+listindex,] <- dyn.SIL[[listindex]]$mod16_x10+dyn.SIL[[listindex]]$mod16_x11
    y.model[17*18+listindex,] <- dyn.SIS10L[[listindex]]$mod16_x10+dyn.SIS10L[[listindex]]$mod16_x11
    y.model[18*18+listindex,] <- dyn.SIS33L[[listindex]]$mod16_x10+dyn.SIS33L[[listindex]]$mod16_x11
    y.model[19*18+listindex,] <- dyn.SIR33L[[listindex]]$mod16_x10+dyn.SIR33L[[listindex]]$mod16_x11
  }
  return(y.model)
}

for(i in 1:20) load(paste("dyn.",models[i],".b4age27.RData",sep=""))
y.model11 <- fill.y1_()
par(mfg=c(1,1))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.05), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, penile")
for (i in 1:360) lines(y.model11[i,] ~ time, col="grey")
lines(apply(y.model11 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model11, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model11, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("Eligibility: MSM <= 26y")

for(i in 1:20) load(paste("dyn.",models[i],".b4age40.RData",sep=""))
y.model12 <- fill.y1_()
par(mfg=c(1,2))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.05), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, penile")
for (i in 1:360) lines(y.model12[i,] ~ time, col="grey")
lines(apply(y.model12 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model12, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model12, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("Eligibility: MSM <= 40y")

for(i in 1:20) load(paste("dyn.",models[i],".b4age80.RData",sep=""))
y.model13 <- fill.y1_()
par(mfg=c(1,3))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.05), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, penile")
for (i in 1:360) lines(y.model13[i,] ~ time, col="grey")
lines(apply(y.model13 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model13, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model13, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("Eligibility: all MSM")

### ANAL
fill.y_1 <- function() {
  y.model <- array(NA,dim=c(360,length(time)))
  for (listindex in 1:18) {
    y.model[0*18+listindex,] <- dyn.SIS[[listindex]]$mod16_x01+dyn.SIS[[listindex]]$mod16_x11
    y.model[1*18+listindex,] <- dyn.SISPS[[listindex]]$mod16_x01+dyn.SISPS[[listindex]]$mod16_x11
    y.model[2*18+listindex,] <- dyn.SISPminRS[[listindex]]$mod16_x01+dyn.SISPminRS[[listindex]]$mod16_x11
    y.model[3*18+listindex,] <- dyn.SISPmaxRS[[listindex]]$mod16_x01+dyn.SISPmaxRS[[listindex]]$mod16_x11
    y.model[4*18+listindex,] <- dyn.SIRS[[listindex]]$mod16_x01+dyn.SIRS[[listindex]]$mod16_x11
    y.model[5*18+listindex,] <- dyn.SIS10RS[[listindex]]$mod16_x01+dyn.SIS10RS[[listindex]]$mod16_x11
    y.model[6*18+listindex,] <- dyn.SIS33RS[[listindex]]$mod16_x01+dyn.SIS33RS[[listindex]]$mod16_x11
    y.model[7*18+listindex,] <- dyn.SIRlocalS[[listindex]]$mod16_x01+dyn.SIRlocalS[[listindex]]$mod16_x11
    y.model[8*18+listindex,] <- dyn.SIS10RlocalS[[listindex]]$mod16_x01+dyn.SIS10RlocalS[[listindex]]$mod16_x11
    y.model[9*18+listindex,] <- dyn.SIS33RlocalS[[listindex]]$mod16_x01+dyn.SIS33RlocalS[[listindex]]$mod16_x11
    y.model[10*18+listindex,] <- dyn.SIRanalS[[listindex]]$mod16_x01+dyn.SIRanalS[[listindex]]$mod16_x11
    y.model[11*18+listindex,] <- dyn.SIS10RanalS[[listindex]]$mod16_x01+dyn.SIS10RanalS[[listindex]]$mod16_x11
    y.model[12*18+listindex,] <- dyn.SIS33RanalS[[listindex]]$mod16_x01+dyn.SIS33RanalS[[listindex]]$mod16_x11
    y.model[13*18+listindex,] <- dyn.SIRpenileS[[listindex]]$mod16_x01+dyn.SIRpenileS[[listindex]]$mod16_x11
    y.model[14*18+listindex,] <- dyn.SIS10RpenileS[[listindex]]$mod16_x01+dyn.SIS10RpenileS[[listindex]]$mod16_x11
    y.model[15*18+listindex,] <- dyn.SIS33RpenileS[[listindex]]$mod16_x01+dyn.SIS33RpenileS[[listindex]]$mod16_x11
    y.model[16*18+listindex,] <- dyn.SIL[[listindex]]$mod16_x01+dyn.SIL[[listindex]]$mod16_x11
    y.model[17*18+listindex,] <- dyn.SIS10L[[listindex]]$mod16_x01+dyn.SIS10L[[listindex]]$mod16_x11
    y.model[18*18+listindex,] <- dyn.SIS33L[[listindex]]$mod16_x01+dyn.SIS33L[[listindex]]$mod16_x11
    y.model[19*18+listindex,] <- dyn.SIR33L[[listindex]]$mod16_x01+dyn.SIR33L[[listindex]]$mod16_x11
  }
  return(y.model)
}

for(i in 1:20) load(paste("dyn.",models[i],".b4age27.RData",sep=""))
y.model21 <- fill.y_1()
par(mfg=c(2,1))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.15), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, anal")
for (i in 1:360) lines(y.model21[i,] ~ time, col="grey")
lines(apply(y.model21 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model21, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model21, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("Eligibility: MSM <= 26y")

for(i in 1:20) load(paste("dyn.",models[i],".b4age40.RData",sep=""))
y.model22 <- fill.y_1()
par(mfg=c(2,2))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.15), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, anal")
for (i in 1:360) lines(y.model22[i,] ~ time, col="grey")
lines(apply(y.model22 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model22, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model22, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("Eligibility: MSM <= 40y")

for(i in 1:20) load(paste("dyn.",models[i],".b4age80.RData",sep=""))
y.model23 <- fill.y_1()
par(mfg=c(2,3))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.15), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, anal")
for (i in 1:360) lines(y.model23[i,] ~ time, col="grey")
lines(apply(y.model23 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model23, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model23, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("Eligibility: all MSM")

### BENCHMARK: 40% PREADOLESCENT BOYS
setwd("../boys.leaky")

for(i in 1:20) load(paste("dyn.",models[i],".age12.RData",sep=""))
y.model14 <- fill.y1_()
par(mfg=c(1,4))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.05), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, penile")
for (i in 1:360) lines(y.model14[i,] ~ time, col="grey")
lines(apply(y.model14 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model14, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model14, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("40% boys at 12y")

for(i in 1:20) load(paste("dyn.",models[i],".age12.RData",sep=""))
y.model24 <- fill.y_1()
par(mfg=c(2,4))
plot.new();
plot(x=c(-5, 75), y=c(0, 0.15), type="n", xlab="Years since start vaccination program", ylab="HPV16 prevalence, anal")
for (i in 1:360) lines(y.model24[i,] ~ time, col="grey")
lines(apply(y.model24 * as.vector(aic.weights), 2, sum) ~ time, col=2)
lines(apply(y.model24, 2, function(x) {quantile(x,0.05)}) ~ time, col=2, lty=3)
lines(apply(y.model24, 2, function(x) {quantile(x,0.95)}) ~ time, col=2, lty=3)
box()
title("40% boys at 12y")
