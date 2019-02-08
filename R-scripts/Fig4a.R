library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### BASE-CASE SCENARIOS: UNPACK base.case.zip

timeOut <- 90
endTime <- 180
vac.start <- 100
time <- timeOut:endTime-vac.start

load("delta.AIC.RData")
aic.weights <- exp(-0.5*delta.aic)/sum(exp(-0.5*delta.aic))

##### MODELS
setwd("msm")
models <- c("SIS","SISPS","SISPminRS","SISPmaxRS","SIRS","SIS10RS","SIS33RS","SIRlocalS","SIS10RlocalS","SIS33RlocalS",
		"SIRanalS","SIS10RanalS","SIS33RanalS","SIRpenileS","SIS10RpenileS","SIS33RpenileS","SIL","SIS10L","SIS33L","SIR33L")


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

for(i in 1:20) load(paste("dyn.",models[i],".b4age40.RData",sep=""))
y.model12 <- fill.y1_()

for(i in 1:20) load(paste("dyn.",models[i],".b4age80.RData",sep=""))
y.model13 <- fill.y1_()

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

for(i in 1:20) load(paste("dyn.",models[i],".b4age40.RData",sep=""))
y.model22 <- fill.y_1()

for(i in 1:20) load(paste("dyn.",models[i],".b4age80.RData",sep=""))
y.model23 <- fill.y_1()

##### BENCHMARK: 40% PREADOLESCENT BOYS
setwd("../boys")

for(i in 1:20) load(paste("dyn.",models[i],".age12.RData",sep=""))
y.model14 <- fill.y1_()
y.model24 <- fill.y_1()

##### COMBINED STRATEGY
setwd("../combi")

for(i in 1:20) load(paste("dyn.",models[i],".b4age80.RData",sep=""))
y.model15 <- fill.y1_()
y.model25 <- fill.y_1()

plot(x=c(0, 50), y=c(0, 0.16), type="n", xlab="Years since start vaccination program", ylab="MSM population prevalence")
lines(apply(y.model24, 2, function(x) {quantile(x,0.05)}) ~ time, lty=3, col=cols[2])
lines(apply(y.model24, 2, function(x) {quantile(x,0.95)}) ~ time, lty=3, col=cols[2])
lines(apply(y.model21, 2, function(x) {quantile(x,0.05)}) ~ time, lty=3, col=cols[9])
lines(apply(y.model21, 2, function(x) {quantile(x,0.95)}) ~ time, lty=3, col=cols[9])
lines(apply(y.model22, 2, function(x) {quantile(x,0.05)}) ~ time, lty=3, col=cols[10])
lines(apply(y.model22, 2, function(x) {quantile(x,0.95)}) ~ time, lty=3, col=cols[10])
lines(apply(y.model23, 2, function(x) {quantile(x,0.05)}) ~ time, lty=3, col=cols[11])
lines(apply(y.model23, 2, function(x) {quantile(x,0.95)}) ~ time, lty=3, col=cols[11])
lines(apply(y.model25, 2, function(x) {quantile(x,0.05)}) ~ time, lty=3, col=cols[4])
lines(apply(y.model25, 2, function(x) {quantile(x,0.95)}) ~ time, lty=3, col=cols[4])

lines(apply(y.model24 * as.vector(aic.weights), 2, sum) ~ time, lwd=2, col=cols[2])
lines(apply(y.model21 * as.vector(aic.weights), 2, sum) ~ time, lwd=2, col=cols[9])
lines(apply(y.model22 * as.vector(aic.weights), 2, sum) ~ time, lwd=2, col=cols[10])
lines(apply(y.model23 * as.vector(aic.weights), 2, sum) ~ time, lwd=2, col=cols[11])
lines(apply(y.model25 * as.vector(aic.weights), 2, sum) ~ time, lwd=2, col=cols[4])

legendexp = c(expression("40% boys at 12y and all MSM"),
              expression("40% boys at 12y"),
              expression(paste("MSM", phantom(0)<=phantom(0), "26y", sep="")),
              expression(paste("MSM", phantom(0)<=phantom(0), "40y", sep="")),
              expression("all MSM"))
legend("bottomleft", legend=legendexp, lty=1, lwd=2,
 col=cols[c(4,2,9,10,11)], title="Vaccine eligibility:", title.adj=0.1, bty="n")
title("Anal HPV16 prevalence, base-case")
