library("foreign")

##### H2M DATA
baseline <- read.csv("baseline.csv",header=T)
select <- !is.na(baseline$hpv16_anal) &  !is.na(baseline$hpv16_penile)

success <- sum(baseline$hpv16_penile[select] & !baseline$hpv16_anal[select])
prev16_penile <- success / sum(select)
prev16_penile_lcl <- binom.test(success, sum(select), conf.level=0.95)$conf.int[1]
prev16_penile_ucl <- binom.test(success, sum(select), conf.level=0.95)$conf.int[2]

success <-  sum(baseline$hpv16_anal[select] & !baseline$hpv16_penile[select])
prev16_anal <- success / sum(select)
prev16_anal_lcl <- binom.test(success, sum(select), conf.level=0.95)$conf.int[1]
prev16_anal_ucl <- binom.test(success, sum(select), conf.level=0.95)$conf.int[2]

success <-  sum(baseline$hpv16_anal[select] & baseline$hpv16_penile[select])
prev16_both <- success / sum(select)
prev16_both_lcl <- binom.test(success, sum(select), conf.level=0.95)$conf.int[1]
prev16_both_ucl <- binom.test(success, sum(select), conf.level=0.95)$conf.int[2]

##### MODELS: UNPACK base.models.zip

models <- c("SIS","SISPS","SISPminRS","SISPmaxRS","SIRS","SIS10RS","SIS33RS","SIRlocalS","SIS10RlocalS","SIS33RlocalS",
		"SIRanalS","SIS10RanalS","SIS33RanalS","SIRpenileS","SIS10RpenileS","SIS33RpenileS","SIL","SIS10L","SIS33L","SIR33L")

for(i in 1:length(models)) load(paste("base.",models[i],".RData",sep=""))

load("delta.AIC.RData")
aic.weights <- exp(-0.5*delta.aic)/sum(exp(-0.5*delta.aic))

minAge <- 9
maxAge <- 80
gridAge <- 1/4
axisAge <- (minAge/gridAge):(maxAge/gridAge)*gridAge
x_a <- round(baseline$age)[select]

##### PENILE-ONLY FIT
y10.model <- array(NA,dim=c(360,sum(select)))
for (listindex in 1:length(base.SIS)) {
 for (pid in 1:sum(select)) {
  y10.model[0*length(base.SIS)+listindex,pid] <- base.SIS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[1*length(base.SIS)+listindex,pid] <- base.SISPS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[2*length(base.SIS)+listindex,pid] <- base.SISPminRS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[3*length(base.SIS)+listindex,pid] <- base.SISPmaxRS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[4*length(base.SIS)+listindex,pid] <- base.SIRS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[5*length(base.SIS)+listindex,pid] <- base.SIS10RS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[6*length(base.SIS)+listindex,pid] <- base.SIS33RS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[7*length(base.SIS)+listindex,pid] <- base.SIRlocalS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[8*length(base.SIS)+listindex,pid] <- base.SIS10RlocalS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[9*length(base.SIS)+listindex,pid] <- base.SIS33RlocalS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[10*length(base.SIS)+listindex,pid] <- base.SIRanalS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[11*length(base.SIS)+listindex,pid] <- base.SIS10RanalS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[12*length(base.SIS)+listindex,pid] <- base.SIS33RanalS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[13*length(base.SIS)+listindex,pid] <- base.SIRpenileS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[14*length(base.SIS)+listindex,pid] <- base.SIS10RpenileS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[15*length(base.SIS)+listindex,pid] <- base.SIS33RpenileS[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[16*length(base.SIS)+listindex,pid] <- base.SIL[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[17*length(base.SIS)+listindex,pid] <- base.SIS10L[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[18*length(base.SIS)+listindex,pid] <- base.SIS33L[[listindex]]$mod16_x10[axisAge==x_a[pid]]
  y10.model[19*length(base.SIS)+listindex,pid] <- base.SIR33L[[listindex]]$mod16_x10[axisAge==x_a[pid]]
 }
}
y10_mean <- sum(aic.weights*apply(y10.model,1,sum)/sum(select))
y10_min <- min(apply(y10.model,1,sum)/sum(select))
y10_max <- max(apply(y10.model,1,sum)/sum(select))

##### ANAL-ONLY FIT
y01.model <- array(NA,dim=c(360,sum(select)))
for (listindex in 1:length(base.SIS)) {
 for (pid in 1:sum(select)) {
  y01.model[0*length(base.SIS)+listindex,pid] <- base.SIS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[1*length(base.SIS)+listindex,pid] <- base.SISPS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[2*length(base.SIS)+listindex,pid] <- base.SISPminRS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[3*length(base.SIS)+listindex,pid] <- base.SISPmaxRS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[4*length(base.SIS)+listindex,pid] <- base.SIRS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[5*length(base.SIS)+listindex,pid] <- base.SIS10RS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[6*length(base.SIS)+listindex,pid] <- base.SIS33RS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[7*length(base.SIS)+listindex,pid] <- base.SIRlocalS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[8*length(base.SIS)+listindex,pid] <- base.SIS10RlocalS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[9*length(base.SIS)+listindex,pid] <- base.SIS33RlocalS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[10*length(base.SIS)+listindex,pid] <- base.SIRanalS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[11*length(base.SIS)+listindex,pid] <- base.SIS10RanalS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[12*length(base.SIS)+listindex,pid] <- base.SIS33RanalS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[13*length(base.SIS)+listindex,pid] <- base.SIRpenileS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[14*length(base.SIS)+listindex,pid] <- base.SIS10RpenileS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[15*length(base.SIS)+listindex,pid] <- base.SIS33RpenileS[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[16*length(base.SIS)+listindex,pid] <- base.SIL[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[17*length(base.SIS)+listindex,pid] <- base.SIS10L[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[18*length(base.SIS)+listindex,pid] <- base.SIS33L[[listindex]]$mod16_x01[axisAge==x_a[pid]]
  y01.model[19*length(base.SIS)+listindex,pid] <- base.SIR33L[[listindex]]$mod16_x01[axisAge==x_a[pid]]
 }
}
y01_mean <- sum(aic.weights*apply(y01.model,1,sum)/sum(select))
y01_min <- min(apply(y01.model,1,sum)/sum(select))
y01_max <- max(apply(y01.model,1,sum)/sum(select))

##### PENILE+ANAL FIT
y11.model <- array(NA,dim=c(360,sum(select)))
for (listindex in 1:length(base.SIS)) {
 for (pid in 1:sum(select)) {
  y11.model[0*length(base.SIS)+listindex,pid] <- base.SIS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[1*length(base.SIS)+listindex,pid] <- base.SISPS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[2*length(base.SIS)+listindex,pid] <- base.SISPminRS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[3*length(base.SIS)+listindex,pid] <- base.SISPmaxRS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[4*length(base.SIS)+listindex,pid] <- base.SIRS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[5*length(base.SIS)+listindex,pid] <- base.SIS10RS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[6*length(base.SIS)+listindex,pid] <- base.SIS33RS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[7*length(base.SIS)+listindex,pid] <- base.SIRlocalS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[8*length(base.SIS)+listindex,pid] <- base.SIS10RlocalS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[9*length(base.SIS)+listindex,pid] <- base.SIS33RlocalS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[10*length(base.SIS)+listindex,pid] <- base.SIRanalS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[11*length(base.SIS)+listindex,pid] <- base.SIS10RanalS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[12*length(base.SIS)+listindex,pid] <- base.SIS33RanalS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[13*length(base.SIS)+listindex,pid] <- base.SIRpenileS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[14*length(base.SIS)+listindex,pid] <- base.SIS10RpenileS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[15*length(base.SIS)+listindex,pid] <- base.SIS33RpenileS[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[16*length(base.SIS)+listindex,pid] <- base.SIL[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[17*length(base.SIS)+listindex,pid] <- base.SIS10L[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[18*length(base.SIS)+listindex,pid] <- base.SIS33L[[listindex]]$mod16_x11[axisAge==x_a[pid]]
  y11.model[19*length(base.SIS)+listindex,pid] <- base.SIR33L[[listindex]]$mod16_x11[axisAge==x_a[pid]]
 }
}
y11_mean <- sum(aic.weights*apply(y11.model,1,sum)/sum(select))
y11_min <- min(apply(y11.model,1,sum)/sum(select))
y11_max <- max(apply(y11.model,1,sum)/sum(select))

##### PLOT MEANS
plot(c(prev16_penile, y10_mean, prev16_anal, y01_mean, prev16_both, y11_mean) ~ c(1,2,4,5,7,8),
  pch=19, col=c(1,2), xaxt="n", xlab="", xlim=c(0,9), ylim=c(0,0.15), ylab="")
text(x=c(1.5,4.5,7.5), y=-0.025, labels=c("penile-only","anal-only","penile+anal"), cex=1.5, xpd=TRUE)

segments(x0=1, x1=1, y0=prev16_penile_lcl, y1=prev16_penile_ucl)
segments(x0=4, x1=4, y0=prev16_anal_lcl,   y1=prev16_anal_ucl)
segments(x0=7, x1=7, y0=prev16_both_lcl,   y1=prev16_both_ucl)

segments(x0=2, x1=2, y0=y10_min, y1=y10_max, col="red")
segments(x0=5, x1=5, y0=y01_min, y1=y01_max, col="red")
segments(x0=8, x1=8, y0=y11_min, y1=y11_max, col="red")

legend("topright", legend=c("H2M baseline data","Model predictions"),
 lty=1, pch=19, cex=1.25, col=c(1,2), bty="n")
title("Site-specific HPV16 prevalence", cex.main=1.5)
