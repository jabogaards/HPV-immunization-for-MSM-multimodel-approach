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


### PENILE PREDICTIONS
fill.y1_ <- function(age) {
  y.model <- array(NA,dim=c(360,length(time)))
  for(i in 1:20) load(paste("dyn.",models[i],age,".RData",sep=""))
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

### ANAL PREDICTIONS
fill.y_1 <- function(age) {
  y.model <- array(NA,dim=c(360,length(time)))
  for(i in 1:20) load(paste("dyn.",models[i],age,".RData",sep=""))
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

### REDUCTION MATRIX

red.matrix <- array(NA,dim=c(20,18,5,2))

y.model11 <- fill.y1_(age=".b4age27")
red11 <- 1-y.model11[,length(time)]/y.model11[,1]
y.model21 <- fill.y_1(age=".b4age27")
red21 <- 1-y.model21[,length(time)]/y.model21[,1]

y.model12 <- fill.y1_(age=".b4age40")
red12 <- 1-y.model12[,length(time)]/y.model12[,1]
y.model22 <- fill.y_1(age=".b4age40")
red22 <- 1-y.model22[,length(time)]/y.model22[,1]

y.model13 <- fill.y1_(age=".b4age80")
red13 <- 1-y.model13[,length(time)]/y.model13[,1]
y.model23 <- fill.y_1(age=".b4age80")
red23 <- 1-y.model23[,length(time)]/y.model23[,1]

setwd("../boys")

y.model14 <- fill.y1_(age=".age12")
red14 <- 1-y.model14[,length(time)]/y.model14[,1]
y.model24 <- fill.y_1(age=".age12")
red24 <- 1-y.model24[,length(time)]/y.model24[,1]

setwd("../combi")

y.model15 <- fill.y1_(age=".b4age80")
red15 <- 1-y.model15[,length(time)]/y.model15[,1]
y.model25 <- fill.y_1(age=".b4age80")
red25 <- 1-y.model25[,length(time)]/y.model25[,1]

for (model in 1:20) {
  for (listindex in 1:18) {
    red.matrix[model,listindex,1,1] <- red11[(model-1)*18+listindex]
    red.matrix[model,listindex,1,2] <- red21[(model-1)*18+listindex]
    red.matrix[model,listindex,2,1] <- red12[(model-1)*18+listindex]
    red.matrix[model,listindex,2,2] <- red22[(model-1)*18+listindex]
    red.matrix[model,listindex,3,1] <- red13[(model-1)*18+listindex]
    red.matrix[model,listindex,3,2] <- red23[(model-1)*18+listindex]
    red.matrix[model,listindex,4,1] <- red14[(model-1)*18+listindex]
    red.matrix[model,listindex,4,2] <- red24[(model-1)*18+listindex]
    red.matrix[model,listindex,5,1] <- red15[(model-1)*18+listindex]
    red.matrix[model,listindex,5,2] <- red25[(model-1)*18+listindex]
  }
}

##### MAKE FIGURES
windows(width = 12, height = 10)
layout(matrix(1:9, 3, 3, byrow = TRUE), respect = F, heights=c(5,5,1))

table.order <- c(1,2,5,7,6,3,4,8,10,9,14,16,15,11,13,12,17,19,18,20)
modcols <- c(rep(5,2),rep(2,5),rep(8,9),rep(10,3),11)

par(mfg=c(1,1), mar=c(2, 3.5, 3, 0.5))
plot.new(); boxplot(t(red.matrix[table.order,,3,2]), xaxt="n",
 col=cols[modcols]); box()
text(x=c(0.5,5,12,19), y=0.0525, c("SIS","systemic\nimmunity","local\nimmunity", "latency"), cex=1.1, xpd=TRUE)
title("Eligibility: all MSM", cex.main=1.5)

par(mfg=c(1,2), mar=c(2, 3.5, 3, 0.5))
plot.new(); boxplot(t(red.matrix[table.order,,4,2]), xaxt="n",
 col=cols[modcols]); box()
text(x=c(0.5,5,12,19), y=0.445, c("SIS","systemic\nimmunity","local\nimmunity", "latency"), cex=1.1, xpd=TRUE)
title("40% boys at 12y", cex.main=1.5)

par(mfg=c(1,3), mar=c(2, 3.5, 3, 0.5))
plot.new(); boxplot(t(red.matrix[table.order,,5,2]), xaxt="n",
 col=cols[modcols]); box()
text(x=c(0.5,5,12,19), y=0.555, c("SIS","systemic\nimmunity","local\nimmunity", "latency"), cex=1.1, xpd=TRUE)
title("40% boys at 12y + all MSM", cex.main=1.5)

modcols <- c(sort(rep(1:3,2),decreasing=T),sort(rep(4:6,2),decreasing=T),sort(rep(9:11,2)))

par(mfg=c(2,1), mar=c(2, 3.5, 3, 0.5))
plot.new(); boxplot(red.matrix[,,3,2], xaxt="n", col=cols[modcols]); box()

par(mfg=c(2,2), mar=c(2, 3.5, 3, 0.5))
plot.new(); boxplot(red.matrix[,,4,2], xaxt="n", col=cols[modcols]); box()

par(mfg=c(2,3), mar=c(2, 3.5, 3, 0.5))
plot.new(); boxplot(red.matrix[,,5,2], xaxt="n", col=cols[modcols]); box()

par(mfg=c(3,1))
plot.new();
text(x=-0.10, y=5, labels="Assort.", cex=1.1, srt=0, xpd=TRUE)
text(x=c((1+2*(0:8))/18), y=5, labels=c("0","1/3","2/3"), cex=1.1, xpd=TRUE)
text(x=c(6.2/18,12.3/18), y=5, labels="|", cex=1.1, xpd=TRUE)
text(x=-0.10, y=2, labels="Activity", cex=1.1, srt=0, xpd=TRUE)
text(x=c(3/18, 9/18, 15.1/18), y=2, labels=c("80-20%","90-10%","60-30-10%"), cex=1.1, xpd=TRUE)
text(x=c(6.2/18,12.3/18), y=2, labels="|", cex=1.1, xpd=TRUE)

par(mfg=c(3,2))
plot.new();
text(x=-0.10, y=5, labels="Assort.", cex=1.1, srt=0, xpd=TRUE)
text(x=c((1+2*(0:8))/18), y=5, labels=c("0","1/3","2/3"), cex=1.1, xpd=TRUE)
text(x=c(6.2/18,12.3/18), y=5, labels="|", cex=1.1, xpd=TRUE)
text(x=-0.10, y=2, labels="Activity", cex=1.1, srt=0, xpd=TRUE)
text(x=c(3/18, 9/18, 15.1/18), y=2, labels=c("80-20%","90-10%","60-30-10%"), cex=1.1, xpd=TRUE)
text(x=c(6.2/18,12.3/18), y=2, labels="|", cex=1.1, xpd=TRUE)

par(mfg=c(3,3))
plot.new();
text(x=-0.10, y=5, labels="Assort.", cex=1.1, srt=0, xpd=TRUE)
text(x=c((1+2*(0:8))/18), y=5, labels=c("0","1/3","2/3"), cex=1.1, xpd=TRUE)
text(x=c(6.2/18,12.3/18), y=5, labels="|", cex=1.1, xpd=TRUE)
text(x=-0.10, y=2, labels="Activity", cex=1.1, srt=0, xpd=TRUE)
text(x=c(3/18, 9/18, 15.1/18), y=2, labels=c("80-20%","90-10%","60-30-10%"), cex=1.1, xpd=TRUE)
text(x=c(6.2/18,12.3/18), y=2, labels="|", cex=1.1, xpd=TRUE)
