library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### POPULATION AGE DISTRIBUTION: USE sisprs.output (in sisprs.zip)
load("sisprs.output.RData")
N_age <- rowSums(output$N[dim(output$N)[1],,])

minAge <- 9
maxAge <- 80
gridAge <- 1/4
axisAge <- (minAge/gridAge):(maxAge/gridAge)*gridAge

##### MODELS: UNPACK base.models.zip

models <- c("SIS","SISPS","SISPminRS","SISPmaxRS","SIRS","SIS10RS","SIS33RS","SIRlocalS","SIS10RlocalS","SIS33RlocalS",
		"SIRanalS","SIS10RanalS","SIS33RanalS","SIRpenileS","SIS10RpenileS","SIS33RpenileS","SIL","SIS10L","SIS33L","SIR33L")

for(i in 1:length(models)) load(paste("base.",models[i],".RData",sep=""))

### PENILE occurrence
plot.new(); plot.window(xlim = c(15, 75), ylim = c(0, 0.00035)); axis(1); box()
title("Penile HPV16 prevalence (weighted)")
mtext("Age", side=1, line=3)
mtext("Density", side=2, line=2)

for (listindex in 1:18) lines(c(N_age*(base.SIS[[listindex]]$mod16_x10+base.SIS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SISPS[[listindex]]$mod16_x10+base.SISPS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SISPminRS[[listindex]]$mod16_x10+base.SISPminRS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SISPmaxRS[[listindex]]$mod16_x10+base.SISPmaxRS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRS[[listindex]]$mod16_x10+base.SIRS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RS[[listindex]]$mod16_x10+base.SIS10RS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RS[[listindex]]$mod16_x10+base.SIS33RS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRlocalS[[listindex]]$mod16_x10+base.SIRlocalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RlocalS[[listindex]]$mod16_x10+base.SIS10RlocalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RlocalS[[listindex]]$mod16_x10+base.SIS33RlocalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRanalS[[listindex]]$mod16_x10+base.SIRanalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RanalS[[listindex]]$mod16_x10+base.SIS10RanalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RanalS[[listindex]]$mod16_x10+base.SIS33RanalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRpenileS[[listindex]]$mod16_x10+base.SIRpenileS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RpenileS[[listindex]]$mod16_x10+base.SIS10RpenileS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RpenileS[[listindex]]$mod16_x10+base.SIS33RpenileS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIL[[listindex]]$mod16_x10+base.SIL[[listindex]]$mod16_x11))~axisAge, col=cols[9])
for (listindex in 1:18) lines(c(N_age*(base.SIS10L[[listindex]]$mod16_x10+base.SIS10L[[listindex]]$mod16_x11))~axisAge, col=cols[9])
for (listindex in 1:18) lines(c(N_age*(base.SIS33L[[listindex]]$mod16_x10+base.SIS33L[[listindex]]$mod16_x11))~axisAge, col=cols[9])
for (listindex in 1:18) lines(c(N_age*(base.SIR33L[[listindex]]$mod16_x10+base.SIR33L[[listindex]]$mod16_x11))~axisAge, col=cols[9])
box()

### ANAL occurrence
plot.new(); plot.window(xlim = c(15, 75), ylim = c(0, 0.0011)); axis(1); box()
title("Anal HPV16 prevalence (weighted)")
mtext("Age", side=1, line=3)
mtext("Density", side=2, line=2)

for (listindex in 1:18) lines(c(N_age*(base.SIS[[listindex]]$mod16_x01+base.SIS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SISPS[[listindex]]$mod16_x01+base.SISPS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SISPminRS[[listindex]]$mod16_x01+base.SISPminRS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SISPmaxRS[[listindex]]$mod16_x01+base.SISPmaxRS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRS[[listindex]]$mod16_x01+base.SIRS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RS[[listindex]]$mod16_x01+base.SIS10RS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RS[[listindex]]$mod16_x01+base.SIS33RS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRlocalS[[listindex]]$mod16_x01+base.SIRlocalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RlocalS[[listindex]]$mod16_x01+base.SIS10RlocalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RlocalS[[listindex]]$mod16_x01+base.SIS33RlocalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRanalS[[listindex]]$mod16_x01+base.SIRanalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RanalS[[listindex]]$mod16_x01+base.SIS10RanalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RanalS[[listindex]]$mod16_x01+base.SIS33RanalS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIRpenileS[[listindex]]$mod16_x01+base.SIRpenileS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS10RpenileS[[listindex]]$mod16_x01+base.SIS10RpenileS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIS33RpenileS[[listindex]]$mod16_x01+base.SIS33RpenileS[[listindex]]$mod16_x11))~axisAge, col=cols[8])
for (listindex in 1:18) lines(c(N_age*(base.SIL[[listindex]]$mod16_x01+base.SIL[[listindex]]$mod16_x11))~axisAge, col=cols[9])
for (listindex in 1:18) lines(c(N_age*(base.SIS10L[[listindex]]$mod16_x01+base.SIS10L[[listindex]]$mod16_x11))~axisAge, col=cols[9])
for (listindex in 1:18) lines(c(N_age*(base.SIS33L[[listindex]]$mod16_x01+base.SIS33L[[listindex]]$mod16_x11))~axisAge, col=cols[9])
for (listindex in 1:18) lines(c(N_age*(base.SIR33L[[listindex]]$mod16_x01+base.SIR33L[[listindex]]$mod16_x11))~axisAge, col=cols[9])
box()
