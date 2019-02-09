library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### SITE-SPECIFIC TRANSMISSIBILITES: UNPACK beta.estimates

load("beta1.matrix.RData")
load("beta2.matrix.RData")

names.arg = c(
  "SIS","SISPS","SISPminRS","SISPmaxRS",
  "SIRS","SIS10RS","SIS33RS",
  "SIRlocalS","SIS10RlocalS","SIS33RlocalS",
  "SIRanalS","SIS10RanalS","SIS33RanalS",
  "SIRpenileS","SIS10RpenileS","SIS33RpenileS",
  "SIL","SIS10L","SIS33L","SIR33L"
)

table.order <- c(1,2,5,7,6,3,4,8,10,9,14,16,15,11,13,12,17,19,18,20)
modcols <- c(rep(5,2),rep(2,5),rep(8,9),rep(10,3),rep(11,1))

### PLOT ESTIMATES
x11()
boxplot(t(beta1.matrix[table.order,]), xaxt="n", ylab="", col=cols[modcols])
text(x=c(1,5,12,18), y=-0.13, c("SIS","systemic\nimmunity","local\nimmunity", "latency"), cex=0.9, xpd=TRUE)
title("Penile-to-anal transmissibility")

x11()
boxplot(t(beta2.matrix[table.order,]), xaxt="n",ylab="", col=cols[modcols])
text(x=c(1,5,12,18), y=-0.033, c("SIS","systemic\nimmunity","local\nimmunity", "latency"), cex=0.9, xpd=TRUE)
title("Anal-to-penile transmissibility")

modcols <- c(sort(rep(1:3,2),decreasing=T),sort(rep(4:6,2),decreasing=T),sort(rep(9:11,2)))

x11()
boxplot(beta1.matrix[1:20,], xaxt="n", ylab="", col=cols[modcols])
text(x=-0.5, y=-0.1, labels="Assortativity", cex=0.8, srt=0, xpd=TRUE)
text(x=c(1.5+2*(0:8)), y=-0.1, labels=c("0","1/3","2/3"), cex=0.8, xpd=TRUE)
text(x=c(6.5,12.5), y=-0.1, labels="|", xpd=TRUE)
text(x=0, y=-0.15, labels="Activity", cex=0.8, srt=0, xpd=TRUE)
text(x=c(3.5, 9.5, 15.5), y=-0.15, labels=c("80-20%","90-10%","60-30-10%"), cex=0.8, xpd=TRUE)
text(x=c(6.5,12.5), y=-0.15, labels="|", xpd=TRUE)
title("Penile-to-anal transmissibility")

x11()
boxplot(beta2.matrix[1:20,], xaxt="n", ylab="", col=cols[modcols])
text(x=-0.5, y=-0.026, labels="Assortativity", cex=0.8, srt=0, xpd=TRUE)
text(x=c(1.5+2*(0:8)), y=-0.026, labels=c("0","1/3","2/3"), cex=0.8, xpd=TRUE)
text(x=c(6.5,12.5), y=-0.026, labels="|", xpd=TRUE)
text(x=0, y=-0.0375, labels="Activity", cex=0.8, srt=0, xpd=TRUE)
text(x=c(3.5, 9.5, 15.5), y=-0.0375, labels=c("80-20%","90-10%","60-30-10%"), cex=0.8, xpd=TRUE)
text(x=c(6.5,12.5), y=-0.0375, labels="|", xpd=TRUE)
title("Anal-to-penile transmissibility")
