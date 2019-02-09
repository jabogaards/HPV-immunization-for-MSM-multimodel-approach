library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### MODEL-SPECIFIC AKAIKE WEIGHTS: EXTRACT delta.AIC (e.g. base.models.zip)

load("delta.AIC.RData")
aic.weights <- exp(-0.5*delta.aic)/sum(exp(-0.5*delta.aic))

midpts = c(9+(0:19)*19)
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

x11()
boxplot(t(delta.aic[table.order,]), xaxt="n", ylab="delta AIC", col=cols[modcols])
text(x=1:20, y=-2, names.arg[table.order], cex=0.8, srt=90, xpd=TRUE)

modcols <- c(sort(rep(1:3,2),decreasing=T),sort(rep(4:6,2),decreasing=T),sort(rep(9:11,2)))

x11()
boxplot(delta.aic[table.order,], xaxt="n", ylab="delta AIC", col=cols[modcols])
text(x=-0.5, y=-1.6, labels="Assortativity", cex=0.8, srt=0, xpd=TRUE)
text(x=c(1.5+2*(0:8)), y=-1.6, labels=c("0","1/3","2/3"), cex=0.8, xpd=TRUE)
text(x=c(6.5,12.5), y=-1.6, labels="|", xpd=TRUE)
text(x=0, y=-2.25, labels="Activity", cex=0.8, srt=0, xpd=TRUE)
text(x=c(3.5, 9.5, 15.5), y=-2.25, labels=c("80-20%","90-10%","60-30-10%"), cex=0.8, xpd=TRUE)
text(x=c(6.5,12.5), y=-2.25, labels="|", xpd=TRUE)
