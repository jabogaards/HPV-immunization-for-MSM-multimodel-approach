library("foreign")

sexage <- read.csv("sexage.csv",header=T)
sexage$sexualage <- sexage$age - sexage$agefasm

scatterhist = function(x, y, xmin, xmax){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, breaks=20, plot=FALSE)
  yhist = hist(y, breaks=20, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3.3,3.3,0,0))
  plot(x,y, xlim=c(xmin, xmax), ylab="", xlab="")
  par(mar=c(1,3.3,1,1))
  barplot(xhist$counts, axes=T, ylim=c(0, top), space=0)
  par(mar=c(3.3,1,1,1))
  barplot(yhist$counts, axes=T, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3.3,3.3,1,1))
}

with(sexage, scatterhist(agefasm, sexualage, xmin=9, xmax=60))
mtext("Sexual age at baseline", side=2, line=1, outer=TRUE, adj=0, 
    at=(.8 * 1/3))
mtext("Age first anal sex with male partner", side=1, line=1, outer=TRUE,
 adj=0, at=(.8 * 1/3))
