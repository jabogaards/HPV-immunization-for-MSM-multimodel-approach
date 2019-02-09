library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

##### PROJECTED REDUCTIONS: see Table 3 in main text

EFFECTIVENESS <- array(NA,dim=c(13,4))
colnames(EFFECTIVENESS) <- c("target","mode","coverage","reduction")

EFFECTIVENESS <- as.data.frame(EFFECTIVENESS)
EFFECTIVENESS[,1] <-
 c(rep("MSM26",3),rep("MSM40",3),rep("MSMall",3),rep("boys12",3),"combi")
EFFECTIVENESS[,2] <-
 c("base","improved","leaky","base","improved","leaky","base","improved","leaky","base","improved","leaky","base")
EFFECTIVENESS[,3] <-
 as.numeric(c(9.4,17.6,9.4,19.2,33.8,19.2,21.2,36.5,21.2,39.7,79.3,39.7,52.4))
EFFECTIVENESS[,4] <-
 as.numeric(c(13.4,24.6,6.1,25.5,43.6,11.0,27.3,46.1,11.7,61.6,97.4,57.9,74.8))

EFFECTIVENESS$symbol <-
 ifelse(EFFECTIVENESS$mode=="leaky",1,
  ifelse(EFFECTIVENESS$mode=="base",19,15))

EFFECTIVENESS$col <-
 ifelse(EFFECTIVENESS$target=="combi",4,
  ifelse(EFFECTIVENESS$target=="boys12",2,
   ifelse(EFFECTIVENESS$target=="MSM26",9,
    ifelse(EFFECTIVENESS$target=="MSM40",10,11))))

plot(reduction ~ coverage, pch=symbol, col=cols[col], cex=1.33,
 xlim=c(0,100), xlab="Vaccine coverage among MSM (%)",
  ylim=c(0,100), ylab="Anal HPV16 prevalence reduction (%)",
   data=EFFECTIVENESS)
abline(a=0,b=1,lty=3)

legendexp = c(expression("    boys at 12y and all MSM"),
              expression("    boys at 12y"),
              expression(paste("    MSM", phantom(0)<=phantom(0), "26y", sep="")),
              expression(paste("    MSM", phantom(0)<=phantom(0), "40y", sep="")),
              expression("    all MSM"))
legend("topleft", legend=legendexp, pch=19,
 col=cols[c(4,2,9,10,11)], title="Vaccine eligibility:", title.adj=0.1, bty="n")
title("Projected reductions versus HPV vaccine coverage")

legend(x=-2, y=104, legend=rep(" ",5), pch=1,
 col=c("white",cols[c(2,9,10,11)]), title="", bty="n")

legend(x=0, y=103.9, legend=rep(" ",5), pch=15,
 col=c("white",cols[c(2,9,10,11)]), title="", bty="n")
