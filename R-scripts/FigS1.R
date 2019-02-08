library(MASS)
library(sn)

sexage <- read.csv("sexage.csv",header=T)
y <- sexage$agefasm[!is.na(sexage$agefasm)]
summary(y)

hist(y, freq=TRUE, breaks=c(0:30)*2)
# distribution skewed to the right, try log transformation
hist(log(y), freq=TRUE, breaks=30)
# consider sexual debut at 4 years to be an outlier
hist(log(y[y>=10]), freq=TRUE, breaks=20)

# binning data for goodness-of-fit test
y.cut <- cut(y, breaks=c(0,16,20,24,30,35,40,60))
table(y.cut)

##### first fit lognormal distribution to the data
lnorm.fit <- fitdistr(log(y[y>=10]), densfun="normal")
lnorm.logy <- rnorm(10^6, mean=lnorm.fit$estimate[1], sd=lnorm.fit$estimate[2])
lnorm.y <- exp(lnorm.logy)
hist(y, breaks=c(0:30)*2, freq=FALSE, col="pink",
 xlab="Age first anal sex with male partner", main="MSM population entrance")
hist(lnorm.y, breaks=c(0:50)*2, freq=FALSE, add=TRUE, dens=10)
# goodness-of-fit?
lnorm.cut <- cut(lnorm.y, breaks=c(0,16,20,24,30,35,40,max(lnorm.y)))
lnorm.exp <- table(lnorm.cut)*10^-6*length(y)
lnorm.X2 <- sum( (table(y.cut) - lnorm.exp)^2 / lnorm.exp )
lnorm.df <- length(lnorm.exp)-1-length(lnorm.fit$estimate)
1 - pchisq(lnorm.X2, lnorm.df)

##### consider lognormal distribution on a log scale
llnorm.fit <- fitdistr(log(y[y>=10]), densfun="lognormal")
llnorm.logy <- rlnorm(10^6, meanlog=llnorm.fit$estimate[1], sdlog=llnorm.fit$estimate[2])
llnorm.y <- exp(llnorm.logy)
hist(y, breaks=c(0:30)*2, freq=FALSE, col="pink",
 xlab="Age first anal sex with male partner", main="MSM population entrance")
hist(llnorm.y, breaks=c(0:50)*2, freq=FALSE, add=TRUE, dens=10)
# goodness-of-fit?
llnorm.cut <- cut(llnorm.y, breaks=c(0,16,20,24,30,35,40,max(llnorm.y)))
llnorm.exp <- table(llnorm.cut)*10^-6*length(y)
llnorm.X2 <- sum( (table(y.cut) - llnorm.exp)^2 / llnorm.exp )
llnorm.df <- length(llnorm.exp)-1-length(llnorm.fit$estimate)
1 - pchisq(llnorm.X2, llnorm.df)

##### try skew-normal distribution on a log scale
x <- rep(1, length(y[y>=10]))
slnorm.fit <- selm.fit(as.matrix(x), log(y[y>=10]))$param
slnorm.logy <- rsn(10^6, xi=slnorm.fit$dp[1], omega=slnorm.fit$dp[2], alpha=slnorm.fit$dp[3])
slnorm.y <- exp(slnorm.logy)
hist(y, breaks=c(0:30)*2, freq=FALSE, col="pink",
 xlab="Age first anal sex with male partner", main="MSM population entrance")
hist(slnorm.y, breaks=c(0:50)*2, freq=FALSE, add=TRUE, dens=10)
# goodness-of-fit?
slnorm.cut <- cut(slnorm.y, breaks=c(0,16,20,24,30,35,40,max(slnorm.y)))
slnorm.exp <- table(slnorm.cut)*10^-6*length(y)
slnorm.X2 <- sum( (table(y.cut) - slnorm.exp)^2 / slnorm.exp )
slnorm.df <- length(slnorm.exp)-1-length(slnorm.fit$estimate)
1 - pchisq(slnorm.X2, slnorm.df)
