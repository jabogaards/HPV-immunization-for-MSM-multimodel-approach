library("foreign")

sexpartners <- read.csv("sexpartners.csv",header=T)

table(sexpartners$analsx6m)
table(sexpartners$nranpa6m)

summary(sexpartners$nranpa6m[sexpartners$analsx6m==1])

catanpa6m <- ifelse(sexpartners$nranpa6m == 0, 0, NA)
catanpa6m[sexpartners$nranpa6m == 1] <- 1
catanpa6m[sexpartners$nranpa6m == 2] <- 2
catanpa6m[sexpartners$nranpa6m == 3] <- 3
catanpa6m[sexpartners$nranpa6m == 4] <- 4
catanpa6m[sexpartners$nranpa6m == 5] <- 5
catanpa6m[sexpartners$nranpa6m >=  6 & sexpartners$nranpa6m < 10] <- 6
catanpa6m[sexpartners$nranpa6m >= 10 & sexpartners$nranpa6m < 15] <- 10
catanpa6m[sexpartners$nranpa6m >= 15 & sexpartners$nranpa6m < 25] <- 15
catanpa6m[sexpartners$nranpa6m >= 25] <- 25
table(catanpa6m)

table(catanpa6m, sexpartners$analsxpos)
chisq.test(table(catanpa6m[catanpa6m>=1], sexpartners$analsxpos[catanpa6m>=1]))
chisq.test(table(catanpa6m[catanpa6m>=2], sexpartners$analsxpos[catanpa6m>=2]))

analsxpos <- array(NA, dim=c(9,3))
index <- 0
for (i in 1:25) {
  temp <- table(sexpartners$analsxpos[catanpa6m==i])
  if (sum(temp)) {
    index <- index + 1
    for (j in 1:3) analsxpos[index, j] <- temp[j]/sum(temp)
  }
}

# likelihood for saturated model
x <- table(sexpartners$nranpa6m[catanpa6m>=1],sexpartners$analsxpos[catanpa6m>=1])
analsxp <- array(NA, dim=dim(x))
for (i in 1:dim(x)[1]) analsxp[i, ] <- x[i, ]/apply(x, 1, sum)[i]
ll_satmod = x * log(analsxp)
sum(ll_satmod, na.rm=T)

# likelihood for homogeneous model
x <- table(sexpartners$nranpa6m[catanpa6m>=1],sexpartners$analsxpos[catanpa6m>=1])
n <- as.numeric(rownames(x))
p = 0.49
log_likelihood = array(NA, dim=dim(x))
log_likelihood[,1] = x[,1] * log(p^n)
log_likelihood[,2] = x[,2] * log(p^n)
log_likelihood[,3] = x[,3] * log(1 - 2 * p^n)
sum(log_likelihood)

# write a function and optimize() to find optimum
grrr <- function(p) {
  log_likelihood[,1] = x[,1] * log(p^n)
  log_likelihood[,2] = x[,2] * log(p^n)
  log_likelihood[,3] = x[,3] * log(1 - 2 * p^n)
  sum(-log_likelihood)
}
optimize(grrr, interval=c(0,0.5))
pi <- optimize(grrr, interval=c(0,0.5))$minimum
ll <- -optimize(grrr, interval=c(0,0.5))$objective

# simulated data
xsim <- x
xsim[,1] <- apply(x, 1, sum) * pi^n
xsim[,2] <- apply(x, 1, sum) * pi^n
xsim[,3] <- apply(x, 1, sum) * (1 - 2 * pi^n)

analsxsim <- array(NA, dim=c(9,3))
for (i in 1:5) analsxsim[i, ] <- xsim[i, ]
analsxsim[6,] <- apply(xsim[n>=6  & n<10, ], 2, sum)
analsxsim[7,] <- apply(xsim[n>=10 & n<15, ], 2, sum)
analsxsim[8,] <- apply(xsim[n>=15 & n<25, ], 2, sum)
analsxsim[9,] <- apply(xsim[n>=25, ], 2, sum)

analsxdat <- array(NA, dim=c(9,3))
index <- 0
for (i in 1:25) {
  temp <- table(sexpartners$analsxpos[catanpa6m==i])
  if (sum(temp)) {
    index <- index + 1
    analsxdat[index, ] <- temp
  }
}

x11()
barplot(t(analsxdat), beside=TRUE,
 names.arg=c("1","2","3","4","5","6-10","10-15","15-25","25+"),
 legend.text=c("only active",
               "only passive",
               "both active and passive"))
title("Self-reported behaviour")
title(xlab="No. anal sex partners/6m", ylab="N (respondents)")

# likelihood for heterogeneous model
x <- table(sexpartners$nranpa6m[catanpa6m>=1],sexpartners$analsxpos[catanpa6m>=1])
n <- as.numeric(rownames(x))
p = 0.2
q = 0.2
log_likelihood = array(NA, dim=dim(x))
log_likelihood[,1] = x[,1] * log(q + (1 - 2*q) * p^n)
log_likelihood[,2] = x[,2] * log(q + (1 - 2*q) * p^n)
log_likelihood[,3] = x[,3] * log(1 - 2*q - 2 * (1 - 2*q) * p^n)
sum(log_likelihood)

# write a function and optim() to find minimum
brrr <- function(pars) {
  p = pars[1]
  q = pars[2]
  log_likelihood[,1] = x[,1] * log(q + (1 - 2*q) * p^n)
  log_likelihood[,2] = x[,2] * log(q + (1 - 2*q) * p^n)
  log_likelihood[,3] = x[,3] * log(1 - 2*q - 2 * (1 - 2*q) * p^n)
  sum(-log_likelihood)
}
optim(c(0.2,0.2), brrr)

qu  <- optim(c(0.2,0.2), brrr)$par[2]
pi2 <- optim(c(0.2,0.2), brrr)$par[1]
ll2 <- -optim(c(0.2,0.2), brrr)$value

# simulated data
xsim <- x
xsim[,1] <- apply(x, 1, sum) * (qu + (1 - 2*qu) * pi2^n)
xsim[,2] <- apply(x, 1, sum) * (qu + (1 - 2*qu) * pi2^n)
xsim[,3] <- apply(x, 1, sum) * (1 - 2*qu - 2 * (1 - 2*qu) * pi2^n)

analsxsim2 <- array(NA, dim=c(9,3))
for (i in 1:5) analsxsim2[i, ] <- xsim[i, ]
analsxsim2[6,] <- apply(xsim[n>=6  & n<10, ], 2, sum)
analsxsim2[7,] <- apply(xsim[n>=10 & n<15, ], 2, sum)
analsxsim2[8,] <- apply(xsim[n>=15 & n<25, ], 2, sum)
analsxsim2[9,] <- apply(xsim[n>=25, ], 2, sum)

x11()
barplot(t(analsxsim2), beside=TRUE,
 names.arg=c("1","2","3","4","5","6-10","10-15","15-25","25+"),
 legend.text=c("only active",
               "only passive",
               "both active and passive"))
title("Mixture model predictions")
title(xlab="No. anal sex partners/6m", ylab="N (predicted respondents)")

# think of a goodness-of-fit test...
deviance_ll <- -2*ll - -2*sum(ll_satmod, na.rm=T)
df_ll  <- dim(ll_satmod)[1] * 2 - 1
1 - pchisq(q=deviance_ll, df=df_ll)

deviance_ll2 <- -2*ll2 - -2*sum(ll_satmod, na.rm=T)
df_ll2  <- dim(ll_satmod)[1] * 2 - 2
1 - pchisq(q=deviance_ll2, df=df_ll2)

# likelihood for restricted model
x <- table(sexpartners$nranpa6m[catanpa6m>=1],sexpartners$analsxpos[catanpa6m>=1])
analsxp2 <- array(NA, dim=dim(x))
for (i in 1:dim(x)[1]) {
  analsxp2[i,1] <- 0.5*(x[i,1] + x[i,2])/apply(x, 1, sum)[i]
  analsxp2[i,2] <- analsxp2[i,1]
  analsxp2[i,3] <- x[i,3]/apply(x, 1, sum)[i]
}
ll_resmod = x * log(analsxp2)
sum(ll_resmod, na.rm=T)

deviance_resmod <- -2*sum(ll_resmod, na.rm=T) - -2*sum(ll_satmod, na.rm=T)
df_resmod  <- dim(ll_resmod)[1]
1 - pchisq(q=deviance_resmod, df=df_resmod)

deviance_dif <- deviance_ll2 - deviance_resmod
df_dif <- df_ll2 - df_resmod
1 - pchisq(q=deviance_dif, df=df_dif)

# alternative: goodness-of-fit on frequency table

analsxdat # observed freq table
analsxsim2 # expected freq table
X2 <- sum((abs(analsxdat-analsxsim2)-0.5)^2/analsxsim2)
1 - pchisq(q=X2, df=dim(analsxdat)[1]*2)
