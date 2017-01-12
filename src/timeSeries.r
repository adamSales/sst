require(foreign)
require(vars)
#### do the VAR analysis
#### modeled after Bernanke and Gertler 1995 "Inside the Black Box..."
#### Fig 1 in B&G
#### data downloaded 1/2/17 from FRED

### assemble data
gdp <- read.csv('data/GDP.csv')
gdpdef <- read.csv('data/GDPDEF.csv')
priceIndex <- read.csv('data/PPIACO.csv')
fedfund <- read.csv('data/FEDFUNDS.csv')

tsdat <- merge(gdp,gdpdef)
tsdat <- merge(tsdat,priceIndex)
names(tsdat)[4] <- 'priceIndex'
tsdat <- merge(tsdat,fedfund)
names(tsdat)[5] <- 'fedfund'

tsdat <- as.ts(tsdat)

mod <- VAR(tsdat[,-1],49)

mods <- lapply(1:48,function(p) VAR(tsdat[,-1],p=p))

ports <- list()
for(mm in mods) ports <- c(ports,serial.test(mm,lags.pt=100,type='PT.adjusted')[[2]]$p.value)

ps <- vapply(mods,function(m) m$p,1)


m(1 + jm) + m(m + 1)/2
s <- function(mod){
    sum(vapply(coef(mod),nrow,1))
}
aic <- function(mod){
    -2*logLik(mod)+2*s(mod)
}

lrs <- NULL
for(i in 1:(length(mods)-1)){
    lrs[i] <- (logLik(mods[[i+1]])-logLik(mods[[i]]))
}
lrps <- vapply(lrs,function(lr) pchisq(lr,16,lower.tail=FALSE),1)

totMod <- VAR(tsdat[,-1],lag.max=48,ic='AIC')

mod <- mods[[5]]
sig <- function(mod) crossprod(resid(mod))/(mod$obs)
aaa <- log(det(sig))+2/mod$obs*16*mod$p


lrTest <- function(modBig,modSmall){
    lr <- 2*(logLik(modBig)-logLik(modSmall))
    pchisq(lr,s(modBig)-s(modSmall),lower.tail=FALSE)
}


lrPs <- vapply(48:2,function(i) lrTest(mods[[i]],mods[[i-1]]),1)
lrPs2 <- vapply(47:1,function(i) lrTest(mods[[48]],mods[[i]]),1)
mods <- lapply(1:48,function(p) VAR(tsdat[,-c(1,2)],p=p))


modsCan <- lapply(1:20,function(p) VAR(Canada,p=p,type='both'))

lrPsCan <- vapply(15:2,function(i) lrTest(modsCan[[i]],modsCan[[i-1]]),1)

stDat <- read.dta('data/lutkepohl2.dta')
dddd <- stDat[-1,c('dln_inv', 'dln_inc', 'dln_consump')]
modsSt <- lapply(1:40,function(p) VAR(dddd,p=p))


makePs <- function(mods,maxp){
    res <- NULL
    for(p in maxp:2){
        res <- try(lrTest(mods[[p]],mods[[p-1]]))
    }
    res
}
