library(urca)
data(npext)
y <- npext$unemploy
y2 <- data.frame(Year=1890:1988,unemp=y)
new <- read.csv('data/unemployment2.csv',header=F)
new <- with(new,data.frame(year=V1,unemp=log(V16)))
y <- c(y,new$unemp[new$year>1988])
y <- ts(y,start=1890)

lrTest <- function(p,edge=TRUE){
    stopifnot(p<20)
    p <- p+1
    big <- ifelse(edge,p+1,21)
    llBig <- logLik(arp[[big]])
    llSmall <- logLik(arp[[p]])
    lrtest <- as.numeric(2*(llBig-llSmall))
    pchisq(lrtest , df = 1, lower.tail = FALSE)
}


arp <- lapply(0:20,function(p) arima(y,order=c(p,0,0)))
pvalsEdge <- vapply(0:19,lrTest,1)
pvalsFull <- vapply(0:19,lrTest,1,edge=FALSE)

aics <- vapply(arp,AIC,1)
bics <- vapply(arp,BIC,1)

np <- length(pvalsEdge)
dhatsEdge <- c(dhatAll(pvalsEdge,backwards=TRUE),AIC=which.min(aics),BIC=which.min(bics))
dhatsTot <- c(dhatAll(pvalsFull,backwards=TRUE),AIC=which.min(aics),BIC=which.min(bics))
