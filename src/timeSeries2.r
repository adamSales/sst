
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
pvalsTS <- vapply(0:19,lrTest,1)
#pvalsFull <- vapply(0:19,lrTest,1,edge=FALSE)

aics <- vapply(arp,AIC,1)
bics <- vapply(arp,BIC,1)

npTS <- length(pvalsTS)
dhatsTS <- c(dhatAll(pvalsTS,backwards=TRUE,alphas=c(0.05,0.15)),AIC=which.min(aics),BIC=which.min(bics))
#dhatsTot <- c(dhatAll(pvalsFull,backwards=TRUE),AIC=which.min(aics),BIC=which.min(bics))


## unemp <- read.csv('data/unemployment.csv')
## unempVec <- NULL
## for(i in 1:nrow(unemp)) unempVec <- c(unempVec,unemp[i,-1])
## unempVec <- unlist(unempVec)



## unemp <- ts(log(unempVec))
## arp <- lapply(0:50,function(p) arima(unemp,order=c(p,0,0)))
## pvalsEdge <- vapply(0:19,lrTest,1)
## pvalsFull <- vapply(0:19,lrTest,1,edge=FALSE)

## aics <- vapply(arp,AIC,1)
## bics <- vapply(arp,BIC,1)

## np <- length(pvalsEdge)
## dhatsEdge <- c(dhatAll(pvalsEdge,backwards=TRUE),AIC=which.min(aics),BIC=which.min(bics))
## dhatsTot <- c(dhatAll(pvalsFull,backwards=TRUE),AIC=which.min(aics),BIC=which.min(bics))
