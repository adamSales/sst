library(parallel)
cl<-makePSOCKcluster(10)





#### simulation 1
#### one covariate x
#### r=+/- 1:20
#### compare x for negative r to positive r within bandwidth bw=1:20
#### for r= +/- 1:10 x~N(0,1)
#### for r= +/- 11:20 x=-10b+br+N(0,1) for varying b
#### also vary sample size at each point n
#### sst techniques: two factors
#### first: edge testing, on or off
#### second: how to choose based on p-values
#### method 1: biggest p> alpha, alpha=0.05, 0.15, 0.25
#### method 2: mallik
#### method 3: ks criterion

r <- seq(-20,20)

makeData <- function(n,beta){
    x <- vapply(r,function(r) rnorm(n)+ifelse(abs(r)>10,(r-sign(r)*10)*beta,0),numeric(n))
    x
}

test <- function(x,bw,edge=TRUE){
    if(edge)
        return(wilcox.test(x[,r==bw],x[,r==-bw])$p.value)
    wilcox.test(x[,r<=bw & r>0],x[,r>= -bw & r<0])$p.value
}

alphaTest <- function(ps,alpha=c(0.05,0.15,0.25)){
    bstars <- NULL
    for(a in alpha)
        bstars <- c(bstars,seq(2,20)[max(which(ps>a))])
    names(bstars) <- alpha
    bstars
}

mallik <- function(ps){
    seq(2,20)[which.max(cumsum(ps-1/4))]
}

ks <- function(ps){
    cdf <- function(x) x
    stat <- function(bw){
        ks.test(ps[1:bw],cdf,alternative='greater')$statistic+
            ks.test(ps[(bw:1):length(ps)],cdf,alternative='less')$statistic
    }
    seq(2,20)[which.min(vapply(2:length(ps),stat,1))]
}

simOne <- function(n,beta){
    x <- makeData(n,beta)
    psEdge <- vapply(2:20,function(bw) test(x=x,bw=bw,edge=TRUE),1)
    bstar <- c(alphaE=alphaTest(psEdge),mallikE=mallik(psEdge),ksE=ks(psEdge))

    psFull <- vapply(2:20,function(bw) test(x=x,bw=bw,edge=FALSE),1)
    bstar <- c(bstar,alpha=alphaTest(psFull),mallik=mallik(psFull),ksFull=ks(psFull))
    bstar
}

simTot <- function(nsims,ns =c(10,50,100),betas=c(0.2,0.5,1)){

    results <- list()
    stuff <- c(ls(envir=.GlobalEnv)[vapply(ls(envir=.GlobalEnv), function(obj) class(get(obj))[1],'a')=='function'],'r')
    print(stuff)
    clusterExport(cl,stuff)
    for(n in ns) for(beta in betas){
        results[[paste0(n,'_',beta)]] <- parLapply(cl,1:nsims,function(i) simOne(n,beta))
    }
    stopCluster(cl)
    results
}

