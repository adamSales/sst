library(tidyverse)
library(tidyLPA)
library(parallel)
cl<-makePSOCKcluster(min(4,detectCores()))

simulateData <- function(n=200,mu1=1){
    nclass <- 4
    p <- 16
    cls <- sample(1:nclass,n,replace=TRUE)
    dat <- matrix(rnorm(n*p),nrow=n)

    means <- matrix(0,nrow=n,ncol=p)

    for(cc in 1:4)
        means[cls==cc,seq(4*(cc-1)+1,4*cc)] <- mu1

    dat <- dat+means
    as.data.frame(dat)
}

simOne <- function(i,n=200,mu1=1){
    require(tidyLPA)
    dat <- simulateData(n)
    mods <- estimate_profiles(dat,1:10)
    save(mods,file=paste0('simResults/fullMods',mu1,'_',i,'.RData'))
    compare_solutions(mods)$fit
}

clusterExport(cl,c('simulateData','simOne'))

print(system.time(res05 <- parLapply(cl,1:500,simOne,mu1=0.5)))

save(res,file='simResults/res.RData')

stopCluster(cl)
