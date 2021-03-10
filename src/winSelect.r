dhatAlpha <- function(ps,alpha){
    D <- length(ps)
    res <- try(max(seq(D)[ps>=alpha]))
    if(class(res)=='try-error') return(NA)
    res
}


dhatM <- function(ps){
    M <- cumsum(ps-0.25)
    res <- which.max(M)
    if(class(res)=='try-error') return(NA)
    res
}

m2 <- function(ps,a,b,d,D)
    sum((ps[1:d]-a)^2)+sum((ps[(d+1):D]-b)^2)

M2 <- function(ps,a,b,D){
    mm <- vapply(seq(D),function(d) m2(ps,a,b,d,D),1)
    c(min(mm),which.min(mm),a,b)
}

dhatM2 <- function(ps,stepsize=0.01){
    D <- length(ps)-1 ### the last value shouldn't be computed--will give NA

    ## minimum p-value given a, over b. b<a, a>=0.5
    mina <- function(a){
        mma <- vapply(seq(0,a-stepsize,by=stepsize),function(b) M2(ps,a,b,D),numeric(4))
        mma[,which.min(mma[1,])]
    }

    mmb <- vapply(seq(0.5,1,by=stepsize),mina,numeric(4))
    setNames(mmb[,which.min(mmb[1,])],c('Mopt','dhat','aopt','bopt'))
}



dhatAlphaBad <- function(ps,alpha){
    D <- length(ps)
    res <- try(min(seq(D)[ps<alpha])-1)
    if(class(res)=='try-error') return(NA)
    res
}

dhatAll <- function(ps, alphas=c(0.05,0.25),backwards=FALSE){
    if(backwards) ps <- ps[length(ps):1]
    res <- setNames(
        c(vapply(alphas, function(alph) dhatAlpha(ps,alph),1),
          vapply(alphas, function(alph) dhatAlphaBad(ps,alph),1),
          dhatM(ps),
          dhatM2(ps)[2]),
        c(paste0('$\\bar{d}_{',alphas,'}$'),
          paste0('$\\tilde{d}_{',alphas,'}$'),
          '$\\dhatm$','$\\dhatmab$')
        )
    if(backwards) res <- vapply(res,function(x) length(ps)-x+1,1)
    res
}

