
dhatAlpha <- function(ps,alpha){
    D <- length(ps)
    res <- try(max(seq(D)[ps>=alpha],na.rm=TRUE))
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

dhatAll <- function(ps, alphas=c(0.05,0.25),backwards=FALSE,m2=TRUE){
    if(backwards) ps <- ps[length(ps):1]
    res <- setNames(
        c(vapply(alphas, function(alph) dhatAlpha(ps,alph),1),
          vapply(alphas, function(alph) dhatAlphaBad(ps,alph),1),
          dhatM(ps)),
#          dhatM2(ps)[2]),
        c(paste0('$\\dhatU_{',alphas,'}$'),
          paste0('$\\dhatB_{',alphas,'}$'),
          '$\\dhatm$'))#,'$\\dhatmab$')

    if(m2) res <- setNames(c(res,dhatM2(ps)[2]),c(names(res),'$\\dhatmab$'))

    if(backwards) res <- vapply(res,function(x) length(ps)-x+1,1)
    res
}

ds <- function(ps,alphas=c(0.05,0.15,0.25)){
    dhat <- list()
    for(a in alphas){
        dhat[[paste0('alpha',a)]] <- dhatAlpha(ps,a)
    }
    dhat[['Mallik']] <- dhatM(ps)
    dhat
}


multLine <- function(dhats,amount=0.1){
    ttt <- table(dhats)
    if(max(ttt)==1) return(dhats)
    for(d in unique(dhats)){
        num <- ttt[paste(d)]
        if(num>1)
            dhats[dhats==d] <- dhats[dhats==d]+seq(-floor(num/2)*amount,by=amount,length=num)
    }
    dhats
}

cols <- colors()[c(red=552,lightblue=636,mag=642,green=255,orange=499,black=153,purple=469,darkblue=477)]
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


smoothplot <- function(x,y,...){


    naFlag <- is.na(x) | is.na(y)
    x <- x[!naFlag]
    y <- y[!naFlag]


    tmp <- vapply(unique(x), function(xval) c(
        av=mean(y[x==xval],na.rm=TRUE),
        sizes = length(y[x==xval])
        ),numeric(2))

    tmp['sizes',] <- 2*tmp['sizes',]/max(tmp['sizes',])
    plot(unique(x),tmp['av',],cex=tmp['sizes',],...)
    abline(v=0,lty='dotted')


}


facet2 <- function(st,curv){
    ns <- sort(unique(vapply(strsplit(names(st),'_'),function(x) as.numeric(x[2]),1)))
    dat <- NULL
    cname <- ifelse(curv,'curved','mono')
    for(n in ns){

        run <- do.call('rbind',lapply(st[[paste0(cname,'_',n)]],function(x) x$full))
        st1 <- do.call('c',as.data.frame(run))
        st1 <- data.frame(est=st1,method=rep(colnames(run),each=nrow(run)))
        st1$selector <- ifelse(grepl('max',st1$method),'$\\bar{d}_\\alpha$',
                               ifelse(grepl('min',st1$method),'$\\tilde{d}_\\alpha$','$\\hat{d}_M$'))
        st1$method <- gsub('\\\\\\hat\\{d\\}\\^\\{m[axin]*\\}_\\{','$\\\\\\alpha=',st1$method)
        st1$method <- gsub('\\}','$',st1$method)
        st1$method <- gsub('\\\\\\hat\\{d\\$\\^\\{flex\\$_M','a,b',st1$method)
        st1$method <- gsub('\\\\\\hat\\{d\\$_M','0.5,1',st1$method)

        st1$method <- factor(st1$method)
        st1$n <- n

        dat <- rbind(dat,st1)
    }

    dat$n <- factor(dat$n)
    levels(dat$n)=paste0('n=',levels(dat$n))

    p <- ggplot(dat, aes(x=selector, y=est,fill=method)) +
        geom_violin(bw=1,position=position_dodge(.5))+
        geom_boxplot(width=0.1,position=position_dodge(.5),outlier.size=.5)+
        geom_hline(yintercept=10,lty=2)+
            labs(title=paste('Imbalance:',ifelse(curv,'Sinusoidal','Linear')),y='$\\hat{d}$',x='')+
                facet_wrap(~ n,ncol=3)
    p
}

sumStats <- function(st){
    plotlist <- list()
    for(curv in c('mono','curved'))
        for(n in c(10,50,100)){
            run <- do.call('rbind',lapply(st[[paste0(curv,'_',n)]],function(x) x$full))
            plotlist[[paste0(curv,n)]] <- cbind(
                lessThan3=apply(run,2,function(x) mean(x[is.finite(x)]<=3,na.rm=TRUE)),
                moreThan13=apply(run,2,function(x) mean(x[is.finite(x)]>13,na.rm=TRUE)),
                avg=apply(run,2,function(x) mean(x[is.finite(x)],na.rm=TRUE)),
                mse=apply(run,2,function(x) mean((x[is.finite(x)]-10)^2,na.rm=TRUE)),
                sd=apply(run,2,function(x) sd(x[is.finite(x)],na.rm=TRUE)))

        }
    plotlist
}
