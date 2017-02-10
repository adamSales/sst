
library(vioplot)
library(parallel)
cl<-makePSOCKcluster(min(10,detectCores()))

source('src/winSelect.r')

beta <- 0.1
ns <- c(10,50,100)

#### simulation 1
#### one covariate x
#### r=+/- 1:30
#### compare x for negative r to positive r within bandwidth bw=1:30
#### for r= +/- 1:10 x~N(0,1)
#### for |r|>10 x=-10b+br+N(0,1) for b=beta defined above in code
#### if curved, for |r|>20 x=20b-br+N(0,1) so that relationship reverses itself at r=20
#### vary sample size at each point n and whether curved
#### sst techniques: two factors
#### first: edge testing, on or off
#### second: how to choose based on p-values
#### method 1: biggest d st p> alpha, alpha=0.05, 0.15, 0.25
#### method 2: biggest d st p[d]> alpha & p[d+1]<alpha
#### method 3: mallik
#### method 4: flexible mallik


r <- seq(-30,30)
r <- r[r!=0] ## zero just confuses things
bws <- sort(unique(abs(r)))

makeData <- function(n,curved,beta=0.1){
    if(curved){
        x <- vapply(r,function(r)
            rnorm(n)+ifelse(abs(r)>10,
                            sin((r-10*sign(r))*pi/10),0),numeric(n))
    } else
        x <- vapply(r,function(r) rnorm(n)+ifelse(abs(r)>10,(r-sign(r)*10)*beta,0),numeric(n))
    x
}

test <- function(x,bw,edge=FALSE){
    if(edge)
        return(t.test(x[,r==bw],x[,r==-bw])$p.value)
    t.test(x[,r<=bw & r>0],x[,r>= -bw & r<0])$p.value
}

simOneOne <- function(x,alphas,edge){
    ps <- vapply(bws,function(bw) test(x=x,bw=bw,edge=edge),1)
    bstar <- setNames(c(vapply(alphas,function(alph) dhatAlpha(ps,alph),1),
                            vapply(alphas,function(alph) dhatAlphaBad(ps,alph),1),
                            mallik=dhatM(ps),mallikFlex=dhatM2(ps)['dhat']),
                          c(paste0('\\hat{d}^{max}_{',alphas,'}'),
                            paste0('\\hat{d}^{min}_{',alphas,'}'),
                            '\\hat{d}_M','\\hat{d}^{flex}_M'))
    bstar
}

simOne <- function(n,curved,beta=0.1,alphas=c(0.05,0.15,0.25)){
    x <- makeData(n,curved, beta)
    list(#edge=simOneOne(x=x,alphas=alphas,edge=TRUE),
         full=simOneOne(x=x,alphas=alphas,edge=FALSE))
}

simTot <- function(nsims,ns =c(10,50,100),beta=beta){

    results <- list()
    stuff <- c(
        ls(envir=.GlobalEnv)[
            vapply(ls(envir=.GlobalEnv), function(obj) class(get(obj,envir=.GlobalEnv))[1],'a')
            =='function'],'r','bws','beta')
    print(stuff)
    clusterExport(cl,stuff)
    for(N in ns){
        print(N)
        results[[paste0('mono_',N)]] <-
            parLapply(cl,1:nsims,function(i) simOne(n=N,curved=FALSE,beta=0.1))
        save(results,file='results.RData')
        print('part 2')
        results[[paste0('curved_',N)]] <-
            parLapply(cl,1:nsims,function(i) simOne(n=N,curved=TRUE,beta=0.1))
        save(results,file='results.RData')
                 }
    stopCluster(cl)
    results
}

fix <- function(st){
    ff <- lapply(st,function(x) do.call('rbind',x))
    for(i in seq(length(st))) attr(ff[[i]],'runName') <- names(st)[i]
    ff
}

summaries <- function(x){
    round(c(mean=mean(x,na.rm=TRUE),
      mse=mean((x-10)^2,na.rm=TRUE),
      msePlus=mean((x[x>10]-10)^2,na.rm=TRUE),
      mseMinus=mean((x[x<10]-10)^2,na.rm=TRUE),
      numNA=sum(is.na(x))),2)
}

violin <- function(run){
    args <- list()
    args[['x']] <- run[!is.na(run[,1]),1]
    for(i in 2:ncol(run))
        args[[i]] <- run[!is.na(run[,i]),i]

    args[['names']] <- colnames(run)
    do.call('vioplot',args)
    abline(h=10,lty=2)
    text(2,5,attr(run,'runName'))
    }



violin2 <- function(run,n,b){
    st1 <- do.call('c',as.data.frame(run))
    st1 <- data.frame(est=st1,method=rep(colnames(run),each=nrow(run)))
    st1$edge <- ifelse(grepl('E',st1$method),'edge','full')
    st1$method <- as.factor(gsub('E','',st1$method))

    p <- ggplot(st1, aes(x=method, y=est,fill=edge)) +
        geom_violin(bw=1,position=position_dodge(.5))+
        geom_boxplot(width=0.1,position=position_dodge(.5))+
        geom_hline(yintercept=10,lty=2)+
        labs(title=paste('n=',n,'b=',b))
    p
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

violinTot <- function(st,ns,bs){
    plots <- list()
    for(n in ns) for(b in bs) plots[[paste0(n,'_',b)]] <- violin2(st[[paste0(n,'_',b)]],n,b)
    multiplot(plotlist=plots,cols=3)
}


aaa <- function(st){


    plotlist <- list()
    for(curv in c('mono','curved'))
        for(n in c(10,50,100)){
            run <- do.call('rbind',lapply(st[[paste0(curv,'_',n)]],function(x) x$full))
            plotlist[[paste0(curv,n)]] <- violin3(run,curv,n)
        }

    multiplot(plotlist=plotlist,  cols=2)

}



violin3 <- function(run,n,curv){
    st1 <- do.call('c',as.data.frame(run))
    st1 <- data.frame(est=st1,method=rep(colnames(run),each=nrow(run)))
    st1$selector <- ifelse(grepl('max',st1$method),'$\\bar{d}_\\alpha$',
                       ifelse(grepl('min',st1$method),'$\\underline{d}_\\alpha$','\\hat{d}_M'))
    st1$method <- gsub('\\\\\\hat\\{d\\}\\^\\{m[axin]*\\}_\\{','$\\\\\\alpha=',st1$method)
    st1$method <- gsub('\\}','$',st1$method)
    st1$method <- gsub('\\\\\\hat\\{d\\$\\^\\{flex\\$_M','a,b',st1$method)
    st1$method <- gsub('\\\\\\hat\\{d\\$_M','0.5,1',st1$method)

    st1$method <- factor(st1$method)


    p <- ggplot(st1, aes(x=selector, y=est,fill=method)) +
        geom_violin(bw=1,position=position_dodge(.5))+
        geom_boxplot(width=0.1,position=position_dodge(.5))+
        geom_hline(yintercept=10,lty=2)+
            labs(title=paste('n=',n,curv),x='$\\hat{d}$')
    #if(curv!='curved' | n!=50) p <- p+guides(fill=FALSE)
    #if(curv & n==50) p <- p+guides(fill=TRUE)
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


simTable <- function(st,curv){
    ss <- sumStats(st)

    ss <- lapply(ss, function(x) x[,-grep('mse',colnames(ss))])

    tab <- do.call("cbind",ss[ifelse(rep(curv,3),4:6,1:3)])
    rownames(tab) <- paste0('$',c("\\dhatU_{0.05}","\\dhatU_{0.15}","\\dhatU_{0.25}","\\dhatB_{0.05}","\\dhatB_{0.15}","\\dhatB_{0.25}","\\dhatm","\\dhatmab"),'$')

    colnames(tab) <- rep(c('avg','sd','$\\le 3$','$>13$'),3)

    xtab <- xtable(tab)
    align(xtab) <- "r|cccc|cccc|cccc|"
    addtorow <- list()
    addtorow$pos <- list(-1)
    addtorow$command <- c("& \\multicolumn{4}{c}{$n=10$}&\\multicolumn{4}{c}{$n=50$} &\\multicolumn{4}{c}{$n=100$}\\\\\n")
    print(xtab,sanitize.text.function=function(x) x,add.to.row=addtorow)
}
