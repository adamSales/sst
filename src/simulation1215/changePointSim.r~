library(vioplot)
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
    bstar <- c(alphaE=alphaTest(psEdge),mallikE=mallik(psEdge))

    psFull <- vapply(2:20,function(bw) test(x=x,bw=bw,edge=FALSE),1)
    bstar <- c(bstar,alpha=alphaTest(psFull),mallik=mallik(psFull))
    bstar
}

simTot <- function(nsims,ns =c(10,50),betas=c(0.05,0.1,.2)){

    results <- list()
    stuff <- c(ls(envir=.GlobalEnv)[vapply(ls(envir=.GlobalEnv), function(obj) class(get(obj))[1],'a')=='function'],'r')
    print(stuff)
    clusterExport(cl,stuff)
    for(n in ns) for(beta in betas){
        results[[paste0(n,'_',beta)]] <- parLapply(cl,1:nsims,function(i) simOne(n,beta))
        save(results,file='results.RData')
                 }
    stopCluster(cl)
    results
}

fix <- function(st)
    lapply(st,function(x) do.call('rbind',x))

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
    multiplot(plotlist=plots,cols=2)
}
