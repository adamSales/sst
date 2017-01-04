fix <- function(st){
    ff <- lapply(st,function(x) do.call('rbind',lapply(x,fixRun)))
    for(i in seq(length(st))) attr(ff[[i]],'runName') <- names(st)[i]
    ff
}

fixRun <- function(run){
    alphE <- run$edge[1:3]
    alphT <- run$full[1:3]
#    names(alphE) <- sub('max','edge',names(alphE))
#    names(alphT) <- sub('max','tot',names(alphT))
    m1e <- run$edge[7]
    m1t <- run$full[7]
    m2e <- run$edge[8]
    m2t <- run$full[8]
    c(alpha05E=alphE[1],alpha05F=alphT[1],alph15E=alphE[2],alph15F=alphT[2],alph25E=alphE[3],alph25F=alphT[3],m1E=m1e,m1F=m1t,m2E=m2e,m2F=m2t)
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



violin2 <- function(run,n){
    st1 <- do.call('c',as.data.frame(run))
    st1 <- data.frame(est=st1,method=rep(colnames(run),each=nrow(run)))
    st1$edge <- ifelse(grepl('E',st1$method),'edge','full')
    st1$method <- as.factor(gsub('E','',st1$method))

    p <- ggplot(st1, aes(x=method, y=est,fill=edge)) +
        geom_violin(bw=1,position=position_dodge(.5))+
        geom_boxplot(width=0.1,position=position_dodge(.5))+
        geom_hline(yintercept=10,lty=2)+
        labs(title=paste('n=',n))
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

violinTot <- function(st,ns){
    plots <- list()
    for(n in ns) for(crv in c(FALSE,TRUE)) plots[[paste0(n,'_',crv)]] <- violin2(st[[paste0(ifelse(crv,'curved_','mono_'),n)]],n)
    multiplot(plotlist=plots,cols=2)
}
