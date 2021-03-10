library(tidyverse)
### simulate test statistics directly

## b is the number of replications per step
## np is the number of steps
## cp is the change point---last null model
## slope is self explanatory


sim1 <- function(b,cp,slope,np){
    n <- b*np
    T <- matrix(rnorm(n),nrow=b)
    ncp <- c(rep(0,cp),slope*seq(np-cp))
    T <- sweep(T,2,ncp,'+')
    2*pnorm(-abs(T))
    }

sim <- function(b){
    set.seed(613)
    conditions <- expand.grid(slope=seq(.5,3,.5),cp=c(2,5,10),np=c(10,20))
    conditions <- subset(conditions,cp<np)
    lapply(1:nrow(conditions),
           function(i)
               list(
                   cond=conditions[i,],
                   ps= with(conditions,
                            sim1(b=b,cp=cp[i],slope=slope[i],np=np[i])
                            )
               )
           )
}

dM <- function(ps) which.max(cumsum(ps-1/4))
dalpha <- function(ps,alpha=0.05) max(which(ps>alpha))
dmin <- function(ps,alpha=0.05) min(which(ps<alpha))-1

proc1 <- function(res1,alt,alpha=0.05){
    dsM <- apply(res1[[2]],1,dM)
    dsAlph <- apply(res1[[2]],1,alt,alpha=alpha)
    d <- seq(res1[[1]]$np)
    out <- map_dfr(d,~tibble(d=.,m=mean(dsM==.),alph=-mean(dsAlph==d)))
    out <- cbind(rbind(res1[[1]]),out)%>%
        pivot_longer(cols=5:6,names_to="Estimator",values_to="dstar")
}

proc <- function(res)
    map_dfr(c(0.05,0.25),function(a)
        bind_rows(
            map_dfr(res,proc1,alt=dalpha,alpha=a)%>%mutate(alt=paste0("dmax",a)),
            map_dfr(res,proc1,alt=dmin,alpha=a)%>%mutate(alt=paste0("dmin",a))
        )
        )

plotRes <- function(res,NP){
    proc(res)%>%
        filter(np==NP,d<=NP)%>%
        ggplot(aes(d,dstar,fill=Estimator))+#,xintercept=cp))+
        geom_col()+
        geom_vline(aes(xintercept=cp))+
        scale_x_continuous(breaks=1:NP)+
        scale_y_continuous(breaks=c(-.75,-.5,-.25,0,.25,.5,.75),labels=c("75","50","25","0","25","50","75"))+
        coord_flip()+
        facet_grid(cp+alt~slope)
}
