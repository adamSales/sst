require(RItools)
require(foreign)
source('src/winSelect.r')
ap <- read.dta('data/AEJApp2008-0202_data/data_for_analysis.dta')
#### downloaded 1/2/17 from https://www.aeaweb.org/aej/app/data/2008-0202_data.zip

ap$r <- round(ap$dist_from_cut,2)
ap <- ap[abs(ap$r)<=0.5,]
ap$r <- round(ap$r+0.005,3)

ap$z <- ap$r<0

bws <- sort(unique(abs(ap$r)))

balform <- z~hsgrade_pct+age_at_entry+totcredits_year1+male+loc_campus1+loc_campus2+english+
    bpl_north_america

balTestT <- function(bw){
    xBalance(balform,data=ap[abs(ap$r)<=bw,], report='chisquare.test')$overall$p.value[1]
}
balTestE <- function(bw){
    res <- try(xBalance(balform,data=ap[abs(ap$r)==bw,],
                        report='chisquare.test')$overall$p.value[1])
    if(class(res)=='try-error') return(NA)
    res
}

psT <- vapply(bws,balTestT,1)
psE <- vapply(bws,balTestE,1)

windowsT <- dhatAll(psT)
windowsE <- dhatAll(psE[-51]) ## the last value is NA

balTestT1 <- function(bw){
    datbw <- ap[abs(ap$r)<=bw,]
    with(datbw,wilcox.test(hsgrade_pct[z],hsgrade_pct[!z]))$p.value
}

balTestE1 <- function(bw){
    datbw <- ap[abs(ap$r)==bw,]
    res <- try(with(datbw,wilcox.test(hsgrade_pct[z],hsgrade_pct[!z]))$p.value)
    if(class(res)=='try-error') return(NA)
    res

}


psT1 <- vapply(bws,balTestT1,1)
psE1 <- vapply(bws,balTestE1,1)

windowsT1 <- dhatAll(psT1)
windowsE1 <- dhatAll(psE1[-51])

dhatm2 <- dhatM2(psE1[-51])


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
