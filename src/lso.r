require(RItools)
ap <- read.csv('~/Dropbox/Sales Hansen/drafts/dataResults/LindoDat.csv')
ap$r <- round(ap$dist_from_cut,2)
ap <- ap[abs(ap$r)<=0.5,]
ap$r <- round(ap$r+0.005,3)

ap$z <- ap$r<0

bws <- sort(unique(abs(ap$r)))

balform <- z~hsgrade_pct+age_at_entry+totcredits_year1+male+loc_campus1+loc_campus2+english
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

balTestTW <- function(bw){
    datbw <- ap[abs(ap$r)<=bw,]
    with(datbw,wilcox.test(hsgrade_pct[z],hsgrade_pct[!z]))$p.value
}

balTestEW <- function(bw){
    datbw <- ap[abs(ap$r)==bw,]
    res <- try(with(datbw,wilcox.test(hsgrade_pct[z],hsgrade_pct[!z]))$p.value)
    if(class(res)=='try-error') return(NA)
    res

}


psTW <- vapply(bws,balTestTW,1)
psEW <- vapply(bws,balTestEW,1)
