#### simulating Pr(Max(sum(p-1/4)_[1:n])=n) when all p~U(0,1)
sim1 <- function(r){
    uu <- runif(r,-1/4,3/4)
    cs <- cumsum(uu)
    vapply(1:r,function(a) which.max(cs[1:a]),integer(1))
}
sims <- vapply(1:100000, function(i) sim1(100),integer(100))


Mode <- function(x) {
     ux <- unique(x)
     ux[which.max(tabulate(match(x, ux)))]
 }


#plot(apply(sims,1,Mode))

prbRigth <- vapply(1:100,function(x) mean(sims[x,]==x),1)

