
#################################################################

# Auxiliary functions for rdlocrand

# !version 0.01 27-Aug-2016

# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare

#################################################################



#################################################################
# rdrandinf observed statistics and asymptotic p-values
#################################################################

rdrandinf.model = function(Y,D,statistic,pvalue=FALSE,kweights,endogtr,delta=''){

  n = length(D)
  n1 = sum(D)
  n0 = n - n1

  Y = as.matrix(Y)

  if (statistic=='ttest'){
    if (all.equal(kweights,rep(1,n))==TRUE){
      Y1 = t(Y[D==1,])
      Y0 = t(Y[D==0,])
      M1 = rowMeans(Y1)
      M0 = rowMeans(Y0)
      stat = M1-M0
      if (pvalue==TRUE){
        V1 = rowMeans((Y1-rowMeans(Y1))^2)/(n1-1)
        V0 = rowMeans((Y0-rowMeans(Y0))^2)/(n0-1)
        se = sqrt(V1+V0)
        t.stat = (M1-M0)/se
        asy.pval = 2*pnorm(-abs(t.stat))
        if (delta!=''){
          asy.power = 1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se)
        } else{asy.power = NA}
      }
    } else {
      stat = NULL
      asy.pval = NULL
      for (k in 1:ncol(Y)){
        lm.aux = lm(Y[,k] ~ D,weights=kweights)
        stat = c(stat,lm.aux$coefficient['D'])
        if (pvalue==TRUE){
          se = sqrt(vcovHC(lm.aux,type='HC2')['D','D'])
          t.stat = stat/se
          asy.pval = c(asy.pval,2*pnorm(-abs(t.stat)))
          if (delta!=''){
            asy.power = 1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se)
          } else {asy.power = NA}
        }
      }
    }
  }

  if (statistic=='ksmirnov'){
    stat = NULL
    asy.pval = NULL
    for (k in 1:ncol(Y)){
      aux.ks = ks.test(Y[D==0,k],Y[D==1,k])
      stat = c(stat,aux.ks$statistic)
      if (pvalue==TRUE){
        asy.pval = c(asy.pval,aux.ks$p.value)
        asy.power = NA
      }
    }
  }

  if (statistic=='ranksum'){
    stat = NULL
    asy.pval = NULL
    for (k in 1:ncol(Y)){
      Ustat = wilcox.test(Y[D==0,k],Y[D==1,k])$statistic
      Tstat = Ustat + n0*(n0+1)/2
      ri = rank(Y[,k])
      s2 = var(ri)
      ETstat = n0*(n+1)/2
      VTstat = n0*n1*s2/n
      se = sqrt(VTstat)
      stat = c(stat,(Tstat-ETstat)/se)
    }
    if (pvalue==TRUE){
      asy.pval = 2*pnorm(-abs(stat))
      sigma = sd(Y)
      if (delta!=''){
        asy.power = pnorm(sqrt(3*n0*n1/((n0+n1+1)*pi))*delta/sigma-1.96)
      } else {asy.power = NA}
    }
  }

  if (statistic=='all'){
    stat1 = mean(Y[D==1])-mean(Y[D==0])
    aux.ks = ks.test(Y[D==0],Y[D==1])
    stat2 = aux.ks$statistic
    Ustat = wilcox.test(Y[D==0],Y[D==1])$statistic
    Tstat = Ustat + n0*(n0+1)/2
    ri = seq(1,n)
    s2 = var(ri)
    ETstat = n0*(n+1)/2
    VTstat = n0*n1*s2/n
    se3 = sqrt(VTstat)
    stat3 = (Tstat-ETstat)/se3
    stat = c(stat1,stat2,stat3)
    if (pvalue==TRUE){
      se1 = sqrt(var(Y[D==1])/n1+var(Y[D==0])/n0)
      t.stat = stat1/se1
      asy.pval1 = 2*pnorm(-abs(t.stat))
      asy.pval2 = aux.ks$p.value
      asy.pval3 = 2*pnorm(-abs(stat3))
      asy.pval = c(asy.pval1,asy.pval2,asy.pval3)
      if (delta!=''){
        asy.power1 = 1-pnorm(1.96-delta/se1)+pnorm(-1.96-delta/se1)
        asy.power2 = NA
        sigma = sd(Y)
        asy.power3 = pnorm(sqrt(3*n0*n1/((n0+n1+1)*pi))*delta/sigma-1.96)
        asy.power = c(asy.power1,asy.power2,asy.power3)
      } else {asy.power = c(NA,NA,NA)}
    }
  }

  if (statistic=='ar'){
    stat = mean(Y[D==1])-mean(Y[D==0])
    if (pvalue==TRUE){
      se = sqrt(var(Y[D==1])/n1+var(Y[D==0])/n0)
      t.stat = stat/se
      asy.pval = 2*pnorm(-abs(t.stat))
      if (delta!=''){
        asy.power = 1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se)
      } else {asy.power = NA}
    }
  }

  if (statistic=='wald'){
    fs = lm(endogtr ~ D)
    rf = lm(Y ~ D)
    stat = rf$coefficients['D']/fs$coefficients['D']
    if (pvalue==TRUE){
      ehat = Y - mean(Y) - stat*(endogtr - mean(endogtr))
      ehat2 = ehat^2
      se = sqrt((mean(ehat2)*var(D))/(n*cov(D,endogtr)^2))
      tstat = stat/se
      asy.pval = 2*pnorm(-abs(tstat))
      if (delta!=''){
        asy.power = 1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se)
      } else {asy.power = NA}
    }
  }


  if (pvalue==TRUE){
    output = list(statistic = as.numeric(stat),p.value = as.numeric(asy.pval),asy.power = as.numeric(asy.power))
  } else{
      output = list(statistic = as.numeric(stat))
  }

  return(output)
}


#################################################################
# Hotelling's T2 statistic
#################################################################

hotelT2 = function(X,D) {

  n = length(D)
  n1 = sum(D)
  n0 = n - n1
  p = ncol(X)

  X1 = X[D==1,]
  X0 = X[D==0,]
  X1bar = as.matrix(colMeans(X[D==1,]))
  X0bar = as.matrix(colMeans(X[D==0,]))
  S1 = cov(X1)
  S0 = cov(X0)
  Spool = (S1*(n1-1)+S0*(n0-1))/(n-2)
  SpoolInv = solve(Spool)

  T2 = (n0*n1/n)*t(X1bar-X0bar)%*%SpoolInv%*%(X1bar-X0bar)
  Fstat = ((n-p-1)/((n-2)*p))*T2
  pval = 1-pf(Fstat,p,n-1-p)

  output = list(statistic = as.numeric(T2),Fstat = as.numeric(Fstat),p.value = as.numeric(pval))
  names(output) = c('statistic','Fstat','p.value')

  return(output)
}


Mabd <- function(ps,a,b,d)
    sum((ps[1:d]-a)^2)+sum((ps[(d+1):length(ps)]-b)^2)

dhatmg <- function(ps,inc=0.01){
    as <- seq(0,1,by=inc)
    np <- length(as)
    params <- c(0,0,0)
    theMin <- Mabd(ps,1,0,1)
    for(a in as[5:np]){
        for(b in as[as<a]){
            for(d in 1:(length(ps)-1)){
                M <- Mabd(ps,a,b,d)
                if(M<=theMin){
                    theMin <- M
                    print(params <- c(a,b,d))
                }
            }
        }
    }
    c(params,M)
}

### for edge testing w/ smoothing
### test H_0r: E[x|R=-r]=E[x|R=r]
### for vector valued X
### use Nadaraya–Watson estimator for E[x|R=r]

smoothCov <- function(X,R,D,h){

    ## get resids: x-mu(x)
    ## mu(x)= sum( x K((R-r)/h))/sum(K((R-r)/h))
    resids <- NULL
    for(trt in c(0,1)){
        Xs <- X[D==trt,]
        Rs <- R[D==trt]
        mu <- matrix(nrow=nrow(Xs),ncol=ncol(Xs))

        for(i in 1:nrow(Xs)){
            kern <- dnorm((Rs[i]-Rs)/h)
            fr <- sum(kern)
            for(j in 1:ncol(Xs)) mu[i,j] <- sum(Xs[,j]*kern,na.rm=TRUE)/fr
        }
        resids <- rbind(resids,Xs-mu)
    }
    cov(resids,use='pairwise.complete.obs')
}

smoothCov1 <- function(X,R,D,h){

    ## get resids: x-mu(x)
    ## mu(x)= sum( x K((R-r)/h))/sum(K((R-r)/h))
    mu <- matrix(nrow=nrow(X),ncol=ncol(X))

    for(i in 1:nrow(X)){
        kern <- dnorm((R[i]-R)/h)
        fr <- sum(kern)
        for(j in 1:ncol(X)) mu[i,j] <- sum(X[,j]*kern,na.rm=TRUE)/fr
    }
    resids <- X-mu

    cov(resids,use='pairwise.complete.obs')
}

xbalsmooth <- function(X,D,R,r,h,smCv){
    require(Matrix)
    n = length(D)
    n1 = sum(D)
    n0 = n - n1
    p = ncol(X)



    R1 <- R[D==1]
    R0 <- R[D==0]

    kern1 <- dnorm((r-R1)/h)
    kern0 <- dnorm((-r-R0)/h)
    n0 <- sum(kern0)
    n1 <- sum(kern1)
    n <- n0+n1
    hm <- (1/n1+1/n0)^(-1)

    X1 = X[D==1,]*kern1
    X0 = X[D==0,]*kern0


    X1bar = colSums(X1,na.rm=TRUE)/sum(kern1)
    X0bar = colSums(X0,na.rm=TRUE)/sum(kern0)

    if(missing(smCv)){
        S = smoothCov(X,R,D,h)/hm
    } else{
        S <- smCv/hm
    }

    Sinv <- ginv(S)
    d <- X1bar-X0bar
    T2 = rbind(d)%*%Sinv%*%cbind(d)

    pval = pchisq(T2,df=rankMatrix(S),lower.tail=FALSE)

    output = list(statistic = as.numeric(T2),p.value = as.numeric(pval))
    names(output) = c('statistic','p.value')

  return(output)
}

hotelT2smooth <- function(X,D,R,r,h,smCv){

    n = length(D)
    n1 = sum(D)
    n0 = n - n1
    p = ncol(X)

    R1 <- R[D==1]
    R0 <- R[D==0]

    kern1 <- dnorm((r-R1)/h)
    kern0 <- dnorm((r-R0)/h)
    n0 <- sum(kern0)
    n1 <- sum(kern1)
    n <- n0+n1

    X1 = X[D==1,]*kern1
    X0 = X[D==0,]*kern0


    X1bar = as.matrix(colSums(X1,na.rm=TRUE)/sum(kern1))
    X0bar = as.matrix(colSums(X0,na.rm=TRUE)/sum(kern0))

    if(missing(smCv)){
        S1 = smoothCov(X1,R1,h)
        S0 <- smoothCov(X0,R0,h)
    } else{
        S1 <- smCv[[1]]
        S0 <- smCv[[2]]
    }

    Spool = (S1*(n1-1)+S0*(n0-1))/(n-2)
    SpoolInv = solve(Spool)

    T2 = (n0*n1/n)*t(X1bar-X0bar)%*%SpoolInv%*%(X1bar-X0bar)
    Fstat = ((n-p-1)/((n-2)*p))*T2
    pval = 1-pf(Fstat,p,n-1-p)

    output = list(statistic = as.numeric(T2),Fstat = as.numeric(Fstat),p.value = as.numeric(pval))
  names(output) = c('statistic','Fstat','p.value')

  return(output)
}


xbal <- function(X,D){
    dat1 <- data.frame(D=D,X)
    mod <- xBalance(D~.,data=dat1,report='chisquare.test',na.rm=TRUE)
    mod$overall$p
}

#################################################################
# Find default window length in rdwinselect
#################################################################

wlength = function(R,D,num){

  X1 = R[D==1]
  X1 = sort(abs(X1))
  X0 = R[D==0]
  X0 = sort(abs(X0))
  m = min(length(X1),length(X0))
  if (num>m){
    num = m
  }
  xt = X1[num]
  xc = X0[num]
  minw = max(xc,xt)
  return(minw)
}


#################################################################
# Find default step in rdwinselect
#################################################################

findstep = function(R,D,obsmin,obsstep,times) {
  S = NULL
  for (i in 1:times){
    U = wlength(R,D,obsmin+obsstep*i)
    L = wlength(R,D,obsmin+obsstep*(i-1))
    Snext = U - L
    S = c(S,Snext)
  }
  step = max(S)
  return(step)
}

