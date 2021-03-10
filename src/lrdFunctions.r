


#################
### SH functions
#################

#' Wrapper function for limitless regression discontinuity testing and estimation
#'
#' In a regression discontinuity design, this function tests for an effect using \code{testSH},
#' and inverts the test to estimate a confidence interval and Hodges-Lehmann point estimate.
#' The test is the test associated with the Z coefficient after fitting
#' a model from the robustbase suite that uses Z and R as independent variables.
#' Covariance adjustment is via robust logistic regression if
#' ytilde is binary, robust linear regression otherwise. In the latter case,
#' since we're using `robustbase::lmrob` with the `MM` method, the SE used
#' in the test is of the Huber-White type, allowing for heteroskedasticity
#' and accounting for  propagation of errors from the preliminary fitting
#' done by the robust regression routine.
#'
#' This is the method recommended in version 2 of Limitless
#' Regression Discontinuity
#'
#' @param dat Data set with columns \code{Z}, \code{R}, an outcome variable, and (optionally)
#' an indicator for treatment received, in the case of partial compliance.
#' @param BW An optional bandwidth \code{b>0}. The test uses data for subjects with |R|<b. When
#' \code{BW} is omitted, the function uses sequential covariate balance tests to choose a
#' data-driven bandwidth.
#' @param outcome A string specifying the name of the outcome variable
#' @param Dvar An optional string specifying the name of the treatment-received variable in the
#' case of partial compliance.
#' @param alpha The level of the confidence interval
#' @param rhs A string specifying the \eqn{\mu_\beta(\cdot)} model relating the outcome to
#' \code{R}. Defaults to \code{"R+Z"}
#'
#' @return a list consisting of
#' \describe{
#'  \item{p.value} The p-value testing no effect
#'  \item{CI} a vector of confidence limits and the HL treatment effect
#'  \item{BW} the RDD bandwidth
#'  \item{bal.pval} a p-value testing for covariate balance
#'  \item{n} the number of subjects in the window of analysis
#'  \item{W} the range of R values in the window of analysis
#' }
#' @import robustbase
#' @export

sh <- function(dat,BW,outcome,Dvar='probation_year1',alpha=0.05,rhs=NULL,...){
    if(missing(BW)) BW <- bwMult(dat,...)
    if(!is.null(Dvar)) dat$D <- dat[[Dvar]]
    bal.pval <- balMult(dat,BW,method='sh',reduced.covars=FALSE,rhs=rhs)
    p.value <- testSH(dat=dat,BW=BW,outcome=outcome,rhs=rhs)
    est <- HLsh(dat=dat,BW=BW,outcome=outcome,rhs=rhs)
    CI <- CIsh(dat=dat,BW=BW,outcome=outcome,est=est,alpha=alpha,rhs=rhs)
    list(p.value=p.value,CI=CI,BW=BW,bal.pval=bal.pval,W=c(min(abs(dat$R)),BW),n=sum(abs(dat$R)<BW))
}

#' Univariate covariance-adjusted balance test, using Huber-White SE
#'
#' The test is the test associated with the Z coefficient after fitting
#' a linear or logistic-linear model with Z and R as independent variables.
#' The Huber-White SE is the one provided by the sandwich package.
#' #'
#'
#' @param dat Data set with columns \code{Z}, \code{R} and a covariate
#' @param BW The bandwidth within which to test balance.
#' @param varb The character name of the variable to test balance on.
#' @param rhs  A string specifying the \eqn{\mu_\beta(\cdot)} model relating the outcome to
#' \code{R}. Defaults to \code{"R+Z"}
#'
#' @return scalar, the p-value associated w/ coefficient on Z
#' @import sandwich
#' @export
balOneSH <- function(dat,BW,varb,rhs=NULL){

    if(sum(dat$R<0 & -BW<dat$R)<2 || sum(dat$R>0 & dat$R<BW)<2) return(NA)
    dat <- dat[abs(dat$R)<BW,]
    dat$ytilde <- dat[[varb]]

    form <- if(is.null(rhs)) as.formula('ytilde~Z+R') else
                   as.formula(paste0('ytilde',rhs))

    if(length(unique(dat$ytilde))>2)
        mod <- lm(form,dat=dat)
    else mod <- glm(form,dat=dat,family=binomial)
    vcv <- sandwich::sandwich(mod)
    Z.pos.c <- pmatch("Z", names(coef(mod)))
    Z.pos.v <- pmatch("Z", colnames(vcv))
    tstat <- coef(mod)[Z.pos.c]/sqrt(vcv[Z.pos.v, Z.pos.v])
    2 * pt(-abs(tstat), df=mod$df.residual, lower.tail=T)
}


#' Robust balance test following robust covariate adjustment
#'
#' The test is the test associated with the Z coefficient after fitting
#' a model from the robustbase suite that uses Z and R as independent variables.
#' Covariance adjustment is via robust logistic regression if
#' ytilde is binary, robust linear regression otherwise. In the latter case,
#' since we're using `robustbase::lmrob` with the `MM` method, the SE used
#' in the test is of the Huber-White type, allowing for heteroskedasticity
#' and accounting for  propagation of errors from the preliminary fitting
#' done by the robust regression routine.
#'
#' This is the method recommended in version 2 of Limitless
#' Regression Discontinuity
#'
#' @param dat Data set with columns \code{Z}, \code{R}, and an outcome variable
#' @param BW A bandwidth \code{b>0}. The test uses data for subjects with |R|<b
#' @param tau A hypothetical treatment effect. Defaults to 0
#' @param outcome A string specifying the name of the outcome variable
#' @param return.coef If TRUE return the Z-coefficient from the model; if FALSE return
#' the p-value
#' @param rhs A string specifying the \eqn{\mu_\beta(\cdot)} model relating the outcome to
#' \code{R}. Defaults to \code{"R+Z"}
#'
#' @return scalar, if return.coef the coefficient on Z; otherwise the p-value associated
#' w/ coefficient on Z
#' @import robustbase
#' @export
testSH <- function(dat,BW,tau=0,outcome='Y',return.coef=FALSE,rhs=NULL){
    if(sum(dat$R<0 & -BW<dat$R)<2 || sum(dat$R>0 & dat$R<BW)<2) return(NA)
    dat <- dat[abs(dat$R)<BW,]

    if('D'%in%names(dat)) dat$ytilde <- dat[[outcome]]-tau*dat$D
    else dat$ytilde <- dat[[outcome]]-tau*dat$Z

    ctl <- if (exists('lmrob_seed')) {
        lmrob.control(seed=lmrob_seed, k.max=500, maxit.scale=500)
    } else lmrob.control(          k.max=500, maxit.scale=500)

    form <- if(is.null(rhs)) as.formula('ytilde~Z+R') else
                   as.formula(paste0('ytilde',rhs))

    mod <- lmrob(form,data=dat,method='MM',control=ctl)

    newcov <- try(LRDsandwich(mod))
    mod$cov <- if (inherits(newcov, "try-error")) NA else newcov
    Z.pos = pmatch("Z", names(coef(mod)))
    if(return.coef) return(coef(mod)[Z.pos])
    coef(summary(mod))[Z.pos,4]

}


#' Hodges-Lehmann estimate of a treatment effect for limitless regression discontinuity
#'
#' Inverts the test in \code{testSH()} to find the Hodges-Lehmann point estimate
#'
#'
#' @param dat Data set with columns \code{Z}, \code{R}, and an outcome variable
#' @param BW A bandwidth \code{b>0}. The test uses data for subjects with |R|<b
#' @param outcome A string specifying the name of the outcome variable
#' @param rhs A string specifying the \eqn{\mu_\beta(\cdot)} model relating the outcome to
#' \code{R}. Defaults to \code{"R+Z"}
#'
#' @return scalar, the Hodges-Lehmann estimate
#' @import robustbase
#' @export

HLsh <- function(dat,BW,outcome='Y',rhs=NULL){

    shortFunction <- function(tau)
        testSH(dat=dat,BW=BW,tau=tau,outcome=outcome,return.coef=TRUE,rhs=rhs)
    out <- try(uniroot(shortFunction,interval=c(-2,2),extendInt='yes')$root)
    ifelse(inherits(out, 'try-error'),NA,out)
}


#' Confidence Interval and Hodges-Lehmann estimate of a treatment effect for limitless
#' regression discontinuity
#'
#' Inverts the test in \code{testSH()} to find the Confidence Interval and Hodges-Lehmann
#' point estimate
#'
#'
#' @param dat Data set with columns \code{Z}, \code{R}, and an outcome variable
#' @param BW A bandwidth \code{b>0}. The test uses data for subjects with |R|<b
#' @param outcome A string specifying the name of the outcome variable
#' @param est A point estimate generated by \code{HLsh}. If not provided it will be estimated
#' automatically.
#' @param alpha the level of the confidence interval.
#' @param rhs A string specifying the \eqn{\mu_\beta(\cdot)} model relating the outcome to \code{R}. Defaults to \code{"R+Z"}
#'
#' @return a vector including
#' \describe{
#'  \item{CI1} and
#'  \item{CI2}, the lower and upper bounds of the confidence interval;
#'  \item{est}, the HL point estimate
#' }
#' @import robustbase
#' @export
CIsh <- function(dat,BW,outcome='Y',est,alpha=0.05,rhs=NULL){

    shortFunction <- function(tau) testSH(dat=dat,BW=BW,tau=tau,outcome=outcome,rhs=rhs)-alpha
    if(missing(est)) est <- HLsh(dat=dat,BW=BW,outcome=outcome)
    CI1 <- uniroot(shortFunction,interval=c(-2,est),extendInt='upX')$root
    CI2 <- uniroot(shortFunction,interval=c(est,2),extendInt='downX')$root
    c(CI1=CI1,CI2=CI2,est=est)
}




###############
### IK method
#############


#' Wrapper for Imbens-Kalyanaraman style RDD analysis
#'
#' Uses OLS to estimate the RDD Local Average Treatment Effect (LATE) as in
#' Imbens, G. and Kalyanaraman, K. (2012). Uses a "rectangular" kernal.
#' Tests balance on a set of pre-treatment covariates using \code{balMult} and either the
#' provided bandwidth or the bandwidth calculated with the IK procedure.
#'
#'
#' @param dat Data set with columns \code{Z}, \code{R}, and an outcome variable
#' @param BW (Optional) A bandwidth \code{b>0}. If not provided it will be estimated from the
#' data
#' @param outcome A string specifying the name of the outcome variable
#'
#' @return a list consisting of
#' \describe{
#'  \item{p.value} The p-value testing no effect
#'  \item{CI} a vector of confidence limits and the HL treatment effect
#'  \item{BW} the RDD bandwidth
#'  \item{bal.pval} a p-value testing for covariate balance
#'  \item{n} the number of subjects in the window of analysis
#'  \item{W} the range of R values in the window of analysis
#' }
#' @import rdd
#' @export

ik <- function(dat,BW=NULL,outcome){
    mod <- ikTest(dat,BW,varb=outcome,justP=FALSE)
    BW <- mod$bw[1]
    bal.pval=balMult(dat=dat,BW=BW,method='sh',reduced.covars=FALSE)
    list(p.value=mod$p[1],CI=c(mod$ci[1,],est=mod$est[1]),bal.pval=bal.pval,
         W=c(min(abs(dat$R)),BW),n=sum(abs(dat$R)<BW))
}


#' Imbens-Kalyanaraman style RDD analysis
#'
#' Uses OLS to estimate the RDD Local Average Treatment Effect (LATE) as in
#' Imbens, G. and Kalyanaraman, K. (2012). Uses a "rectangular" kernal.
#' Tests balance on a set of pre-treatment covariates using \code{balMult} and either the
#' provided bandwidth or the bandwidth calculated with the IK procedure.
#'
#'
#' @param dat Data set with columns \code{Z}, \code{R}, and an outcome variable
#' @param BW (Optional) A bandwidth \code{b>0}. If not provided it will be
#' estimated from the data
#' @param varb A string specifying the name of the outcome or covariate variable
#' @param rhs Ignored
#' @param justP if TRUE only returns the p-value
#'
#' @return if justP=TRUE, returns the p-value for an effect. if justP=FALSE,
#' returns the fitted RDD model (an object of class "RD").
#' @import rdd
#' @export

ikTest <- function(dat,BW=NULL,varb,rhs=NULL,justP=TRUE,cutpoint=-0.005){
    if(!missing(varb)) dat$Y <- -dat[[varb]]
    if(missing(BW) | is.null(BW))
        mod <- try(RDestimate(Y~R,kernel='rectangular',
                              data=dat,cutpoint=cutpoint))
    else  mod <- try(RDestimate(Y~R,kernel='rectangular',
                                data=dat,bw=BW,cutpoint=cutpoint))
    if(class(mod)=='try-error') return(rep(NA,ifelse(justP,1,5)))
    if(justP) return(mod$p[1])
    mod
}

ikMultBal <- function(dat,BW,xvars,int=FALSE){
    require(systemfit)
    xeq <- if(int) lapply(xvars,function(xx) as.formula(paste0(xx,'~R+Z+R:Z'))) else
    lapply(xvars,function(xx) as.formula(paste0(xx,'~R+Z')))
    names(xeq) <- gsub('_','',xvars)

    sur <- systemfit(xeq,data=subset(dat,abs(R)<BW))

    rest <- paste(names(xeq),'Z',sep='_')
    rest <- paste(rest,'=0')

    linearHypothesis(sur,rest,test='Chisq')$Pr[2]
}

###############
### CCT method
#############



#' @export

cct <- function(dat,BW=NULL,outcome){
    mod <- cctTest(dat,BW,varb=outcome,justP=FALSE)
    BW <- mod$bws['h','left']
    bal.pval=balMult(dat=dat,BW=BW,method='sh',reduced.covars=FALSE)
    list(p.value=mod$pv['Robust',1],CI=c(mod$ci['Robust',],est=mod$coef['Robust','Coeff']),bal.pval=bal.pval,
         W=c(min(abs(dat$R)),BW),n=sum(mod$Nh))
}


#' @export

cctTest <- function(dat,BW=NULL,varb,rhs=NULL,justP=TRUE){
    if(!missing(varb)) dat$Y <- -dat[[varb]]
    if(missing(BW) | is.null(BW))
        mod <- try(with(dat,rdrobust(Y,R,-0.005,kernel='uniform')))
    else  mod <- try(with(dat,rdrobust(Y,R,-0.005,kernel='uniform',h=BW)))
    if(class(mod)=='try-error') return(rep(NA,ifelse(justP,1,5)))
    if(justP) return(mod$pv['Robust',1])
    mod
}

cctMultBal <- function(dat,BW,xvars,int=FALSE){
    balPs <- lapply(xvars,function(xx) rdrobust(xx,dat$R,-0.005,h=BW))
    names(xeq) <- gsub('_','',xvars)

    sur <- systemfit(xeq,data=subset(dat,abs(R)<BW))

    rest <- paste(names(xeq),'Z',sep='_')
    rest <- paste(rest,'=0')

    linearHypothesis(sur,rest,test='Chisq')$Pr[2]
}

################
### CFT Method
##############
##' RDD analysis method modeled on Cattaneo et al 2014, J Causal Inference
##'
##'
##' @title CFT analysis method
##' @param dat Data set with columns \code{Z}, \code{R}, and an outcome variable
##' @param BW (Optional) A bandwidth \code{b>0}. If not provided it will be estimated from the data
##' @param outcome A string specifying the name of the outcome or covariate variable
##' @param alpha 1-confidence interval level
##' @param rhs Ignored
##' @return list of
##' \describe{
##'  \item{p.value} of Fisher-style no effect hypothesis
##'  \item{CI} vector of confidence limits, estimate
##'  \item{bal.pval} balance p-value
##'  \item{n} number of observations inselected window
##'  \item{W}  window
##' }
##' @author lrd author 1
##' @export
cft <- function(dat,BW,outcome,alpha=0.05,rhs=NULL){
    if(missing(BW) |is.null(BW))
        BW <- bwMult(dat,balMult.control=list(method='cft',reduced.covars=FALSE))

    if(sum(dat$R<0 & -BW<dat$R)<2 || sum(dat$R>0 & dat$R<BW)<2) return(NA)

    bal.pval=balMult(dat=dat,BW=BW,method='cft',reduced.covars=FALSE)

    estAll <- cftTest(dat=dat,BW=BW,varb=outcome,alpha=alpha,justP=FALSE)

    list(estAll$p.value,CI=c(estAll$conf.int,est=estAll$estimate),bal.pval=bal.pval,
         n=sum(abs(dat$R)<BW),W=c(min(abs(dat$R)),BW))
}

cftTest <- function(dat,BW,varb,alpha=0.05,rhs=NULL,justP=TRUE){
    if(sum(dat$R<0 & -BW<dat$R)<2 || sum(dat$R>0 & dat$R<BW)<2) return(NA)
    dat <- dat[abs(dat$R)<BW,]

    if(!missing(varb)) dat$Y <- dat[[varb]]
    if(length(unique(dat$Y))==2){
        wilcox <- with(dat,fisher.test(Y,Z,conf.int=!justP,conf.level=1-alpha))
    }
    else wilcox <- with(dat,wilcox.test(Y[Z==1],Y[Z==0],conf.int=!justP,conf.level=1-alpha))
    if(justP) return(wilcox$p.value)
    return(wilcox[c('p.value','conf.int','estimate')])
}


##############################
### data-driven bandwidths
############################




balMult <- function(dat,BW,method='sh',reduced.covars=TRUE,rhs=NULL,bonferonni=TRUE){
    ps <- NULL
    balTest <- switch(method, "sh"=balOneSH,"cft"=cftTest,"ik"=ikTest)

    xvars <- c('lhsgrade_pct','totcredits_year1','male','loc_campus1','loc_campus2',
               'bpl_north_america','english')
    xvars <- if (reduced.covars) c(xvars, 'age') else c(xvars, 'age_at_entry')

    if(method=='ik') return(ikMultBal(dat,BW,xvars))

    for(varb in xvars){
        ps <- c(ps,balTest(dat=dat,BW=BW,varb=varb,rhs=rhs))

    }
    if(bonferonni) ps <- ps*length(ps)
    min(1,min(ps))
}



#' Choose a bandwidth using either Limitless or Permutations method with one covariate, or the IK method
#'
#' @param dat a data.frame. must contain variables `R` (the running variable) and `Z`
#' (treatment) and, if method='ik', `Y` (outcome)
#' @param covname a character string, the name of the covariate to test balance on
#' (unless method='ik')
#' @param method a character string. either 'sh' for Limitless, 'ik' for IK, or 'cft' for
#' Permutations
#'
#' @return bandwidth choice (scalar)
#' @export
bw <- function(dat,covname='x',method='sh',alphax=0.15,plt=FALSE){

    if(method=='ik') return(IKbandwidth(dat$R,dat$Y,kernel='rectangular'))
    else test <- switch(method,"sh"= balOneSH, "cft"=cftTest)
    short <- function(w) test(dat=dat,BW=w,varb=covname)# - alphax
    #uniroot(short,interval=c(min(max(dat$R),-min(dat$R))-0.1,min(max(dat$R),-min(dat$R))),
    #                  extendInt='downX')$root
    bws <- seq(0.02,min(max(dat$R),-min(dat$R)),0.01)
    ps <- vapply(bws,short,1)
    if(plt) plot(ps)
    bws[max(which(ps>alphax))]
}



#' Choose a bandwidth using one of several methods
#'
#' Method options include Limitless, a permutation-testing method, or the IK method
#'
#' @param dat a data.frame. must contain variables `R` (the running variable) and `Z`
#' (treatment) and, if method='ik', `Y` (outcome) along with these covariates: 'lhsgrade_pct',
#' 'totcredits_year1','male','loc_campus1','loc_campus2','bpl_north_america','english',
#' and 'age'
#' @param alpha the level of the balance test
#' @param newbal.control list of optional arguments to `newBal`
#' @param bsw numeric, a sequence of bandwidths to try (in increasing order)
#'
#' @return bandwidth choice (scalar)
#' @export
bwMult <- function(dat,alpha=0.15,balMult.control=list(method='sh',reduced.covars=FALSE), bws = seq(0,2,0.01)){
    stopifnot(is.list(balMult.control),
              !any(c('dat', 'BW') %in% names(balMult.control)),
              is.numeric(bws) && all(bws >=0),
              length(bws)==1 || !is.unsorted(bws))
    short <- function(w) {
        nbargs <- c(list(dat=dat, BW=w), balMult.control)
        try(do.call(balMult, nbargs))
    }

    p <- 0
    i <- length(bws)+1
    while(p<alpha){
        i <- i-1
        if(i==0) return(i)
        bw <- bws[i]
        pval <- short(bw)
        p <- ifelse(is.numeric(pval)&is.finite(pval),pval,0)
    }
    bws[i]
}


##############################
### Frandsen Manipulation Test
############################
frandsenK <- function(R,BW){
    if(!missing(BW)) R <- R[abs(R)<BW]
    rvals <- sort(unique(R))
    delta <- min(c(rvals,NA)-c(NA,rvals),na.rm=TRUE)
    delta <- delta/sd(R,na.rm=TRUE)

    delta^3*dnorm(delta/2)/(2*(pnorm(3*delta/2)-pnorm(delta/2)))
}

frandsenTest <- function(N0,Nplus,Nminus,k){
    m <- N0+Nplus+Nminus
    results <- data.frame(k=k,p=NA)
    for(k in results$k){
        results$p[results$k==k] <- 2*min(pbinom(N0,m,(1-k)/(3-k)),1-pbinom(N0,m,(1+k)/(3+k)))
    }
    results
}



#' Frandsen Manipulation Test
#'
#'
#'
#' @param R a discrete running variable
#' @param cutoff the value of the running variable at the boundary between treatment and control---in Frandsen (2016), the maximum (minimum) R at which subjects are assigned to treatment; alternatively, the minimum (maximum) R at which subjects are assigned to control.
#' @param BW an optional bandwidth, helpful for choosing a benchmark for k
#' @param k a vector of possible values of k, a tuning parameter for the test. The function automatically computes the benchmark value Frandsen (2016) suggests based on normal theory, and adds it to this list.
#'
#' @return a p-value
#'
#' @export
frandsen <- function(R,cutoff,BW,k=c(0,0.01,0.02,0.1)){
    if(missing(BW)) BW <- max(abs(R),na.rm=TRUE)
    k.bench <- frandsenK(R,BW)
    k <- sort(c(k,k.bench))

    N0 <- sum(R==cutoff,na.rm=TRUE)
    rplus <- min(R[R>cutoff],na.rm=TRUE)
    rminus <- max(R[R<cutoff],na.rm=TRUE)
    Nplus <- sum(R==rplus,na.rm=TRUE)
    Nminus <- sum(R==rminus,na.rm=TRUE)

    frandsenTest(N0,Nplus,Nminus,k)
}


##' Extract bread matrix from an lmrob fit
##'
##' This is part of a workaround for an issue in the robustbase code
##' affecting sandwich covariance estimation.
##' The issue in question is issue #6471, robustbase project on R-Forge.
##'
##' @title Bread method for lmrob objects
##' @param x an lmrob object produced using an MM/SM estimator chain
##' @param ...
##' @return k by (k+1) matrix, with first column for scale estimate and rows, remaining cols for coefficients
##'
##' @author lrd author 2
##' @export
bread.lmrob <- function(x, ...)
{
  stopifnot(is.list(ctrl <- x$control))
  if (!(!is.null(ctrl$method) && nchar(ctrl$method)<=2 &&
      substr(ctrl$method, nchar(ctrl$method),nchar(ctrl$method))=="M") )
      stop("bread.lmrob() supports only SM or MM estimates")

       psi <- chi <- ctrl$psi
    if (is.null(psi))
        stop("parameter psi is not defined")
    stopifnot(is.numeric(c.chi <- ctrl$tuning.chi), is.numeric(c.psi <- ctrl$tuning.psi))
    r0 <- x$init$resid
    r <- resid(x)
    scale <- x$scale
    xmat <- model.matrix(x)
       bb <- 1/2
    n <- length(r)
    stopifnot(n == length(r0), is.matrix(xmat), n == nrow(xmat))
    p <- ncol(xmat)
    r.s <- r/scale
    r0.s <- r0/scale
    w <- robustbase::Mpsi(r.s, cc = c.psi, psi = psi, deriv = 1)
    w0 <- robustbase::Mchi(r0.s, cc = c.chi, psi = chi, deriv = 1)
    x.wx <- crossprod(xmat, xmat * w)
    if (inherits(A <- tryCatch(solve(x.wx) * scale, error = function(e) e),
        "error")) {
        A <- tryCatch(solve(x.wx, tol = 0) * scale, error = function(e) e)
        if (inherits(A, "error"))
            { stop("X'WX is singular.") } else warning("X'WX is almost singular.")
    }
    ## At this point A has no sample size scaling, as in robustbase:::.vcov.avar1
    ## The lack of scaling there precisely compensates for the lack of scaling of the crossproduct
    a <- A %*% (crossprod(xmat, w * r.s)/mean(w0 * r0.s))
    colnames(a) <- "sigma"
    ## Now we restore sample size scaling to A
    A <- n * A

    cbind(a, A)
}
##'
##'
##' Only SM or MM estimates supported
##'
##' @title Estfun method for lmrob objects
##' @param x a fitted lmrob
##' @param ...
##' @return an estfun object, as in the sandwich package
##' @author lrd author 2
##' @export
estfun.lmrob <- function(x, ...)
{
  stopifnot(is.list(ctrl <- x$control))
  if (!(!is.null(ctrl$method) && nchar(ctrl$method)<=2 &&
      substr(ctrl$method, nchar(ctrl$method),nchar(ctrl$method))=="M") )
      stop("estfun.lmrob() supports only SM or MM estimates")

  xmat <- model.matrix(x)
    xmat <- naresid(x$na.action, xmat)
       psi <- chi <- ctrl$psi
    if (is.null(psi))
        stop("parameter psi is not defined")
    stopifnot(is.numeric(c.chi <- ctrl$tuning.chi), is.numeric(c.psi <- ctrl$tuning.psi))
    r0 <- x$init$resid
    r <- resid(x)
    scale <- x$scale

    n <- length(r)
    stopifnot(n == length(r0), is.matrix(xmat), n == nrow(xmat))
    p <- ncol(xmat)
    r0.s <- r0/scale
    w0 <- robustbase::Mchi(r0.s, cc = c.chi, psi = chi)
    Usigma <- scale(w0, center=TRUE, scale=FALSE)

    r.s <- r/scale
    w <- robustbase::Mpsi(r.s, cc = c.psi, psi = psi)

  Ubeta <- w * xmat
    rval <- cbind("sigma"=Usigma, Ubeta)
       attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    rval
    }

##' Overloading of sandwich::sandwich to accommodate non-square bread
##'
##' The sandwich package's sandwich function presumes the bread matrix
##' to be symmetric. Obviously this won't do if the bread is rectangular but not
##' square.
##'
##' This is part of a workaround for an issue in the robustbase code
##' affecting sandwich covariance estimation.
##' The issue in question is issue #6471, robustbase project on R-Forge.
##'
##' @title Sandwich estimate of covariance
##' @param x a fitted model object, as in sandwich::sandwich
##' @param bread. function or matrix,
##' @param meat. function or matrix, as in sandwich::sandwich
##' @param ... additional arguments to downstream methods, as in sandwich::sandwich
##' @return matrix, bread %*% meat %*% t(bread)
##' @author lrd author 2
LRDsandwich <- function (x, bread. = sandwich::bread, meat. = sandwich::meat, ...)
{
    if (is.list(x) && !is.null(x$na.action))
        class(x$na.action) <- "omit"
    if (is.function(bread.))
        bread. <- bread.(x)
    if (is.function(meat.))
        meat. <- meat.(x, ...)
    n <- NROW(sandwich::estfun(x))
    ## the t() in the below is the only difference from sandwich::sandwich()
    return(1/n * (bread. %*% meat. %*% t(bread.)))
}


.vcov.avar2 <- function(obj, x=obj$x, posdef.meth = c("posdefify","orig"))
{ ## was .vcov.MM
    stopifnot(is.list(ctrl <- obj$control))
    ## works only for MM & SM estimates:
###    if (!is.null(ctrl$method) && !ctrl$method %in% c('SM', 'MM'))
### I replaced line above w/ 2 lines below for own reasons unrelated to bug fix -BH
    if (!is.null(ctrl$method) && !(nchar(ctrl$method)==2 &&
      substr(ctrl$method, nchar(ctrl$method),nchar(ctrl$method))=="M") )
    stop('.vcov.avar2() supports only SM or MM estimates')
    ## set psi and chi constants
    psi <- chi <- ctrl$psi
    if (is.null(psi)) stop('parameter psi is not defined')
    stopifnot(is.numeric(c.chi <- ctrl$tuning.chi),
	      is.numeric(c.psi <- ctrl$tuning.psi))

    ## need (r0, r, scale, x, c.psi,c.chi, bb)
    r0 <- obj$init$resid
    r <- resid(obj)
    scale <- obj$scale
    if (is.null(x))  x <- model.matrix(obj)
    bb <- 1/2 ## this is always 1/2 for S estimates by convention
### --- start code from .vcov.MM ---
    ## scaled residuals
    n <- length(r)
    stopifnot(n == length(r0), is.matrix(x), n == nrow(x))
    p <- ncol(x)
    ## Next 2 lines added post-.vcov.MM, addressing #6471.
    ## This assumes initial S-estimate solved sum( loss )/(n-p) == bb,
    ## not mean( loss ) == bb as assumed in .vcov.avar1
    adj <- (n-p)/n
    bb <- bb * adj
    r.s	 <- r / scale # final   scaled residuals
    r0.s <- r0 / scale # initial scaled residuals
    w  <- Mpsi(r.s, cc = c.psi, psi = psi, deriv = 1)
    w0 <- Mchi(r0.s, cc = c.chi, psi = chi, deriv = 1)
    ## FIXME for multivariate y :
    x.wx <- crossprod(x, x * w)
    if(inherits(A <- tryCatch(solve(x.wx) * scale,
			      error=function(e)e), "error")) {
	warning("X'WX is almost singular. Consider rather using cov = \".vcov.w\"")
	A <- tryCatch(solve(x.wx, tol = 0) * scale, error=function(e)e)
	if(inherits(A, "error"))
	    stop("X'WX is singular. Rather use cov = \".vcov.w\"")
    }
    a <- A %*% (crossprod(x, w * r.s) / mean(w0 * r0.s))
    w <- Mpsi( r.s, cc = c.psi, psi = psi)

    ## 3) now the standard part  (w, x, r0.s,  n, A,a, c.chi, bb)
    w0 <- Mchi(r0.s, cc = c.chi, psi = chi)
    Xww <- crossprod(x, w*w0)
    u1 <- A %*% crossprod(x, x * w^2) %*% (n * A)
    u2 <- a %*% crossprod(Xww, A)
    u3 <- A %*% tcrossprod(Xww, a)
    u4 <- mean(w0^2 - bb^2) * tcrossprod(a)

    ## list(cov = matrix((u1 - u2 - u3 + u4)/n, p, p),
    ##      wt = w / r.s, a = a)
### --- end code from .vcov.MM ---
    ret <- (u1 - u2 - u3 + u4)/n

    ## this might not be a positive definite matrix
    ## check eigenvalues (symmetric: ensure non-complex)
    ev <- eigen(ret, symmetric = TRUE)
    if (any(neg.ev <- ev$values < 0)) { ## there's a problem
	posdef.meth <- match.arg(posdef.meth)
	if(ctrl$trace.lev)
	    message("fixing ", sum(neg.ev),
		    " negative eigen([",p,"])values")
	Q <- ev$vectors
	switch(posdef.meth,
	       "orig" = {
		   ## remove negative eigenvalue:
		   ## transform covariance matrix into eigenbasis
		   levinv <- solve(Q)
		   cov.eb <- levinv %*% ret %*% Q
		   ## set vectors corresponding to negative ev to zero
		   cov.eb[, neg.ev] <- 0
		   ## cov.eb[cov.eb < 1e-16] <- 0
		   ## and transform back
		   ret <- Q %*% cov.eb %*% levinv
	       },
	       "posdefify" = {
		   ## Instead of using	require("sfsmisc") and
		   ## ret <- posdefify(ret, "someEVadd",eigen.m = ev,eps.ev = 0)
		   lam <- ev$values
		   lam[neg.ev] <- 0
		   o.diag <- diag(ret)# original one - for rescaling
		   ret <- Q %*% (lam * t(Q)) ## == Q %*% diag(lam) %*% t(Q)
		   ## rescale to the original diagonal values
		   ## D <- sqrt(o.diag/diag(ret))
		   ## where they are >= 0 :
		   D <- sqrt(pmax.int(0, o.diag)/diag(ret))
		   ret[] <- D * ret * rep(D, each = p) ## == diag(D) %*% m %*% diag(D)
	       },
	       stop("invalid 'posdef.meth': ", posdef.meth))
    }
    attr(ret,"weights") <- w / r.s
    attr(ret,"eigen") <- ev
    ret
}## end{.vcov.avar2}

