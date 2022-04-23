#' @title The Cressie-Read test for continuous data.
#'
#' @description
#' Goodnes of fit test for distributions: normal, t-Student, cauchy, laplace, gamma, beta and betap (aka Beta Prime).
#' @param x a numeric vector of values.
#' @param dist names distributions: norm (deflaut), t-Student, cauchy, laplace, gamma, beta and betap (aka Beta Prime).
#' @param lambda numeric parametr
#' @usage chi2.gof(x, dist="norm", lambda=1)
#' @return An object of class "htest" containing the following components:
#' \describe{
#' \item{estimate}{estimates distribution parameters and loglik}
#' \item{parameter}{the degrees of freedom and lambda parameter}
#' \item{statistic}{statistic test value}
#' \item{p.value}{p-value}
#'          }
#' @note estimate parameter in \code{\link[fitdistrplus]{fitdist}}.
#' @keywords chi2.gof
#' @author
#' Krzysztof Trajkowski
#' @examples
#' set.seed(2305)
#' g <- rnorm(120)
#' chi2.gof(g)
#' @rdname chi2.gof
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dcauchy pcauchy dgamma dpois dnbinom dnorm dbeta nlminb pchisq pgamma ppois pnbinom pnorm pbeta sd var IQR median integrate dt pt
#'
#' @seealso \code{\link[nortest]{pearson.test}} \code{\link{ks.test}} \code{\link[ADGofTest]{ad.test}} \code{\link[vsgoftest]{vs.test}}
#' @export

chi2.gof=function(x, dist="norm", lambda=1)
{
  DNAME <- deparse(substitute(x))
  lambda <- lambda
  n <- length(x)
  n.classes <- ceiling(2*n^(2/5))
  prob <- rep(1/n.classes, n.classes)
  num1 <- floor(1 + n.classes * pnorm(x, mean(x), sd(x)))
  count1 <- tabulate(num1, n.classes)
  ST1 <- (2/(lambda*(lambda+1)))*sum(count1*( (count1/(n*prob))^(lambda)-1 ))
  pvalue1 <- pchisq(ST1, n.classes - 2 - 1, lower.tail = FALSE)
  Df <- length(prob)-2-1

  if (dist=="norm") {
    RVAL <- list(statistic = c(CR=ST1), parameter = c(df=Df,lambda=lambda), p.value = pvalue1, estimate =c(mu=mean(x),std=sd(x),loglik=sum(dnorm(x,mean(x), sd(x),log=T))),
                 method = "Cressie-Read Goodness-of-Fit Test - normal distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="gamma") {
    num2 <- floor(1 + n.classes * pgamma(x, shape=mean(x)^2/var(x), scale=var(x)/mean(x)))
    count2 <- tabulate(num2, n.classes)
    ST2 <- (2/(lambda*(lambda+1)))*sum(count2*( (count2/(n*prob))^(lambda)-1 ))
    pvalue2 <- pchisq(ST2, n.classes - 2 - 1, lower.tail = FALSE)
    LL <- with(function(par,data) {
      shape= par[1]
      scale= par[2]
      -sum(dgamma(x=x,shape=shape,scale=scale,log=TRUE))
    },data = data.frame(x=x))
    res <- nlminb(start = c(shape=mean(x)^2/var(x),scale=var(x)/mean(x)),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    RVAL <- list(statistic = c(CR=ST2), parameter = c(df=Df,lambda=lambda), p.value = pvalue2, estimate =c(shape=G1, scale=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - gamma distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="beta") {
    LL <- with(function(par,data) {
      shape1= par[1]
      shape2= par[2]
      -sum(dbeta(x=x,shape1=shape1,shape2 = shape2,log=TRUE))
      },data = data.frame(x=x))
    res <- nlminb(start = c(shape1=-(mean(x)*var(x)+mean(x)^3-mean(x)^2)/var(x),shape2=((mean(x)-1)*var(x)+mean(x)^3-2*mean(x)^2+mean(x))/var(x)),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    num3 <- floor(1 + n.classes * pbeta(x, shape1=G1, shape2=G2))
    count3 <- tabulate(num3, n.classes)
    ST3 <- (2/(lambda*(lambda+1)))*sum(count3*( (count3/(n*prob))^(lambda)-1 ))
    pvalue3 <- pchisq(ST3, n.classes - 2 - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(CR=ST3), parameter = c(df=Df,lambda=lambda), p.value = pvalue3, estimate =c(shape1=G1, shape2=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - beta distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="betap") {
    fB <- function(x,shape1,shape2) log(x^(shape1-1))-(shape1+shape2)*log(1+x)-lbeta(shape1,shape2)
    pB <- function(shape1,shape2,h) integrate(function(x) (x^(shape1-1)*(1+x)^(-shape1-shape2))/beta(shape1,shape2),0,h)
    LL <- with(function(par,data) {
      shape1= par[1]
      shape2= par[2]
      -sum(fB(x = x, shape1 = shape1, shape2 = shape2))
    },data = data.frame(x=x))
    MU = mean(x)
    VR = var(x)
    res <- nlminb(start = c(shape1=(MU*VR+MU^3+MU^2)/VR ,shape2=(2*VR+MU^2+MU)/VR ),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    num4 <- floor(1 + n.classes * sapply(1:length(x),function(i) pB(shape1=G1, shape2=G2,x[i])$value))
    count4 <- tabulate(num4, n.classes)
    ST4 <- (2/(lambda*(lambda+1)))*sum(count4*( (count4/(n*prob))^(lambda)-1 ))
    pvalue4 <- pchisq(ST4, n.classes - 2 - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(CR=ST4), parameter = c(df=Df,lambda=lambda), p.value = pvalue4, estimate =c(shape1=G1, shape2=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - aka Beta Prime distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="t") {
    dt_ls <- function(x, df=df, mu=mu, sigma=sigma) 1/sigma * dt((x - mu)/sigma, df)
    pt_ls <- function(q, df=df, mu=mu, sigma=sigma)  pt((q - mu)/sigma, df)
    LL <- with(function(par,data) {
      df= par[1]
      mu= par[2]
      sigma = par[3]
      -sum(log(dt_ls(x = x, df=df, mu=mu, sigma=sigma)))
    },data = data.frame(x=x))
    res <- nlminb(start = c(df=3 ,mu=mean(x),sigma=sd(x) ),lower=c(1,-1, 0.00001),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    G3 <- as.numeric(res$par[3])
    GG <- -1*as.numeric(res$objective)
    num5 <- floor(1 + n.classes * pt_ls(x, df=G1, mu=G2, sigma=G3))
    count5 <- tabulate(num5, n.classes)
    ST5 <- (2/(lambda*(lambda+1)))*sum(count5*( (count5/(n*prob))^(lambda)-1 ))
    pvalue5 <- pchisq(ST5, n.classes - 2 - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(CR=ST5), parameter = c(df=Df,lambda=lambda), p.value = pvalue5, estimate =c(df=G1, mu=G2, std=G3,loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - t-Student distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="cauchy") {
    LL <- with(function(par,data) {
      location= par[1]
      scale= par[2]
      -sum(dcauchy(x=x, location = location, scale = scale, log=TRUE))
    },data = data.frame(x=x))
    res <- nlminb(start = c(location=median(x), scale=IQR(x)/2),lower=c(-Inf, 0.00001),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    num6 <- floor(1 + n.classes * pcauchy(x, location=G1, scale=G2))
    count6 <- tabulate(num6, n.classes)
    ST6 <- (2/(lambda*(lambda+1)))*sum(count6*( (count6/(n*prob))^(lambda)-1 ))
    pvalue6 <- pchisq(ST6, n.classes - 2 - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(CR=ST6), parameter = c(df=Df,lambda=lambda), p.value = pvalue6, estimate =c(location=G1, scale=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - Cauchy distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="laplace") {
    dL <- function(x,location,scale) exp(-abs(x-location)/scale)/(2*scale)
    pL <- function(h,location,scale) integrate(function(x) exp(-abs(x-location)/scale)/(2*scale),-Inf,h)
    LL <- with(function(par,data) {
      location= par[1]
      scale= par[2]
      -sum(log(dL(x=x, location = location, scale = scale)))
    },data = data.frame(x=x))
    res <- nlminb(start = c(location=median(x), scale=IQR(x)/2),lower=c(-Inf, 0.00001),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    num7 <- floor(1 + n.classes * sapply(1:length(x),function(i) pL(x[i],location=G1, scale=G2)$value))
    count7 <- tabulate(num7, n.classes)
    ST7 <- (2/(lambda*(lambda+1)))*sum(count7*( (count7/(n*prob))^(lambda)-1 ))
    pvalue7 <- pchisq(ST7, n.classes - 2 - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(CR=ST7), parameter = c(df=Df,lambda=lambda), p.value = pvalue7, estimate =c(location=G1, scale=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - Laplace distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
}





