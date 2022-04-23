#' @title The Cressie-Read test for continuous data.
#'
#' @description
#' Goodnes of fit test for distributions: normal, gamma, beta and betap (aka Beta Prime).
#' @param x a numeric vector of values.
#' @param dist names distributions: norm (deflaut), gamma, beta and betap (aka Beta Prime).
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
#' @importFrom stats dgamma dpois dnbinom dnorm dbeta nlminb pchisq pgamma ppois pnbinom pnorm pbeta sd var
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
    LL <- with(function(par,data) {
      shape1= par[1]
      shape2= par[2]
      -sum(fB(x=x,shape1=shape1,shape2 = shape2))
    },data = data.frame(x=x))
    MU = mean(x)
    VR = var(x)
    res <- nlminb(start = c(shape1=(MU*VR+MU^3+MU^2)/VR ,shape2=(2*VR+MU^2+MU)/VR ),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    num4 <- floor(1 + n.classes * pbeta(x, shape1=G1, shape2=G2))
    count4 <- tabulate(num4, n.classes)
    ST4 <- (2/(lambda*(lambda+1)))*sum(count4*( (count4/(n*prob))^(lambda)-1 ))
    pvalue4 <- pchisq(ST4, n.classes - 2 - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(CR=ST4), parameter = c(df=Df,lambda=lambda), p.value = pvalue4, estimate =c(shape1=G1, shape2=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - aka Beta Prime distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
}





