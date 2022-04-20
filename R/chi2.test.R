#' @title The Cressie-Read test for discrete data.
#'
#' @description
#' Goodnes of fit test for two distributions: Poisson and negative binomial.
#' @param x a numeric vector of values.
#' @param dist names distributions: pois (deflaut) or nbinom
#' @param lambda parameter lambda test Cressie-Read.
#' @usage chi2.test(x, dist = "pois", lambda = 1)
#'
#' @return An object of class "htest" containing the following components:
#' \describe{
#' \item{estimate}{estimates distribution parameters and loglik}
#' \item{parameter}{the degrees of freedom and lambda parameter}
#' \item{statistic}{statistic test value}
#' \item{p.value}{p-value}
#'          }
#' @keywords chi2.test
#' @author
#' Krzysztof Trajkowski
#' @examples
#' set.seed(2305)
#' g <- rpois(120,3)
#' chi2.test(g)
#' @rdname chi2.test
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dnbinom dnorm nlminb pchisq pgamma ppois pnbinom pnorm sd var
#'
#' @seealso \code{\link[vcd]{goodfit}}
#' @export

chi2.test=function(x,dist="pois", lambda = 1)
{
  DNAME <- deparse(substitute(x))
  lambda <- lambda
  mu <- mean(x)
  size <- mean(x)^2/(var(x)-mean(x))
  f <- as.numeric(table(x))
  v <- as.numeric(names(table(x)))
  d <- dpois(v, mu)
  p <- replace(d,length(d),ifelse(sum(d)==1, d[length(d)],d[length(d)]+1-sum(d)))
  e <- p*sum(f)
  c <- (2/(lambda*(lambda+1)))*sum(f*( (f/e)^(lambda)-1 ))
  Df <- length(f)-1-1
  pvalue <- 1-pchisq(c, df=Df)
  L <- sum(dpois(x,mu,log=T))

  if (dist=="pois") {
    RVAL <- list(statistic = c(CR=c), parameter = c(df=Df,lambda=lambda), p.value = pvalue, estimate =c(mu=mean(x),logLik=L),
                 method = "Cressie-Read Goodness-of-Fit Test - Poisson distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="nbinom") {
    LL <- with(function(par,data) {
      mu= par[1]
      size= par[2]
      -sum(dnbinom(x=x,mu=mu,size=size,log=TRUE))
    },data = data.frame(x=x))
    res <- nlminb(start = c(mu=mu,size=size),lower = c(-Inf,0), upper = c(Inf,Inf),objective=LL)
    G1 <- as.numeric(res$par[1])
    G2 <- as.numeric(res$par[2])
    GG <- -1*as.numeric(res$objective)
    d <- dnbinom(v, mu=G1,size=G2)
    p <- replace(d,length(d),ifelse(sum(d)==1, d[length(d)],d[length(d)]+1-sum(d)))
    e <- p*sum(f)
    c <- (2/(lambda*(lambda+1)))*sum(f*( (f/e)^(lambda)-1 ))
    Df <- length(f)-1-1
    pvalue2 <- 1-pchisq(c, df=Df)
    RVAL <- list(statistic = c(CR=c), parameter = c(df=Df,lambda=lambda), p.value = pvalue2, estimate =c(mu=G1, size=G2, loglik=GG),
                 method = "Cressie-Read Goodness-of-Fit Test - nbinomial distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
}
