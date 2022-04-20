#' @title The Cressie-Read test for given probabilities
#'
#' @description
#' Computes D-squared Cressie-Read test for given probabilities.
#' @param x a numeric vector of values.
#' @param lambda lambda parameter lambda test Cressie-Read.
#' @param p probabilities.
#' @details
#' Statistics Cressie-Read:
#'
#'  D = 2/(lambda*(lambda+1)) * sum( O * ((O/np)^lambda - 1) )
#'
#'  > O - the number of observed
#'
#'  > np - n*p
#'
#'  Lambda parameter to be different from 0, and -1.
#'
#'  > lambda=  1 (Pearson test)
#'
#'  > lambda= -2 (Neyman's modified Pearson test)
#'
#'  > lambda~  0, lambda is close to 0 (Likelihood Ratio test)
#'
#'  > lambda~ -1 lambda is close to -1 (Kullback-Leibler modified Likelihood Ratio test)
#'
#'  > lambda= -1/2 (Freeman-Tukey's test)
#' @usage cr.gof(x, lambda=2/3, p=rep(1/length(x), length(x)))
#' @return An object of class "htest" containing the following components:
#' \describe{
#' \item{parameter}{the degrees of freedom and lambda parameter}
#' \item{statistic}{statistic test value}
#' \item{p.value}{p-value}
#'          }
#' @keywords cr.gof
#' @author
#' Krzysztof Trajkowski
#' @examples
#' # data:
#' s <- c(12,15,20,25,17,28)
#' p <- 1/6
#' n <- sum(s)*p
#' np <- n*p
#' # Cressie-Read:
#' cr.gof(s)
#'
#' # Pearson:
#' cr.gof(s,lambda=1)
#' @rdname cr.gof
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dnbinom dnorm nlminb pchisq pgamma ppois pnbinom pnorm sd var
#'
#' @seealso \code{\link{chisq.test}}
#' @export
#'

cr.gof=function(x, lambda=2/3, p=rep(1/length(x), length(x)) )
{
  DNAME <- deparse(substitute(x))
  Df <- length(x)-1
  E= sum(x)*p
  D <- (2/(lambda*(lambda+1)))*sum(x*( (x/E)^(lambda)-1 ))
  pvalue <- pchisq(D,Df,lower.tail=FALSE)

  if (lambda==2/3 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda==1 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test (Pearson) for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda <= 1e-05 & lambda > 0  & sum(p)==1 | lambda < 0 && lambda >= -1e-05 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test (G-squared) for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda <= -0.99999 & lambda > -1  & sum(p)==1 | lambda < -1 & lambda >= -1.00001 & sum(p)==1 ) {
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test (Kullback-Leibler) for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda == -1/2 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test (Freeman-Tukey's) for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda == -2 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test (Neyman's) for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(sum(p)==1){
    RVAL <- list(statistic = c(D = D), parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                 method = "D-squared Cressie-Read test for given probabilities",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  stop("error")
}
