#' @title The Cressie-Read test for discrete data.
#'
#' @description
#' Goodnes of fit test for two distributions: Poisson, binomial and negative binomial.
#' @param x a numeric vector of values.
#' @param dist names distributions: pois (deflaut), binom or nbinom
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
#'
#' @rdname chi2.test
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dbinom dnbinom dnorm nlminb pchisq pgamma ppois pbinom pnbinom pnorm sd var
#'
#' @seealso \code{\link[vcd]{goodfit}}
#' @export

chi2.test=function(x,dist="pois", lambda = 1)
{
  DNAME <- deparse(substitute(x))
  lambda <- lambda
  mu <- mean(x)
  DF <- as.data.frame(table(x))
  DF$x <- as.numeric(as.character(DF$x))
  tab <- merge(DF,data.frame(x=0:max(DF$x)),by.x = "x",all = TRUE)
  tab$Freq <- ifelse(is.na(tab$Freq), 0, tab$Freq)
  f <- tab$Freq
  v <- tab$x
  # Poisson:
  dP <- dpois(v, mu)
  pP <- replace(dP,length(dP),ifelse(sum(dP)==1, dP[length(dP)],dP[length(dP)]+1-sum(dP)))
  eP <- pP*sum(f)
  # Negative Binomial:
  size <- mean(x)^2/(var(x)-mean(x))
  dNB <- sapply(1:length(v),function(i) ifelse(size<0,NA,dnbinom(v, mu=mu,size=size)[i]))
  # dNB <- dnbinom(v, mu=mu,size=size)
  pNB <- replace(dNB,length(dNB),ifelse(sum(dNB)==1, dNB[length(dNB)],dNB[length(dNB)]+1-sum(dNB)))
  eNB <- pNB*sum(f)
  # Binomial:
  vr <- var(x)
  sizeB <- floor(-mu^2/(vr-mu))
  probB <- -(vr-mu)/mu
  dB <- sapply(1:length(v),function(i) ifelse(sizeB<0 | probB<0,NA,dbinom(v, size=sizeB,prob=probB)[i]))
  # dB <- dbinom(v, size=sizeB,prob=probB)
  pB <- replace(dB,length(dB),ifelse(sum(dB)==1, dB[length(dB)],dB[length(dB)]+1-sum(dB)))
  eB <- pB*sum(f)
  if (dist=="pois") {
    c <- (2/(lambda*(lambda+1)))*sum(f*( (f/eP)^(lambda)-1 ))
    Df <- length(f)-1-1
    pvalue <- 1-pchisq(c, df=Df)
    L <- sum(dpois(x,mu,log=T))
    RVAL <- list(statistic = c(CR=c), parameter = c(df=Df,lambda=lambda), p.value = pvalue, estimate =c(mu=mean(x),logLik=L),
                 method = "Cressie-Read Goodness-of-Fit Test - Poisson distributions",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(dist=="nbinom") {
    if (size>0) {
      c <- (2/(lambda*(lambda+1)))*sum(f*( (f/eNB)^(lambda)-1 ))
      Df <- length(f)-1-1
      pvalue2 <- 1-pchisq(c, df=Df)
      GG <- sum(dnbinom(x,mu=mu,size=size,log=T))
      RVAL <- list(statistic = c(CR=c), parameter = c(df=Df,lambda=lambda), p.value = pvalue2, estimate =c(mu=mu, size=size, loglik=GG),
                   method = "Cressie-Read Goodness-of-Fit Test - nbinomial distributions",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else list(warning=c(mu=mu,size=size))
  }
  else if(dist=="binom") {
    if (sizeB>0 | probB>0) {
      c <- (2/(lambda*(lambda+1)))*sum(f*( (f/eB)^(lambda)-1 ))
      Df <- length(f)-1-1
      pvalue2 <- 1-pchisq(c, df=Df)
      GG <- sum(dbinom(x,size=sizeB,prob=probB,log=T))
      RVAL <- list(statistic = c(CR=c), parameter = c(df=Df,lambda=lambda), p.value = pvalue2, estimate =c(size=sizeB, prob=probB, loglik=GG),
                   method = "Cressie-Read Goodness-of-Fit Test - binomial distributions",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else list(warning=c(size=sizeB,prob=probB))
  }
}
