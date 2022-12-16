#' @title The observed and expected counts.
#'
#' @description
#' The observed and expected counts under the null hypothesis: Poisson, binomial and negative binomial.
#' @param x a numeric vector of values.
#' @param dist names distributions: pois (deflaut), binom or nbinom
#' @usage obs(x, dist = "pois")
#'
#' @return An object containing the following components:
#' \describe{
#' \item{ob}{observed counts}
#' \item{ex}{expected counts}
#'          }
#' @keywords chi2.test
#' @author
#' Krzysztof Trajkowski
#' @examples
#' set.seed(2305)
#' g <- rpois(120,3)
#' obs(g)
#' @rdname obs.test
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dpois dbinom dnbinom nlminb pchisq ppois pbinom pnbinom sd var
#'
#' @seealso \code{\link[stats]{chisq.test}}
#' @export


obs = function(x,dist="pois")
{
  DNAME <- deparse(substitute(x))
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
    return(list(dist="pois",ob=f,ex=eP))
  }
  else if(dist=="nbinom") {
    if (size>0) {
      return(list(dist="nbinom",ob=f,ex=eNB))
    }
    else list(dist="nbinom",warning=c(mu=mu,size=size))
  }
  else if(dist=="binom") {
    if (sizeB>0 | probB>0) {
      return(list(dist="binom",ob=f,ex=eB))
    }
    else list(dist="binom",warning=c(size=sizeB,prob=probB))
  }
}
