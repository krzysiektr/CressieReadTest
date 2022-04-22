#' @title The Cressie-Read test for contingency table.
#'
#' @description
#' Calculates statistics X-squared, G-squared, F-squared for contingency table.
#'
#' @param x contingency table.
#' @param test "gw","g", "p", "n", "ft", "kl"
#' @usage cgf.test(x, test="pc")
#' @return An object of class "htest" containing the following components:
#' \describe{
#' \item{parameter}{the degrees of freedom}
#' \item{statistic}{statistic test value}
#' \item{p.value}{p-value}
#'          }
#' @keywords cgf.test
#' @author
#' Krzysztof Trajkowski
#' @examples
#' m <- matrix(c(23,32,45,26),2,2)
#'
#' # LR test with correction:
#' cgf.test(m, test="gw")
#'
#' # LR test without correction:
#' cgf.test(m, test="g")
#'
#' # test Pearson:
#' cgf.test(m, test="p")
#'
#' # test Pearson with correction:
#' cgf.test(m)
#'
#' @rdname cgf.test
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dnbinom dnorm nlminb pchisq pgamma ppois pnbinom pnorm sd var
#'
#' @seealso \code{\link[vcd]{assocstats}}
#' @export
#'

cgf.test=function(x, test="pc")
{
  DNAME <- deparse(substitute(x))
  mr <- apply(x,1,sum); mk <- apply(x,2,sum); sm <- sum(mr)
  R <- (sum(1/mr))*sm-1; K <- (sum(1/mk))*sm-1
  W <- 1+(R*K)/(6*sm*(length(mr)-1)*(length(mk)-1))
  E <- outer(mr,mk,"*")/sm
  G <- 2*(sum(x*log(x/E)))
  P <- sum((x-E)^2/E)
  C1 <- sqrt(P/sm) # Y-Yule'a
  C2 <- sqrt(P/(sm+P)) # C-Pearson
  C3 <- sqrt((P)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) ) # V-Cramer
  C4 <- sqrt((P)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) ) # T-Tschuprow
  Pc <- sum((abs(x - E) - 0.5)^2/E)
  FT <- 4*sum((sqrt(x)-sqrt(E))^2)
  N <- sum((x-E)^2/x)
  KL <- 2*sum(E*log(E/x))
  Df <- (nrow(x)-1)*(ncol(x)-1)
  Gw <- G/W
  pvalueA <- pchisq(G,Df,lower.tail=FALSE)
  pvalueB <- pchisq(Gw,Df,lower.tail=FALSE)
  pvalueP <- pchisq(P,Df,lower.tail=FALSE)
  pvalueFT <- pchisq(FT,Df,lower.tail=FALSE)
  pvalueKL <- pchisq(KL,Df,lower.tail=FALSE)
  pvalueN <- pchisq(N,Df,lower.tail=FALSE)
  pvaluePc <- pchisq(Pc,Df,lower.tail=FALSE)
  if (test=="gw") {
    RVAL <- list(statistic = c(Gw = Gw), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4), parameter = c(df = Df), p.value = pvalueB,
                 method = "G-squared Likelihood Ratio test with Williams correction",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="g") {
    RVAL <- list(statistic = c(G = G), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4), parameter = c(df = Df), p.value = pvalueA,
                 method = "G-squared Likelihood Ratio test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="p") {
    RVAL <- list(statistic = c(P = P), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4),  parameter = c(df = Df), p.value = pvalueP,
                 method = "X-squared Pearson test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
    else if(test=="pc") {
      RVAL <- list(statistic = c(Pc = Pc), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4),  parameter = c(df = Df), p.value = pvaluePc,
                   method = "X-squared Pearson test with Yates' continuity correction",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
  }
  else if(test=="ft") {
    RVAL <- list(statistic = c(FT = FT), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4),  parameter = c(df = Df), p.value = pvalueFT,
                 method = "F-squared Freeman-Tukey's test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="kl") {
    RVAL <- list(statistic = c(KL = KL), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4),  parameter = c(df = Df), p.value = pvalueKL,
                 method = "G-squared Kullback-Leibler test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="n") {
    RVAL <- list(statistic = c(N = N), estimate=c(Y_Yulea = C1,C_Pearson = C2,V_Cramer = C3,T_Tschuprow = C4), parameter = c(df = Df), p.value = pvalueN,
                 method = "X-squared Neyman's test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
}
