#' @title The Cressie-Read test for contingency table.
#'
#' @description
#' Computes D-squared Cressie-Read test for contingency table.
#' @param x contingency table.
#' @param lambda lambda parameter lambda test Cressie-Read.
#' @usage cr.test(x, lambda=2/3)
#' @return An object of class "htest" containing the following components:
#' \describe{
#' \item{estimate}{estimates contingency coefficient}
#' \item{parameter}{the degrees of freedom and lambda parameter}
#' \item{statistic}{statistic test value}
#' \item{p.value}{p-value}
#'          }
#' @keywords cr.test
#' @author
#' Krzysztof Trajkowski
#' @examples
#' # data:
#' m <- matrix(c(23,32,45,26),2,2)
#'
#' # Cressie-Read:
#' cr.test(m)
#'
#' # Pearson:
#' cr.test(m,lambda=1)
#'
#' # Likelihood Ratio:
#' cr.test(m,lambda=1e-05)
#'
#' # Freeman-Tukey's:
#' cr.test(m,lambda=-1/2)
#'
#' # Neyman's:
#' cr.test(m,lambda=-2)
#'
#' # Kullback-Leibler:
#' cr.test(m,lambda=-0.99999)
#' @rdname cr.test
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dnbinom dnorm nlminb pchisq pgamma ppois pnbinom pnorm sd var
#'
#' @seealso \code{\link{chisq.test}}   \code{\link{fisher.test}}
#' @export
#'

cr.test=function(x, lambda=2/3)
{
  DNAME <- deparse(substitute(x))
  mr <- apply(x,1,sum); mk <- apply(x,2,sum); sm <- sum(mr)
  E <- outer(mr,mk,"*")/sm
  Df <- (nrow(x)-1)*(ncol(x)-1)
  D <- (2/(lambda*(lambda+1)))*sum(x*( (x/E)^(lambda)-1 ))
  CHI2 <- (2/(1*(1+1)))*sum(x*( (x/E)^(1)-1 ))
  C1 <- sqrt(CHI2/sm) # phi 2x2
  C2 <- sqrt(CHI2/(sm+CHI2)) # C-Pearson
  wp <- min(dim(x))
  Cc <- C2/sqrt((wp-1)/wp) # Cadj-Pearson
  C3 <- sqrt(CHI2/(sm*(wp-1)) ) # V-Cramer
  C4 <- sqrt((CHI2)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) ) # T-Tschuprow
  pvalue <- pchisq(D,Df,lower.tail=FALSE)

  if (lambda==2/3) {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
  else if(lambda==1) {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Pearson)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Pearson)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
  else if(lambda <= 1e-05 && lambda > 0 || lambda < 0 && lambda >= -1e-05) {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (G-squared)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (G-squared)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
  else if(lambda <= -0.99999 && lambda > -1 || lambda < -1 && lambda >= -1.00001 ) {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Kullback-Leibler)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Kullback-Leibler)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
  else if(lambda == -1/2) {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Freeman-Tukey's)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Freeman-Tukey's)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
  else if(lambda == -2) {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Neyman's)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test (Neyman's)",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
  else {
    if (dim(x)[1] == 2 & dim(x)[2] == 2) {
      RVAL <- list(statistic = c(D = D), estimate=c(Q_Yulea = (x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
                                                    Phi = C1,C_Pearson = C2, C_adj = Cc,V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
    else {
      RVAL <- list(statistic = c(D = D), estimate=c(C_Pearson = C2, C_adj = Cc, V_Cramer = C3,T_Tschuprow = C4),
                   parameter = c(lambda = lambda, df = Df), p.value = pvalue,
                   method = "D-squared Cressie-Read test",
                   data.name = DNAME)
      class(RVAL) <- "htest"
      return(RVAL)
    }
  }
}
