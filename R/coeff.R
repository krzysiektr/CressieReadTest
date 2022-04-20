#' @title The contingency coefficient
#'
#' @description
#' Computes contingency coefficient: C-Pearson, T-Czupurow, Fi-Yule'a, V-Cramer.
#' @param x contingency table.
#' @usage coeff(x)
#'
#' @keywords coeff
#' @author
#' Krzysztof Trajkowski
#' @examples
#' # data:
#' m <- matrix(c(23,32,45,26),2,2)
#' # contingency coefficient:
#' coeff(m)
#' @rdname coeff
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dnbinom dnorm nlminb pchisq pgamma ppois pnbinom pnorm sd var
#'
#' @seealso \code{\link{chisq.test}} \code{\link[vcd]{assocstats}}
#' @export
#'

coeff <- function(x)
{
  mr <- apply(x,1,sum); mk <- apply(x,2,sum); sm <- sum(mr)
  R <- (sum(1/mr))*sm-1; K <- (sum(1/mk))*sm-1
  W <- 1+(R*K)/(6*sm*(length(mr)-1)*(length(mk)-1))
  E <- outer(mr,mk,"*")/sm
  C <- sum((x-E)^2/E)

  cat("\n");
  cat("Coefficient:","\n");
  cat("\n");
  cat("fi-Yule'a: ",formatC(sqrt(C/sm)),     ",  T-Czupurow: ",formatC(sqrt((C)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) )),"\n");
  cat("C-Pearson: ",formatC(sqrt(C/(sm+C))), ",  V-Cramer:   ",formatC(sqrt((C)/(sm*min((ncol(x)-1),(nrow(x)-1 ) )) )),"\n");
  cat("\n");

}
