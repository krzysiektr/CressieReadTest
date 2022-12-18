#' @title The contingency coefficient
#'
#' @description
#' Computes contingency coefficient: Q_Yulea, Phi, C-Pearson, C_adj, V-Cramer, T-Tschuprow for 2x2 table or
#' C-Pearson, C_adj, V-Cramer, T-Tschuprow for rxc table where r>2 or c>0.
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
  wp <- min(dim(x))
  C <- sum((x-E)^2/E)
  RES <- c((x[1,1]*x[2,2]-x[1,2]*x[2,1])/(x[1,1]*x[2,2]+x[1,2]*x[2,1]),
           sqrt(C/sm),
           sqrt(C/(sm+C)),
           (sqrt(C/(sm+C)))/sqrt((wp-1)/wp),
           sqrt(C/(sm*(wp-1)) ),
           sqrt((C)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) ) )
  names(RES) <- c("Q_Yulea","Phi","C-Pearson","C_adj","V-Cramer","T-Tschuprow")
  if (dim(x)[1] == 2 & dim(x)[2] == 2) {
    return(RES)
  }
  else {
    return(RES[-c(1,2)])
  }
}
