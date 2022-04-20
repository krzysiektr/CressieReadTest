#' @title All tests of the family Cressie-Read for discrete data.
#'
#' @description
#' All tests of the family Cressie-Read and contingency coefficient.
#'
#' @param x contingency table.
#' @param lambda parameter lambda test Cressie-Read.
#'
#' @usage allcr.test(x, lambda = 2/3)
#' @return The whole family of tests Cressie-Read and contingency coefficient: C-Pearson, T-Czupurow, Fi-Yule'a, V-Cramer.
#' @keywords allcr.test
#' @author
#' Krzysztof Trajkowski
#' @examples
#' # data:
#' m <- matrix(c(23,32,45,26),2,2)
#' # all test:
#' allcr.test(m)
#' @rdname allcr.test
#'
#' @references Noel Cressie and Timothy R. C. Read (1984).
#' Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
#'
#' @importFrom stats dgamma dpois dnbinom dnorm nlminb pchisq pgamma ppois pnbinom pnorm sd var
#'
#' @seealso \code{\link[vcd]{goodfit}}
#' @export
#'

allcr.test=function(x,lambda=2/3)
{
  DNAME <- deparse(substitute(x))
  mr=apply(x,1,sum); mk=apply(x,2,sum); sm=sum(mr)
  E= outer(mr,mk,"*")/sm
  R <- (sum(1/mr))*sm-1; K <- (sum(1/mk))*sm-1
  W <- 1+(R*K)/(6*sm*(length(mr)-1)*(length(mk)-1))
  Df=(nrow(x)-1)*(ncol(x)-1)
  CR1= sum((x-E)^2/E) # Pearson
  CR2= 2*sum(x*log(x/E)) # G2
  WIL=CR2/W
  CR3= (2/(lambda*(lambda+1)))*sum(x*( (x/E)^(lambda)-1 )) # CR
  CR4= 4*sum((sqrt(x)-sqrt(E))^2) # F-T
  CR5= 2*sum(E*log(E/x)) # Kullback-Leibler
  CR6= sum((x-E)^2/x) # N
  C1=sqrt(CR1/sm) # Y-Yule'a
  C2=sqrt(CR1/(sm+CR1)) # C-Pearson
  C3=sqrt((CR1)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) ) # V-Cramer
  C4= sqrt((CR1)/(sm*sqrt((ncol(x)-1)*(nrow(x)-1 ) )) ) # T-Czupurow
#    if(full) {
  return(list(
test=
    matrix(c(
      1,CR1,Df,pchisq(CR1,Df,lower.tail=FALSE),
      0,CR2,Df,pchisq(CR2,Df,lower.tail=FALSE),
      0,WIL,Df,pchisq(WIL,Df,lower.tail=FALSE),
      lambda,CR3,Df,pchisq(CR3,Df,lower.tail=FALSE),
      -0.5,CR4,Df,pchisq(CR4,Df,lower.tail=FALSE),
      -2,CR6,Df,pchisq(CR6,Df,lower.tail=FALSE),
      -1,CR5,Df,pchisq(CR5,Df,lower.tail=FALSE)
    ),nrow=7,ncol=4,byrow=T,
           dimnames=list(
             c("Pearson's Chi-squared","Log likelihood ratio","* Williams correction","Cressie-Read",
               "Freeman-Tukey's","Neyman's","Kullback-Leibler"),
             c("lambda","Statistic","df","p-val"))),
coefficient=
    matrix(c(
      C1,C2,C3,C4
    ),nrow=1,ncol=4,byrow=T,
           dimnames=list(
             c("Coefficient:"),
             c("Y-Yule'a","C-Pearson","V-Cramer","T-Czupurow")))
  ))
    }
