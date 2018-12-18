cr.test=function(x, lambda=2/3) 
{
  DNAME <- deparse(substitute(x))
  mr <- apply(x,1,sum); mk <- apply(x,2,sum); sm <- sum(mr)
  E <- outer(mr,mk,"*")/sm
  Df <- (nrow(x)-1)*(ncol(x)-1)
  D <- (2/(lambda*(lambda+1)))*sum(x*( (x/E)^(lambda)-1 ))
  pvalue <- pchisq(D,Df,lower.tail=FALSE)
  
  if (lambda==2/3) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda==1) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Pearson)", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda <= 1e-05 && lambda > 0 || lambda < 0 && lambda >= -1e-05) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (G-squared)", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda <= -0.99999 && lambda > -1 || lambda < -1 && lambda >= -1.00001 ) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Kullback-Leibler)", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda == -1/2) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Freeman-Tukey's)", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda == -2) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Neyman's)", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
}

cr.gof=function(x, lambda=2/3, p=rep(1/length(x), length(x)) ) 
{
  DNAME <- deparse(substitute(x))
  Df <- length(x)-1
  E= sum(x)*p
  D <- (2/(lambda*(lambda+1)))*sum(x*( (x/E)^(lambda)-1 ))
  pvalue <- pchisq(D,Df,lower.tail=FALSE)
  
  if (lambda==2/3 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda==1 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Pearson) for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda <= 1e-05 & lambda > 0  & sum(p)==1 | lambda < 0 && lambda >= -1e-05 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (G-squared) for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda <= -0.99999 & lambda > -1  & sum(p)==1 | lambda < -1 & lambda >= -1.00001 & sum(p)==1 ) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Kullback-Leibler) for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda == -1/2 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Freeman-Tukey's) for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(lambda == -2 & sum(p)==1) {
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test (Neyman's) for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(sum(p)==1){ 
    RVAL <- list(statistic = c(D = D, lambda = lambda, df = Df), p.value = pvalue, 
                 method = "D-squared Cressie-Read test for given probabilities", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  stop("error")
}

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


cgf.test=function(x, test="gw") 
{
  DNAME <- deparse(substitute(x))
  mr <- apply(x,1,sum); mk <- apply(x,2,sum); sm <- sum(mr)
  R <- (sum(1/mr))*sm-1; K <- (sum(1/mk))*sm-1
  W <- 1+(R*K)/(6*sm*(length(mr)-1)*(length(mk)-1))
  E <- outer(mr,mk,"*")/sm
  G <- 2*(sum(x*log(x/E)))
  P <- sum((x-E)^2/E)
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
  if (test=="gw") {
    RVAL <- list(statistic = c(Gw = Gw, df = Df), p.value = pvalueB, 
                 method = "G-squared Likelihood Ratio test with Williams correction", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="g") {
    RVAL <- list(statistic = c(G = G, df = Df), p.value = pvalueA, 
                 method = "G-squared Likelihood Ratio test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="p") {
    RVAL <- list(statistic = c(P = P, df = Df), p.value = pvalueP, 
                 method = "X-squared Pearson test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="ft") {
    RVAL <- list(statistic = c(FT = FT, df = Df), p.value = pvalueFT, 
                 method = "F-squared Freeman-Tukey's test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="kl") {
    RVAL <- list(statistic = c(KL = KL, df = Df), p.value = pvalueKL, 
                 method = "G-squared Kullback-Leibler test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  else if(test=="n") {
    RVAL <- list(statistic = c(N = N, df = Df), p.value = pvalueN, 
                 method = "X-squared Neyman's test", 
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }
}


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



