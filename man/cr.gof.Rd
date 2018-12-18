\encoding{utf8}
\name{cr.gof}
\alias{cr.gof}

\title{
The D-squared Cressie-Read test for given probabilities
}
\description{
Computes D-squared Cressie-Read test for given probabilities.
}
\usage{
cr.gof(x, lambda = 2/3, p=rep(1/length(x), length(x)))
}

\arguments{
  \item{x}{vector.}
  \item{lambda}{parametr lambda.}
  \item{p}{probabilities}
}

\details{
Statistics Cressie-Read:

D = 2/(lambda*(lambda+1)) * sum( O * ((O/np)^lambda - 1) )

> O - the number of observed

> np - n*p

Lambda parameter to be different from 0, and -1.

> lambda=  1 (Pearson test)

> lambda= -2 (Neyman's modified Pearson test)

> lambda~  0, lambda is close to 0 (Likelihood Ratio test)

> lambda~ -1 lambda is close to -1 (Kullback-Leibler modified Likelihood Ratio test)

> lambda= -1/2 (Freeman-Tukey's test)
}

\value{
\item{statistic[1]}{the Cressie-Read statistic.}
\item{statistic[2]}{the lambda parametr.}
\item{statistic[3]}{the degrees of freedom.}
\item{p.value}{the p-value.}
}
\references{
Noel Cressie and Timothy R. C. Read (1984).
Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
}

\seealso{
\code{\link{chisq.test}}
}

\examples{
# data:
s=c(12,15,20,25,17,28)
p= 1/6
n= sum(s)*p
np= n*p
# Cressie-Read:
cr.gof(s)

# Pearson:
cr.gof(s,lambda=1)
}
\keyword{cr.gof}

