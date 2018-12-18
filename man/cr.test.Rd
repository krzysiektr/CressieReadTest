\encoding{utf8}
\name{cr.test}
\alias{cr.test}

\title{
The D-squared Cressie-Read test
}
\description{
Computes D-squared Cressie-Read test.
}
\usage{
cr.test(x, lambda = 2/3)
}

\arguments{
  \item{x}{contingency table.}
  \item{lambda}{parametr lambda.}
}

\details{
Statistics Cressie-Read:

D = 2/(lambda*(lambda+1)) * sum( O * ((O/E)^lambda - 1) )

> O - the number of observed

> E - the number of expected

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
\code{\link{chisq.test}}   \code{\link{fisher.test}}
}

\examples{
# data:
m= matrix(c(23,32,45,26),2,2)

# Cressie-Read:
cr.test(m)

# Pearson:
cr.test(m,lambda=1)

# Likelihood Ratio:
cr.test(m,lambda=1e-05)

# Freeman-Tukey's:
cr.test(m,lambda=-1/2)

# Neyman's:
cr.test(m,lambda=-2)

# Kullback-Leibler:
cr.test(m,lambda=-0.99999)
}
\keyword{cr.test}


