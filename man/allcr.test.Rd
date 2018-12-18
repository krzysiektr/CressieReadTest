\encoding{utf8}
\name{allcr.test}
\alias{allcr.test}
\title{
All tests of the family Cressie-Read.
}
\description{
All tests of the family Cressie-Read and contingency coefficient.
}
\usage{
allcr.test(x, lambda = 2/3)
}
\arguments{
  \item{x}{contingency table.}
  \item{lambda}{parameter lambda test Cressie-Read.}
}
\details{
The whole family of tests Cressie-Read and contingency coefficient: C-Pearson, T-Czupurow, Fi-Yule'a, V-Cramer.
}
\value{
\item{allcr.test(x)$test}{results for tests.}
\item{allcr.test(x)$coeff}{results for coefficient.}
}

\references{
Noel Cressie and Timothy R. C. Read (1984).
Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
}

\examples{
# data:
m= matrix(c(23,32,45,26),2,2)

# all test:
allcr.test(m)
}
\keyword{allcr.test}
