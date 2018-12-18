\encoding{utf8}
\name{coeff}
\alias{coeff}
\title{
The contingency coefficient
}
\description{
Computes contingency coefficient.
}
\usage{
coeff(x)
}
\arguments{
  \item{x}{contingency table.}
}
\details{
Contingency coefficient: C-Pearson, T-Czupurow, Fi-Yule'a, V-Cramer.
}

\seealso{
\code{\link{chisq.test}}
}
\examples{
# data:
m= matrix(c(23,32,45,26),2,2)

# contingency coefficient:
coeff(m)
}
\keyword{coeff}
