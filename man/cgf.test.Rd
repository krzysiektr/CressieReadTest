% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cgf.test.R
\name{cgf.test}
\alias{cgf.test}
\title{The Cressie-Read test for contingency table.}
\usage{
cgf.test(x, test="pc")
}
\arguments{
\item{x}{contingency table.}

\item{test}{"gw","g", "p", "n", "ft", "kl"}
}
\value{
An object of class "htest" containing the following components:
\describe{
\item{parameter}{the degrees of freedom}
\item{statistic}{statistic test value}
\item{p.value}{p-value}
         }
}
\description{
Calculates statistics X-squared, G-squared, F-squared for contingency table.
}
\examples{
m <- matrix(c(23,32,45,26),2,2)

# LR test with correction:
cgf.test(m, test="gw")

# LR test without correction:
cgf.test(m, test="g")

# test Pearson:
cgf.test(m, test="p")

# test Pearson with correction:
cgf.test(m)

}
\references{
Noel Cressie and Timothy R. C. Read (1984).
Multinomial Goodness-of-Fit Test. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 46, No. 3 (1984), 440-464.
}
\seealso{
\code{\link[vcd]{assocstats}}
}
\author{
Krzysztof Trajkowski
}
\keyword{cgf.test}
