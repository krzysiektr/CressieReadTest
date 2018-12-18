\encoding{utf8}
\name{cgf.test}
\alias{cgf.test}
\title{
The X-squared, G-squared, F-squared.
}
\description{
calculates statistics X-squared, G-squared, F-squared.
}
\usage{
cgf.test(x, test = "gw")
}
\arguments{
  \item{x}{contingency table.}
  \item{test}{"gw","g", "p", "n", "ft", "kl"}
}
\details{
test without correction test="g": 
G=2*(sum(O*log(O/E)))

Williams correction for test="gw" (the default):
Gw=G/W, where correct - W

test="p" - X-squared Pearson test: 
P=sum((O-E)^2/E)

test="n" - Neyman's test: 
N=sum((O-E)^2/O)

test="ft" - Freeman-Tukey's test: 
FT=4*sum((sqrt(O)-sqrt(E))^2)

test="kl" - Kullback-Leibler test: 
KL= 2*sum(E*log(E/O))
}
\value{
\item{statistic[1]}{the test statistic.}
\item{statistic[2]}{the degrees of freedom.}
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

# LR test with correction:
cgf.test(m)

# LR test without correction:
cgf.test(m, test="g")

# test Pearson:
cgf.test(m, test="p")
}
\keyword{cgf.test}
