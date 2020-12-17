\name{ars}

\alias{ars}

\title{Adaptive Rejection Sampling}

\usage{
ars(g = dnorm, bounds = c(-Inf, Inf), n = 1000, initial = NULL)
}

\arguments{
\item{g}{function from which to sample; log(g) must be a log concave function}
\item{bounds}{length 2 vector defining the lower and upper bounds for the distribution g}
\item{n}{the number of observations needed to sample}
\item{initial}{the initial value of x abscissae which can be defined by users}
}


\value{
a length n vector of samples from function g 
}

\description{
\code{ars} generate n samples from distribution g.
}

\details{
Adaptive Rejection Sampling method for sampling from univarite and log-concave function.
}

\examples{
n <- 1000
x <- ars(dnorm, bounds = c(-Inf, Inf), n, initial = NULL)
x <- ars(dexp, bounds = c(1, Inf), n, initial = NULL)
x <- ars(dgamma, bounds = c(0, Inf), n, initial = NULL)
  
## Plot histogram to check
## hist(x, freq=FALSE)
}

\references{
W. R. Gilks and P. Wild (1992). "AAdaptive Rejection Sampling for Gibbs Sampling".
    \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics).} \strong{41} (2): 337-348.
}

\seealso{
Wikipedia \href{https://en.wikipedia.org/wiki/Rejection_sampling}{Rejection Sampling}.
}

\author{
Katherine Kempfert, Eric Chu, Jiahui Zhao
}