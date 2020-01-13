\name{cuby}
\alias{cuby}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Step length}
\description{
Internal function that computes step length of the estimation algorithm in \code{\link{extRC}} by fitting a cubic polynomial.
}
\usage{
cuby(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{vector of likelihood values at different step lengths}
}

\value{
\item{comp1 }{optimal length}
}

\author{Francesco Bartolucci, Antonio Forcina}