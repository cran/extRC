\name{dfm}
\alias{dfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First difference matrix
}
\description{
It creates a matrix of first differences of order \emph{k}. 
}
\usage{
dfm(k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{size of the matrix}
}
\value{
\item{D}{first difference matrix}
}

\author{Francesco Bartolucci, Antonio Forcina}

\examples{
D = dfm(5)
x = runif(5)
(D\%*\%x)
}
