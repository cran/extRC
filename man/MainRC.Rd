\name{MainRC}
\alias{MainRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of extended RC models}
\description{
Internal function that implements the Aitchinson-Silvey algorithm to estimate extended RC models.
}
\usage{
MainRC(y, Model, the0 = NULL, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{row vectorized vector of frequencies of the contingency table}
  \item{Model}{list of model components}
  \item{the0}{initial vector of canonical parameters (optional)}
  \item{output}{to require full output (optional)}
}

\value{
\item{dev}{final deviance}
\item{df}{degrees of freedom}
\item{pj}{vector of joint probabilities}
\item{it}{number of iterations}
\item{dis}{final discrepancy}
}

\author{Francesco Bartolucci, Antonio Forcina}
