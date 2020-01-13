\name{PraD}
\alias{PraD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Discrepancy with respect to inequality constraints}
\description{
Internal function that, given a vector of joint probabilities from an \emph{I}x\emph{J} table (vectorized by row) and the RC model specification in list Model, computes vector of discrepancies and matrix of its derivatives with respect to the canonical parameters.
}
\usage{
PraD(the, Model, der = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{the}{vector of canonical parameters}
  \item{Model}{list of model components}
  \item{der}{to require the derivative (optional)}
}

\value{
\item{hdis}{vector of discrepancies}
\item{Hdis}{matrix of derivatives of discrepancies with respect to the canonical parameter (optional)}
}

\author{
Francesco Bartolucci, Antonio Forcina
}
