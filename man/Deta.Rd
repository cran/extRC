\name{Deta}
\alias{Deta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of marginal parameters
}
\description{
Given a vector of canonical parameters coding distribution for an \emph{I}x\emph{J} contingency table and the RC model specification in list \code{Model}, it computes vector of marginal and joint parameters and matrix of its derivatives with respect the canonical parameters.
}
\usage{
Deta(the, Model, der = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{the}{vector of canonical parameters}
  \item{Model}{list specifying all model components}
  \item{der}{to require derivative computation (optional)}
}

\value{
\item{eta}{vector of marginal parameters}
\item{Der}{derivative matrix with respect to canonical parameters}
}

\author{Francesco Bartolucci, Antonio Forcina}
