\name{MatIn}
\alias{MatIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Aggregation matrices}
\description{
Computation of aggregation matrices for generalized interactions that are used in code{\link{extRC}} to estimate extended RC models.
}
\usage{
MatIn(lev, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lev}{vector number of rows and columns}
  \item{mod}{type of logit for each dimension}
}

\value{
\item{R0}{aggregation matrix for the row margin upper level}
\item{R1}{aggregation matrix for the row margin lower level}
\item{C0}{aggregation matrix for the column margin upper level}
\item{C1}{aggregation matrix for the column margin lower level}
\item{J00}{aggregation matrix for the left upper quadrant}
\item{J01}{aggregation matrix for the rigth upper quadrant}
\item{J10}{aggregation matrix for the left lower quadrant}
\item{J11}{aggregation matrix for the right lower quadrant}
}

\author{Francesco Bartolucci, Antonio Forcina}
