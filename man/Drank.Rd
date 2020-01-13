\name{Drank}
\alias{Drank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check matrix rank}
\description{
Given the row vectorized matrix, it computes the vector of discrepancies with respect to a certain rank and its derivative.}
\usage{
Drank(ga, lev, k, der = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ga}{row vectorized matrix of interaction}
  \item{lev}{vector of the number of row and column categories in the original table (the numbers of rows and columns of the input matrix must be increased by 1)}
  \item{k}{matrix rank}
  \item{der}{to require derivative}
}

\value{
  \item{fr}{vector of discrepancies with respect to the rank}
  \item{Dfr}{derivative of fr}
}

\references{
Bartolucci, F. and Forcina, A. (2002). Extended RC association models allowing for order restrictions and marginal modeling. \emph{Journal of the American Statistical Association}, \bold{97}, 1192-1199.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
A = matrix(rnorm(12),4) # matrix the rank of which must be checked
a = as.vector(t(A))
out = Drank(a,c(5,4),1,der=TRUE)
(out$fr)
(out$Dfr)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{algebra}% use one of  RShowDoc("KEYWORDS")
