\name{tril}
\alias{tril}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lower triangular matrix}

\description{Given a square matrix, it provides the lower triangular part, including the main diagonal.}

\usage{
tril(M)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{square matrix
}
}
\value{
\item{N}{transformed matrix}
}

\author{Francesco Bartolucci, Antonio Forcina
}

\examples{
M = matrix(1:9,3)
N = tril(M)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}% use one of  RShowDoc("KEYWORDS")