\name{extRC}
\alias{extRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extended RC model}
\description{
Main function that fits extended RC models based on different types of aggregation (continuation, local, global) and different divergence functions defined by a suitable value of lambda.
}
\usage{
extRC(N, mod, k, la, marg.cons = c("free","equal","shift"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{observed contingency table}
  \item{mod}{vector indicating the types of aggregation for row and column variables ("c" for continuation, "l" for local, "g" for global)}
  \item{k}{rank required for the matrix of interaction parameters}
  \item{la}{value of lambda parameter}
  \item{marg.cons}{type of constraint on the marginal distributions}
}

\value{
\item{la}{vector of lambda values (when a vector is in input)}
\item{dev}{deviance of the fitted model (when only one lambda value is in input) or vector of deviances (when a vector of lambda values is in input)}
\item{df}{degrees of freedom (when  only one lambda value is in input)}
\item{it}{number of iterations (when  only one lambda value is in input)}
\item{dis}{final discrepancy (when  only one lambda value is in input)}
\item{pj}{vector of joint probabilities under the fitted model (when  only one lambda value is in input)}
\item{eta}{full vector of marginal parameters (when only one lambda value is in input)}
\item{etaX}{vector of row marginal parameters (when only one lambda value is in input)}
\item{etaY}{vector of column marginal parameters (whenonly one lambda value is in input)}
\item{Eta}{matrix of association parameters (when only one lambda value is in input)}
\item{la}{vector of lambda values (when more lambda values are in input)}
\item{dev}{vector of deviance values (when more lambda values are in input)}
}

\author{Francesco Bartolucci, Antonio Forcina}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# load data
data(mobility)

# for a single value of lambda, fit model with constraints of rank 1 on
# local-local logits and without constraints on the marginal distributions
out = extRC(mobility,mod=c("l","l"),k=1,la=0.6)
summary(out)

# for a single value of lambda, fit model with constraints of rank 1 on
# local-local logits and under constrain of equal marginal distributions
out = extRC(mobility,mod=c("l","l"),k=1,la=0.6,marg.cons="equal")
summary(out)

# for a single value of lambda, fit model with constraints of rank 2 on
# global-global logits and under constraint that marginal distributions
# are equal up to a constant shift
out = extRC(mobility,mod=c("g","g"),k=2,la=0.6,marg.cons="shift")
summary(out)

# fit model for a vector of lambdas
la = seq(-1.8,0.6,length.out=10)
out1 = extRC(mobility,mod=c("l","l"),k=1,la=la)
plot(out1)
}
