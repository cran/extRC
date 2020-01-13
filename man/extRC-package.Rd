\name{extRC-package}
\alias{extRC-package}

\docType{package}
\title{Overview of the Package extRC}
\description{
Estimation of extended RC models, which are formulate by constraining different types of association parameters to have a reduced rank.
}
\details{
The package contains functions for maximum likelihood (ML) estimation of an extended class of row-column (RC) association models for two-dimensional contingency tables, as described in Forcina (2019). These models are formulated by a condition of reduced rank on a matrix of extended association parameters, which are defined by choosing the logit type for the row and column variables among four different options and a transformation derived from Cressie and Read (1984). Among the available alternatives, it is possible to use log-odds ratio based on different types of aggregation of the joint probabilities. The class of models generalizes that proposed in Kateri and Papaioannou (1994), Bartolucci and Forcina (2002), and Espendiller (2017), and includes the original RC association models of Goodman (1979) and the correspondence analysis model, as formulated in Goodman (1981) and Gilula et al. (1988). Maximum likelihood estimation is based on an algorithm that is an adaptation of the Aitchison and Silvey (1958) algorithm for constrained ML estimation and is related to the algorithm described in Evans and Forcina (2013) for fitting constrained marginal models.

The main function in the package is \code{\link{extRC}} that provides an output that may be shown by usual R commands \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}}.
}

\author{
Francesco Bartolucci, Antonio Forcina

Maintainer: Francesco Bartolucci <francesco.bartolucci@unipg.it>
}
\references{

Aitchison, J. and Silvey (1958). Maximum-likelihood estimation of parameters subject to restraints. \emph{The Annals of Mathematical Statistics}, \bold{29}, 813-828.

Bartolucci, F. and Forcina, A. (2002). Extended RC association models allowing for order restrictions and marginal modeling. \emph{Journal of the American Statistical Association}, \bold{97}, 1192-1199.

Cressie, N. and Read, T.R., 1984. Multinomial goodness-of-fit tests. \emph{Journal of the Royal Statistical Society: Series B}, \bold{46}, 440-464.

Espendiller, M., 2017. \emph{Association in contingency tables}. Ph.D. thesis.

Evans, R.J. and Forcina, A. (2013). Two algorithms for fitting constrained marginal models. \emph{Computational statistics & Data analysis}, \bold{66}, 1-7.

Forcina (2019), An extended class of RC association models: definition and estimation, arXiv:1910.13848.

Gilula, Z., Krieger, A.M., and Ritov, Y., 1988. Ordinal association in contingency tables: Some interpretive aspects. \emph{Journal of the American Statistical Association}, \bold{83}, 540-545.

Goodman, L.A. (1979). Simple models for the analysis of association in cross-classifications having ordered categories. \emph{Journal of the American Statistical Association}, \bold{74}, 537-552.

Goodman, L.A., 1981. Association models and canonical correlation in the analysis of cross-classifications having ordered categories. \emph{Journal of the American Statistical Association}, \bold{76}, 320-334.

Kateri, M. and Papaioannou, T. (1994). \emph{f-divergence Association Models}. University of Ioannina.
}

\keyword{multivariate}

\examples{
# load data
data(mobility)

# fit model for a single la
out = extRC(mobility,mod=c("l","l"),k=1,la=0.6)
summary(out)
}
