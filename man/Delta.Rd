\name{Delta}
\alias{Delta}
\title{
	Discrete laplacian
}
\description{
	Computes the laplacians of a discrete function
}
\usage{
Delta(k, L, p)
}
\arguments{
  \item{k}{
  Maximum order of the laplacian
}
  \item{L}{
  Support of the function
}
  \item{p}{
  Discrete function represented as a vector
}
}
\value{
	Returns a matrix with the laplacians \eqn{(-1)^j\Delta^j (p(l))} of vector \eqn{p} for \eqn{j} in \eqn{1,\ldots,k} and \eqn{l} in \eqn{0,\ldots,L}. 
}

\references{
  Knopp K. (1925), <DOI:10.1007/BF01479598>
  Mehrfach monotone Zahlenfolgen,
  \emph{Mathematische Zeitschrift}, \bold{22}, 75--85

  Giguelay, J., (2016),
  Estimation of a discrete distribution under k-monotony constraint,
  \emph{in revision}, (arXiv:1608.06541)
}

\author{Jade Giguelay
	\email{jade.giguelay@ens-paris-saclay.fr}
	\url{http://maiage.jouy.inra.fr/jgiguelay}
}
\seealso{
\link{kKnot}
}
\examples{
p=dSpline(k=3, supp=20)
M=Delta(3, 20, p)
}
