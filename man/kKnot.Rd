\name{kKnot}
\alias{kKnot}
\title{
	k-Knot	
}
\description{
	k-Knots of a discrete function.
}
\usage{
kKnot(p, k)
}
\arguments{
  \item{p}{
  Vector
}
  \item{k}{
  Degree of the knots
}
}
\details{
	An integer i is a k-knot of p if \eqn{\Delta^k p(i) >0}, where \eqn{\Delta^k} is the k-th Laplacian of the sequence p.
}
\value{
	Vector with the k-knots of p.
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
\link{Delta}
}
\examples{
p=dmixSpline(c(5, 10, 20), k=3, c(0.5, 0.25, 0.25))
knots=kKnot(p, 3) 
}
