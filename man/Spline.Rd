\name{Spline}
\alias{rSpline}
\alias{dSpline}
\alias{rmixSpline}
\alias{dmixSpline}
\title{
Random generation and distribution function of k-monotone densities
}
\description{
Random generation and distribution function for the spline of the basis from Lefevre and Loisel (2013), and mixtures of splines.
}
\usage{
rSpline(n=1, supp, k)
dSpline(supp, k)
rmixSpline(n=1, supp, k,prob)
dmixSpline(supp, k, prob)
}
\arguments{
  \item{supp}{
  Support of the spline, or vector of the supports of the splines for the mixture of splines
}
  \item{n}{
  Number of random values to return
}
  \item{k}{
  Degree of monotony
}
	\item{prob}{
	Vector of probabilities for the mixture of splines
}
}
\details{
See \link{BaseNorm} for details on the spline basis.
}
\value{
	rSpline and rmixSpline generates random deviates from the splines and mixtures of splines.

	dSpline and dmixSpline gives the distribution function.
}
\references{
  Giguelay, J., (2016),
  Estimation of a discrete distribution under k-monotony constraint,
  \emph{in revision}, (arXiv:1608.06541)


	Lefevre C., Loisel S. (2013) <DOI:10.1239/jap/1378401239>
	On multiply monotone distributions, continuous or discrete, with applications,
	\emph{Journal of Applied Probability}, \bold{50}, 827--847.
}
\author{Jade Giguelay
}
\seealso{
\link{pEmp}
}
\examples{
x=rSpline(n=100, 20, 3)
p=dSpline(20, 3)
xmix=rmixSpline(n=100, c(5, 20), 3, c(0.5, 0.5))
pmix=dmixSpline(c(5, 20), 3, c(0.5, 0.5))
par(mfrow=c(1, 2))
hist(x, freq=FALSE)
lines(p, col="red")
hist(xmix, freq=FALSE)
lines(pmix, col="red")
}
