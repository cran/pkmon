\name{BaseNorm}
\alias{BaseNorm}
\title{
Normalized spline basis
}
\description{
Computes the k-monotone discrete splines from Lefevre and Loisel (2013).
}
\usage{
BaseNorm(k, J)
}
\arguments{
  \item{k}{
  Degree of monotony
}
  \item{J}{
  maximum support of the splines
}
}
\value{
matrix \eqn{Q} with J+1 rows and J+1 columns with \eqn{Q(i,j)=Q_j^k(i-1)=C_{j-i+k-1}^{k-1}}, where \eqn{C} represents the binomial coefficient.
}
\references{Giguelay, J., (2016),
  Estimation of a discrete distribution under k-monotony constraint,
  \emph{in revision}, (arXiv:1608.06541)


		Lefevre C., Loisel S. (2013) <DOI:10.1239/jap/1378401239>
	On multiply monotone distributions, continuous or discrete, with applications,
	\emph{Journal of Applied Probability}, \bold{50}, 827--847.
}
\author{Jade Giguelay
}
\seealso{
\link{rSpline}, \link{dSpline}, \link{rmixSpline}, \link{dmixSpline}
}
\examples{
# Computing 3-monotone splines with maximum support 8
Q=BaseNorm(3, 8)
matplot(Q, type="l", main="3-monotone splines with maximum support 8");
}
