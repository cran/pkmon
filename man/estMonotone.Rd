\name{estMonotone}
\alias{pMonotone}
\alias{fMonotone}
\title{
Estimators of discrete probabilities under k-monotony constraint
}
\description{
Estimators of discrete probabilities under k-monotony constraint. Estimation can be done on the set of k-monotone functions or on the set of k-monotone probabilities.
}
\usage{
pMonotone(ptild, t.zero = 1e-10, t.P = 1e-08, max.sn = NULL, k, verbose = FALSE)
fMonotone(ptild, t.zero = 1e-10, t.P = 1e-08, max.sn = NULL, k, verbose = FALSE)
}
\arguments{
  \item{ptild}{
  Empirical estimator
}
  \item{t.zero}{
  Threshold for the precision of the directionnal derivatives. (see OUTPUT below)
}
  \item{t.P}{
  Threshold for the precision on the stopping criterion. (see OUTPUT below)
}
  \item{max.sn}{
  The maximum support for the evaluation of the estimator
}
  \item{k}{
  Degree of monotony
}
  \item{verbose}{
  if TRUE, print for each iteration on the  maximum support : pi, Psi and sumP (see OUTPUT below)
}
}
\details{
The thresholds t.P and t.zero are used for the precision in the algorithm : in Step one (See REFERENCES below) the algorithm computes the directionnal derivatives of the current estimator and stops if all the directionnal derivarives are null that is to say if they are smaller than t.zero. In Step two (See REFERENCES below) the algorithm computes a stopping criterion and stops if and only if the stopping criterion is verified that is to say if some quantities are non-negative that is to say bigger than -t.P.
}
\value{
\item{cvge}{cvge = 0 if the criterion Psi decreases with the support of pi.
cvge = 1 if Psi increases.
cvge = 2 if maximum number of iterations reached}
\item{Spi}{Support of the positive measure pi at the last iteration}
\item{pi}{Values of the positive measure pi at the last iteration}
\item{p}{Values of pHat}
\item{Psi}{Scalar value of the criterion to be minimised}
\item{sumP}{\code{sum(pHat)} at convergence}
\item{history}{Data frame with components}
\item{L}{The maximum of the support of pi}
\item{Psi}{Value of the criterion for the value L}
\item{SumP}{Value of \code{sum(pHat)}}
}
\references{Giguelay, J., (2016),
  Estimation of a discrete distribution under k-monotony constraint,
  \emph{in revision}, (arXiv:1608.06541)

  Groeneboom P., Jongbloed G. Wellner J. A. (2008) <DOI:10.1111/j.1467-9469.2007.00588.x>
	The Support Reduction Algorithm for Computing Non-Parametric Function Estimates in Mixture Models,
	\emph{Scandinavian Journal of Statistics}, \bold{35}, 385--399
}
\author{Jade Giguelay
	\email{jade.giguelay@ens-paris-saclay.fr}
	\url{http://maiage.jouy.inra.fr/jgiguelay}
}
\seealso{
\link{pEmp}, \link{BaseNorm}
}
\examples{
x=rSpline(n=50, 20, k=4)
ptild=pEmp(x);
res=pMonotone(ptild$freq, k=4)
}
