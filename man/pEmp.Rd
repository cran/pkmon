\name{pEmp}
\alias{pEmp}
\title{
	Empirical estimator of a discrete function
}
\description{
	Empirical estimator of a discrete function
}
\usage{
pEmp(X)
}
\arguments{
  \item{X}{
  A random sample from a discrete probability.
}
}
\details{
The empirical estimator is defined as \eqn{p(j)=\Sigma_{i=1}^n \bold{1}_{x_j=j}}.
}
\value{
	\item{support}{The points of the support of the estimator}
	\item{count}{The counts of the sample}
	\item{freq}{The normalized counts}
}
\author{Jade Giguelay
}

\examples{
x=rpois(100, lambda=0.3)
ptild=pEmp(x)
}
