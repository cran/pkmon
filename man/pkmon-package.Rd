\name{pkmon-package}
\alias{pkmon-package}
\alias{pkmon}
\docType{package}
\title{
Least-squares estimator under k-monotony constraint for discrete functions
}
\description{
Description: This package implements two least-squares estimators under k-monotony constraint using a method based on the Support Reduction Algorithm from Groeneboom et al (2008). The first one is a projection estimator on the set of k-monotone discrete functions. The second one is a projection on the set of k-monotone discrete probabilities. This package provides functions to generate samples from the spline basis from Lefevre and Loisel (2013), and from mixtures of splines.
}
\details{
\tabular{ll}{
Package: \tab pkmon\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2016-07-28\cr
License: \tab CC BY 4.0\cr
}
}
\author{Jade Giguelay
	\email{jade.giguelay@ens-paris-saclay.fr}
	\url{http://maiage.jouy.inra.fr/jgiguelay}
}
\seealso{
\link{pMonotone}, \link{fMonotone}
}
\references{Giguelay J. (2016) 
  Estimation of a discrete distribution under k-monotony constraint,
  \emph{in revision}, (arXiv:1608.06541)


	Lefevre C., Loisel S. (2013) <DOI:10.1239/jap/1378401239>
	On multiply monotone distributions, continuous or discrete, with applications,
	\emph{Journal of Applied Probability}, \bold{50}, 827--847.

	Groeneboom P., Jongbloed G. Wellner J. A. (2008) <DOI:10.1111/j.1467-9469.2007.00588.x>
	The Support Reduction Algorithm for Computing Non-Parametric Function Estimates in Mixture Models,
	\emph{Scandinavian Journal of Statistics}, \bold{35}, 385--399
}
\keyword{ package }
\examples{

####################
# Example 1
# one triangular function T_j=Q_j^2, for j=supp=20 and for k=2 and k=3
n=30;
k1=2;
k2=3;
l=2;
supp=20;
p=dSpline(supp, k=l);
X=rSpline(n=n, supp, k=l);
ptilde=pEmp(X);
phat1=pMonotone(ptilde$freq, k=k1);
phat2=pMonotone(ptilde$freq, k=k2);

x.limits=c(0, max(supp+1, phat1$Spi+1, phat2$Spi+1));
y.limits=range(p, ptilde$freq, phat1$p, phat2$p);

plot(NULL, xlim=x.limits, ylim=y.limits, xlab="Counts", ylab="Frequencies");
points(0:supp, p, pch=16, col=1, lwd=2); 
points(ptilde$supp, ptilde$freq, pch=4, col=2, lwd=2);
points(0:max(phat1$Spi), phat1$p, pch=8, col=3, lwd=2);
points(0:max(phat2$Spi), phat2$p, pch=2, col=4, lwd=2);
legend("topright", pch=c(16, 4, 8, 2), col=c(1, 2, 3, 4),
	legend=c("p", expression(tilde("p")), expression(hat("p")*" - k = 2"), 
	expression(hat("p")*" - k = 3")));

####################
# Example 2
# mixture of 3 splines Q_j^3 and for k=4 and k=3
n=30;
k1=4;
k2=3;
l=3;
supp=c(5, 10, 20);
prob=c(0.5, 0.3, 0.2);
p=dmixSpline(supp, k=l, prob=prob);
X=rmixSpline(n=n, supp, k=l, prob=prob);
ptilde=pEmp(X);
phat1=pMonotone(ptilde$freq, k=k1);
phat2=pMonotone(ptilde$freq, k=k2);

x.limits=c(0, max(supp+1, phat1$Spi+1, phat2$Spi+1));
y.limits=range(p, ptilde$freq, phat1$p, phat2$p);

plot(NULL, xlim=x.limits, ylim=y.limits, xlab="Counts", ylab="Frequencies");
points(0:max(supp), p, pch=16, col=1, lwd=2); 
points(ptilde$supp, ptilde$freq, pch=4, col=2, lwd=2);
points(0:max(phat1$Spi), phat1$p, pch=8, col=3, lwd=2);
points(0:max(phat2$Spi), phat2$p, pch=2, col=4, lwd=2);
legend("topright", pch=c(16, 4, 8, 2), col=c(1, 2, 3, 4),
	legend=c("p", expression(tilde("p")), expression(hat(p)* " - k = 4"), 
	expression(hat(p)* " - k = 3")));

####################
# Example 3
# Poisson density
n=30;
k1=2;
k2=3;
supp=10;
p=dpois(0:supp, lambda=1);
X=rpois(n, lambda=1);
ptilde=pEmp(X);
phat1=pMonotone(ptilde$freq, k=k1);
phat2=pMonotone(ptilde$freq, k=k2);

x.limits=c(0, max(supp, phat1$Spi+1, phat2$Spi+1));
y.limits=range(p, ptilde$freq, phat1$p, phat2$p);

plot(NULL, xlim=x.limits, ylim=y.limits, xlab="Counts", ylab="Frequencies");
points(0:max(supp), p, pch=16, col=1, lwd=2);
points(ptilde$supp, ptilde$freq, pch=4, col=2, lwd=2);
points(0:max(phat1$Spi), phat1$p, pch=8, col=3, lwd=2);
points(0:max(phat2$Spi), phat2$p, pch=2, col=4, lwd=2);
legend("topright", pch=c(16, 4, 8, 2), col=c(1, 2, 3, 4),
       legend=c("p", expression(tilde("p")), expression(hat(p)* " - k = 2"),
                expression(hat(p)* " - k = 3")));

\dontrun{
####################
# Simulation for comparing ptilde and pHat (p is 3-monotone, k=3)
#
#OUTPUT
#
# cvge : percentage of non-convergence of the algorithm
# r.emp : L2-risk for the empirical estimator
# r.Hat : L2-risk for the estimator under k-monotony constraint


nSim=500;

n=30;
k=3;
l=3;
supp=20;
p=dSpline(supp, k=l);

result <- matrix(nrow=nSim,ncol=3);
dimnames(result)[[2]] <- c("cvge","r.emp","r.Hat");
for (i in 1:nSim) {
  X=rSpline(n=n, supp, k=l);
  ptilde=pEmp(X);
  phat=pMonotone(ptilde$freq, k=k);
  m <- max(supp+1,length(ptilde$freq)+1,phat$Spi+1)
  pV=c(p,rep(0,m-length(p)))
  pHat=c(phat$p,rep(0,m-length(phat$p)))
  ptilde=c(ptilde$freq,rep(0,m-length(ptilde$freq)))
  result[i,] <- c(phat$cvge,sum((pV-ptilde)**2),
                  sum((pV-pHat)**2))
}

apply(result,2,mean)

#Example with set.seed(0)
#cvge       r.emp       r.Hat 
#0.000000000 0.030682552 0.004984899
}
}
