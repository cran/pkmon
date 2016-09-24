################################################################################
# Empirical estimator
pEmp=function(X){
	supp=seq(min(X),max(X),by=1);
	count=rep(0,length(supp));
	for(i in 1:length(supp)){
		count[i]=sum(X==supp[i]);
	}
	ptilde=count/length(X);
	return(list(support=supp,count=count,freq=ptilde));
}

################################################################################
# Spline
dSpline=function(supp, k){
  # calcul of the coefficients of a spline 
  #
  # INPUT
  # supp : integer for the maximum of the support of the spline
  # k : integer for the degree of the k-monotony of the spline
  #
  # OUTPUT
  # vector with length supp+1
  
  # Checking input
  if(abs(supp - round(supp)) >= .Machine$double.eps^0.5 | supp<0){
    stop("supp should be a positive integer")
  }
  
  if(abs(k - round(k)) >= .Machine$double.eps^0.5 | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
	Q=BaseNorm(k, supp);	
	return(Q[, supp+1]);
}

rSpline=function(n=1, supp, k){
  # random generation for a spline
  #
  # INPUT
  # supp : integer for the maximum of the support of the spline
  # n : the size of the sample
  # k : integer for the degree of the k-monotony of the spline
  #
  # OUTPUT
  # vector of length n
  
  # Checking input
  if(abs(supp - round(supp)) >= .Machine$double.eps^0.5 | supp<0){
    stop("supp should be a positive integer")
  }
  
  if(abs(k - round(k)) >= .Machine$double.eps^0.5 | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  
	p=dSpline(supp, k=k);
	return(sample(0:supp, size=n, replace=TRUE, prob=p));
}

################################################################################
# Mixture of splines
dmixSpline=function(supp, k, prob){
  # calcul of the coefficients of a mixture of splines
  #
  # INPUT
  # supp : vector containing the support of the mixture
  # prob : vector of the same length as supp, containing the
  #      probabilities of each spline of the mixture
  # k : integer for the degree of the k-monotony of the splines
  #
  # OUTPUT
  # vector with length max(supp)+1
  
  # Checking input
  if(abs(supp - round(supp)) >= .Machine$double.eps^0.5 | supp<0){
    stop("supp should be a positive integer")
  }
  
  if(abs(k - round(k)) >= .Machine$double.eps^0.5 | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  if(!is.numeric(prob)){
    stop("prob should be a numeric vector.")
  }
  
	L=max(supp);
	proba=prob/sum(prob);
	Q=BaseNorm(k, L)[, supp+1];
	return(apply(diag(prob)%*%t(Q), 2, sum));
}

rmixSpline=function(n=1, supp, k, prob){
  # random generation for a mixture of splines
  #
  # INPUT
  # supp : vector containing the support of the mixture
  # prob : vector of the same length as supp, containing the
  #      probabilities of each spline of the mixture
  # n : the size of the sample
  # k : integer for the degree of the k-monotony of the splines
  #
  # OUTPUT
  # vector of length n
  
  # Checking input
  if(abs(supp - round(supp)) >= .Machine$double.eps^0.5 | supp<0){
    stop("supp should be a positive integer")
  }
  
  if(abs(k - round(k)) >= .Machine$double.eps^0.5 | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  if(!is.numeric(prob)){
    stop("prob should be a numeric vector.")
  }
  
	if(length(supp)!=1){
		prob=prob/sum(prob);
		j=sample(supp, size=n, replace=TRUE, prob=prob);
		X=unlist(lapply(j, rSpline, n=1, k=k));
	} else {
		X=rSpline(n=n, k=k);
	}
	return(X);
}
