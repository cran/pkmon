####################################################################################################
##Auxiliary functions
####################################################################################################


Base=function(k,J){
  # calculation of the matrix of (J+1)-th splines for
  # j varying from 0 to J
  #
  # INPUT
  # J : scalar
  # k : integer for the degree of the k-monotony
  #
  # OUTPUT
  # Q : matrix with J+1 rows and J+1 columns
  # Q(i,j)=Q_j^k(i-1)
  Q=matrix(1,J+1,J+1)
  for (j in 1:J){
    for (i in 0:(j-1)){
      Q[j-i,j]=Q[j-i+1,j]*(i+k)/(i+1)
    }
    if (j<J) {for (i in (j+1):J){
      Q[i+1,j]=0
    }
    }
  }
  dimnames(Q)=list(paste("i=",0:J,sep=""),paste("j=",0:J,sep=""))
  Q[,2:(J+1)]<-Q[,1:J]
  Q[,1]<-rep(0,J+1)
  Q[1,1]<-1
  return(Q)
}

BaseNorm=function(k,J){
  # Checking input
  if(J!=floor(J) | J<0){
    stop("J should be a positive integer")
  }
  
  if(k!=floor(k) | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  #Normalized version of Base
  Q=Base(k,J)
  R=apply(Q,2,sum)
  for (j in (1:(J+1))){
    Q[,j]=Q[,j]/R[j]
  }
  return(Q)
}

Delta=function(k,L,p){
  # Checking input
  if(!is.numeric(p)){
    stop("p should be a numeric vector.")
  }
  
  if(L!=floor(L) | L<0){
    stop("L should be a positive integer")
  }
  
  if(k!=floor(k) | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  #Return the j-th Laplacians (-1)^j*delta^j (p(l)) of vector p for j in 1:k and l in 0:L
  M=matrix(0,k,L+k)
  q=rep(0,L+k+1)
  n=length(p)
  q[1:n]=p
  for(l in 1:(L+k)){
    M[1,l]=q[l]-q[l+1]
  }
  for(j in 2:k){
    for(l in 1:(L+k-j+1)){
      M[j,l]=M[j-1,l]-M[j-1,l+1]
      if(abs(M[j,l])<10^(-15)){
        M[j,l]=0
      }
    }
  }
  M=M[1:k,1:(L+1)]
  dimnames(M) <- list(paste("j=",1:k,sep=""),paste("l=",0:L,sep=""))
  return(M)
}


Beta <- function(p,ptild){
  #Return the coefficient beta(p)=<p,p-ptild>
  j <- min(c(length(p),length(ptild)))
  Beta=sum(p**2) - sum(p[1:j]*ptild[1:j])
  return(Beta)
}


MasseBase <- function(k,j){
  #Return the mass of Q_j^k
  Q=Base(k,j)
  masse=apply(Q,2,sum)
  return(masse)
}


kKnot <- function(p,k){
  # Checking input
  if(!is.numeric(p)){
    stop("p should be a numeric vector.")
  }
  
  if(k!=floor(k) | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  #Return the k-knots of a vector p
  H=Delta(k,length(p),p)[k,]
  noeud <- which(H>0)-1
  return(noeud)
}
