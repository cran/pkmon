###############################################################################################
###Estimation on a finite support
###############################################################################################


pStep1 <- function(emp,S,pi,Q,t.zero){
  # calculates the directional derivatives D
  # and add a point to the support of pi if necessary
  #
  # INPUT
  # emp : vector of dimension L, whose M+1 first coordinates equal ptild
  # Spi : vector. Support of pi. max(Spi)=L
  # pi : vector with the same length as Spi. Values for the positive measure
  # Q : matrix with dim (L+1)*(L+1). Q=BaseNorm(k,L)
  # t.zero : non-zero minimum value
  #
  # OUTPUT
  # D : vector with L components. Directional derivatives of Psi
  # Addj : Scalar. Addj = min(D). If  Addj <0, the point Addj must be
  #        added to the support of pi
  
  m=dim(Q)[2] 
  D <- rep(0,m)
  # calculation of the directional derivatives
  # D[j], for j varying from 1 to L
  for(j in 0:(m-1)){
    d <- rep(0,j+1) 
    for(l in 0:(m-1)){
      Sl <- sum(pi*Q[l+1,S+1])
      d[l+1] <- (Q[l+1,j+1]-Sl)*(Sl-emp[l+1])
    }
    D[j+1] <- sum(d)                        
  }
  # zeroing values below the threshold
  if(sum(abs(D) < t.zero) != 0 ){ D[abs(D) < t.zero] <- 0}
  # If there exists at least one D[j] negative, then
  # add to the support Spi the value of j that minimizes D. Warning : it begins on 0.
  minD <- min(D)
  if (minD >=0) Addj=-1
  if (minD <0) Addj <- which.min(D)-1
  return( list(Addj=Addj, D=D) )
}


pStep2 <- function(ptild, S, pi, Addj,Q) {
  #  updates the support Sk if necessary
  #  S=Sk; pi=pik; Addj=St1$Addj;
  #
  # INPUT
  # ptild : ptild : vector of empirical probabilities
  # S : vector. Support of pi.
  # pi : vector with the same length as Spi. Values for the positive
  #      measure
  # Addj : point added to the support at Step1
  # Q : matrix with dim (L+1)*(L+1). Q=BaseNorm(k,L) 
  # k: integer for the degree of the k-monotony, k=2,3,4
  #
  # OUTPUT
  # code : scalar
  #  code = 0. The support Sk and the values pik has been successfully
  #            updated
  #
  #  code =1. It was not possible to find a new support.
  #
  # If code=0,
  #  NewS : the new support
  #  Newpi : values of the positive measure
  #
  # If code=1,
  #  Trace : nature of the problem
  
  cond=0
  Sk1=S
  pik1=pi
  M <- length(ptild)-1
  m<-max(c(Sk1,Addj)) 
  if (m>=M){
    emp<-rep(0,m+1)
    emp[1:(M+1)]<-ptild
  }
  if (m<M){
    emp<-ptild[1:(m+1)]
  }
  while (cond==0) {
    Sprime=c(Sk1,Addj)
    QS=Q[1:(m+1),Sprime+1]
    HS=QS%*%solve(t(QS)%*%QS)%*%t(QS)
    lambdaS=(sum(HS%*%emp)-1)/sum(HS)
    np <- length(c(Sk1, Addj)) 
    reg <- lm(emp-lambdaS~QS-1)$coef
    
    if (reg[np] < 0) {
      cond=1
      trace=paste("Step1 : added",Addj,
                  "Step2 : its coefficient is negative:",reg[indNeg])
      return(list(code=1,trace=trace))
    }
    if (reg[np] >= 0) {
      indNeg <- reg[1:(np-1)]<0
      # if the coefficients corresponding to Sk1 are all positive
      # then return to pEstim to continue the algorithm
      if (sum(indNeg) ==  0) {
        cond=1
        aux=list(code=0,NewS=c(Sk1,Addj),Newpi= reg)
        return(aux)
      }
      # if some coefficients corresponding to Sk1 are negative
      # then suppress one of them
      if (sum(indNeg) >= 1) {
        if (np==2) {
          Sk1 <- NULL
          pik1 <- NULL
        }
        if (np>2) {
          amoinsb  <- pik1-reg[-np]
          ii <- which(amoinsb>0)
          l <- which.min(pik1[ii]/amoinsb[ii])
          Sk1 <- Sk1[-ii[l]]
          pik1 <- pik1[-ii[l]]
        }
      }
      ## and calculate the projection on this new support
    } ### end : if(reg[np] >= 0){
  } ### end : while (cond==0){
}


pEstim <- function(L, ptild, t.zero,k){
  # calculation of pHat when the max(Supp(pi))=L
  # 
  # INPUT
  # L : scalar : max of the support of pi
  # ptild : vector of empirical probabilities
  # t.zero : non-zero minimum value
  # k: integer for the degree of the k-monotony, k=2,3,4
  #
  # OUTPUT
  # code : integer
  #        code=0 : the convergence is reached
  #        code=1 : the convergence is not reached
  #                 (no positive measure is found)
  #  Spi : vector of integers.  Support of pi
  #  pi : vector with the same length as Spi.
  #       Values of the positive measure
  #  p : vector with L components. Values of the convex function p
  #  SumP : scalar. Sum of the components of p
  #  Psi : scalar. value of the criterion Psi
  #  trace : if code=1, nature of the problem
  
  emp <- rep(0,L+1)
  M <- length(ptild)-1
  emp[1:(M+1)] <- ptild
  iter <- 0
  # Initial value for (pi, Spi) 
  
  Q <- BaseNorm(k,L)
  pi0 <-1
  Spi0 <- c(L)
  p0 <- pi0*Q[,L+1]
  Psi0 <- sum(p0^2)/2-sum((emp[1:length(p0)])*p0)
  SumP0 <- sum(p0)
  # The function pStep1 calculates the directional derivatives D
  # and add a point to the support of pi if necessary
  St1 <- pStep1(emp,Spi0,pi0,Q,t.zero)
  # If St1$Addj==-1 then the convergence is reached
  # when the max(Supp(pi))=L
  if (St1$Addj==-1) {
    return(list(code=0, Spi=Spi0, pi=pi0, p=p0, Psi=Psi0, SumP=SumP0))
  }
  # If St1$Addj !=-1
  # Add one point to the support, and update the support in order to
  # get a probability
  if (St1$Addj!=-1) {
    Sk = Spi0
    pik = pi0
    pk=p0
    code=1
    while(code==1 ) {
      # the function pStep2 updates the support Sk if necessary
      St2 <- pStep2(ptild, Sk, pik, St1$Addj, Q)
      # if St2$code=1, then return to the main function for continuing
      # with a new value of L
      if (St2$code >=1 ) {
        return(list(code=St2$code, trace=St2$trace,
                  Spi=Spi0, pi=pi0, p=p0, Psi=Psi0, SumP=SumP0))
      }
      # if not, then go to step 2 with the new positive measure
      if (St2$code==0) {
        Sk=St2$NewS
        pik=St2$Newpi
        St1 <- pStep1(emp,Sk,pik,Q,t.zero)
        # if St1$Addj=0, all the derivative are positive and the minimum
        # of Psi is reached (on the support with maximum smaller than L)
        if (St1$Addj==-1) {
          code=0
          p = matrix(pik,nrow=1,ncol=length(pik))%*%t(Q[1:(max(Sk)+1),Sk+1])
          Psi <- sum(p^2)/2-sum((emp[1:length(p)])*p)
          SumP <- sum(p)
          return(list(code=code, Spi=Sk, pi=pik, p=p, Psi=Psi, SumP=SumP))
        }
      }
    } ### end : while(code==1){
  } ### end : if(St1$Addj!=-1){
}


######################################################################################################
###The Stopping Criterion
######################################################################################################


pStoppingCriterion <- function(ptild,t.P,p,Spi,k){
  #The test for the Stopping Criterion
  #
  # INPUT
  # ptild : vector of empirical probabilities with M+1 components
  #         where M=max(Support(ptild))
  # t.zero : non-zero minimum value
  # p : a vector which is to be tested if it is equal to pHat (p = resL$p)
  # k: integer for the degree of the k-monotony, k=2,3,4
  #
  # OUTPUT
  #
  # a list of 5 booleens. p=pHat if all the booleens are true
   
  s <- max(length(ptild)-1,length(p)-1)
  aux=T
  B=(k==3)+(k==4)
  
  aux2=s+1
  if (B){
    Fp=matrix(0,k+1,aux2+1) #Matrix with the primitives of the function p
    Fptild=matrix(0,k+1,aux2+1) #Matrix with the primitives of the function ptild
    Fp[1,1:length(p)]=p
    Fptild[1,1:length(ptild)]=ptild
    for (j in 2:(k+1)){
        Fp[j,] <- cumsum(Fp[j-1,])
        Fptild[j,] <- cumsum(Fptild[j-1,])
    }
    C4 <- Beta(p,ptild)
    Diff1 <- Fp[k+1,] - Fptild[k+1,]-C4*MasseBase(k,aux2)
    C1 <- min(Diff1)
    if (C1 < -t.P) aux=FALSE
    Diff3 <- Fp[3,aux2+1]-Fptild[3,aux2+1]-C4*MasseBase(2,aux2)[aux2+1]
    Diffk <- Fp[k,aux2+1]-Fptild[k,aux2+1]-C4*MasseBase(k-1,aux2)[aux2+1]
    C2 <- min(Diff3,Diffk)
    if (C2 <  -t.P) aux=FALSE
    Diff <- abs(Diff1[Spi+1])
    C3 <- max(Diff)
    if (C3 > t.P) aux=FALSE
    if (C4 > t.P) aux=FALSE
  }
return(list(Cvge=aux,C1=C1,C2=C2,C3=C3,C4=C4))
}


#########################################################################################################
###Final function
#########################################################################################################


pMonotone <- function(ptild,  t.zero=1e-10, t.P=1e-8, max.sn=NULL, k, verbose=FALSE ){
  # Calculation of pHat with support 0:LHat
  #
  # INPUT
  # ptild : vector of empirical probabilities with M+1 components
  #         where M=max(Support(ptild))
  # t.zero : non-zero minimum value
  # t.P : threshold for testing that pHat  sums to one
  # max.sn : maximum of the support of pHat
  #          By default max.sn=2*k*M
  # k: integer for the degree of the k-monotony, k=2,3,4
  # verbose : if TRUE, print for each iteration
  #           on the  maximum support : pi, Psi and sumP (see OUTPUT below)
  
  
  #
  # OUTPUT
  #   res : list with components
  #   cvge : cvge = 0 if the criterion Psi decreases with the support of pi
  #          cvge = 1 if Psi increases
  #          cvge = 2 if maximum number of iterations reached
  #   Spi : support of the positive measure pi at the last iteration
  #   pi : values of the positive measure pi at the last iteration
  #   p : values of pHat
  #   Psi : scaler value of the criterion to be minimised
  #   sumP : scalar = sum(pHat) at convergence
  #   history : data frame with components
  #                   L : the maximum of the support of pi
  #                   Psi : value of the criterion for the value L
  #                   SumP : value of sum(pHat)
  
  # Checking input
  if(!is.numeric(ptild)){
    stop("ptild should be a numeric vector.")
  }
  
  if(k!=floor(k) | k<2){
    stop("k should be an integer strictly larger than 1")
  }
  
  M <- length(ptild)
  if (is.null(max.sn)) max.sn=2*k*M
  L1=M+1
  Psi <- NULL
  SumP <- NULL
 
  C1 <- NULL
  C2 <- NULL
  C3 <- NULL
  C4 <- NULL
    
  code <- NULL
  i=0
  valL <- round(seq(L1,max.sn,length=20))
  suiteL <- NULL
 
  for (L in valL) {
    i=i+1
    suiteL <- c(suiteL,L)
    resL <- pEstim(L, ptild, t.zero,k)
    if ((k==3)|(k==4)) Cvge <- pStoppingCriterion(ptild,t.P,resL$p,resL$Spi,k)
     
    if (verbose==TRUE) {
      print(paste("maximum support =",L))
      print(resL$pi[order(resL$Spi)])
      print(paste("Crit Psi=",resL$Psi,"Sum Proba=",resL$SumP))
      if (k>2) print(paste("C1=",Cvge$C1,"C2=",Cvge$C2,"C3=",Cvge$C3,"C4"=Cvge$C4))
    }
    
    if (resL$code>0) {
      code=c(code,resL$code)
      Psi <- c(Psi,resL$Psi)
      SumP <- c(SumP,resL$SumP )
      
      if ((k==3)|(k==4)) {
        C1 <-  c(C1, Cvge$C1)
        C2 <-  c(C2, Cvge$C2)
        C3 <-  c(C3, Cvge$C3)
        C4 <-  c(C4, Cvge$C4)
      }
     
      if (L==tail(valL,1)) {
             trace=paste("The maximum support, max.sn=",
                    tail(valL,1)," is reached. The convergence is not reached with Code=",  
                    resL$code)
       
        if (k==2)  history <- data.frame(L=suiteL,code=code,Psi=Psi,SumP=SumP)
        if ((k==3)|(k==4))   { history <-
          data.frame(L=suiteL,code=code,Psi=Psi,SumP=SumP,C1=C1,C2=C2,C3=C3,C4=C4)}
      
        pi <- resL$pi
        names(pi) <- paste("j=",resL$Spi,sep="")
        pi <- pi[order(resL$Spi)]
        Spi <- sort(resL$Spi)
        
        if (k==2) res <- list(cvge=2,trace=trace,history=history,
                              Spi=Spi,pi=pi,p=resL$p, Psi=resL$Psi,SumP=resL$SumP)
        if ((k==3)|(k==4)) res <- list(cvge=2,trace=trace,history=history,
                                       Spi=Spi,pi=pi,p=resL$p, Psi=resL$Psi, SumP=resL$SumP, Cvge$C1,Cvge$C2,Cvge$C3,Cvge$C4)
        return(res)
      }
    }
    ##### end if(resL$code>0)  
    
    if (resL$code ==0) {
      code=c(code,resL$code)
      Psi <- c(Psi,resL$Psi)
      SumP <- c(SumP,resL$SumP)
      
      if ((k==3)|(k==4)) {
        C1 <- c(C1,Cvge$C1)
        C2 <- c(C2,Cvge$C2)
        C3 <- c(C3,Cvge$C3)
        C4 <- c(C4,Cvge$C4)
      }
     
      if (i==1) {
        if ((k==2)&(abs(resL$SumP-1)<t.P)) {
          history <- data.frame(L=L1,code=code,Psi=Psi,SumP=SumP)
          trace="Convergence is reached"
          res <- list(cvge=0,trace=trace,history=history,
                      Spi=resL$Spi,pi=resL$pi,p=resL$p,
                      Psi=resL$Psi,SumP=resL$SumP)
          return(res)
        }
        
        if ((k==3)|(k==4)) {
          if (Cvge$Cvge) {
            history <- data.frame(L=L1,code=code,Psi=Psi,SumP=SumP,C1=C1,C2=C2,C3=C3,C4=C4)
            trace="Convergence is reached"
            res <- list(cvge=0,trace=trace,history=history,
                        Spi=resL$Spi,pi=resL$pi,p=resL$p,
                        Psi=resL$Psi,SumP=resL$SumP,C1=Cvge$C1,C2=Cvge$C2,C3=Cvge$C3,C4=Cvge$C4)
            return(res)
          }
        }
        
      }  #### end if (i==1) 
      
      if (i>1) {
        
        if ((k==2)&(abs(resL$SumP-1)<t.P)) {
          history <- data.frame(L=suiteL,code=code,Psi=Psi,SumP=SumP)
          trace="Convergence is reached"
          pi <- resL$pi
          names(pi) <- paste("j=",resL$Spi,sep="")
          pi <- pi[order(resL$Spi)]
          Spi <- sort(resL$Spi)
          res <- list(cvge=0,trace=trace,history=history,
                      Spi=Spi,pi=pi,p=resL$p,
                      Psi=resL$Psi,SumP=resL$SumP)
          return(res)
        }
        
        if ((k==3)|(k==4)) {
          if (Cvge$Cvge) {
            history <- data.frame(L=suiteL,code=code,Psi=Psi,SumP=SumP)
            trace="Convergence is reached"
            pi <- resL$pi
            names(pi) <- paste("j=",resL$Spi,sep="")
            pi <- pi[order(resL$Spi)]
            Spi <- sort(resL$Spi)
            res <- list(cvge=0,trace=trace,history=history,
                        Spi=Spi,pi=pi,p=resL$p,
                        Psi=resL$Psi,SumP=resL$SumP,C1=Cvge$C1,C2=Cvge$C2,C3=Cvge$C3,C4=Cvge$C4)
            return(res)
          }
        }
      } #### end : if (i>1) {
       
      
      if (k==2) history=data.frame(L=suiteL,code=code,Psi=Psi,SumP=SumP)
                  
      if ((k==3)|(k==4)){ history <- data.frame(L=suiteL,code=code,Psi=Psi,SumP=SumP,C1=C1,C2=C2,C3=C3,C4=C4)}
                  
      pi <- resL$pi
      names(pi) <- paste("j=",resL$Spi,sep="")
      pi <- pi[order(resL$Spi)]
      Spi <- sort(resL$Spi)
      
      if (L==tail(valL,1)) {  
      if (k==2) {res <- list(cvge=1,trace=trace,history=history,
                            Spi=Spi,pi=pi,p=resL$p,
                            Psi=resL$Psi,SumP=resL$SumP)
      }
      if ((k==3)|(k==4)) {res <- list(cvge=1,trace=trace,history=history,
          Spi=Spi,pi=pi,p=resL$p, Psi=resL$Psi,SumP=resL$SumP,C1=Cvge$C1,C2=Cvge$C2,C3=Cvge$C3,C4=Cvge$C4) 
      }
                  
      return(res)
      }          
      
    } #### end : if (resL$code ==0) {
  } #### end : for (L in valL) {
}



