rBetaUniformFixedGamma<-function(n,alpha,beta,gamma){
  m0=round(n*gamma,digits = 0)
  m1=n-m0
  return(c(runif(m0), rbeta(m1,alpha,beta),m0))
}

rBetaUniform<-function(n,alpha,beta,gamma){
  m0=sum(rbinom(n,1,gamma))
  m1=n-m0
  return(c(runif(m0), rbeta(m1,alpha,beta),m0))
}

dBetaUniform<-function(x,alpha,beta,gamma){
  return(gamma+(1-gamma)*dbeta(x, alpha, beta))
}

pBetaUniform<-function(x,alpha,beta,gamma){
  return(gamma*x+(1-gamma)*pbeta(x, alpha, beta))
}
############################ MLEs ############################################
EMequations<-function(alphaBeta,xData,gammaT,cT){
  eqAlpha=sum((log(xData)-digamma(alphaBeta[1])+digamma(alphaBeta[1]+alphaBeta[2]))*cT/(gammaT+cT))
  eqBeta=sum((log(1-xData)-digamma(alphaBeta[2])+digamma(alphaBeta[1]+alphaBeta[2]))*cT/(gammaT+cT))
  return(c(eqAlpha,eqBeta))
}

updateEstimates<-function(xData,alphaT,betaT,gammaT,equationError,maxEquationIterations){
  cT=(1-gammaT)*(xData^(alphaT-1))*((1-xData)^(betaT-1))/beta(alphaT,betaT)
  gammaT1=mean(gammaT/( gammaT + (1-gammaT)*(xData^(alphaT-1))*((1-xData)^(betaT-1))/beta(alphaT,betaT)))
  EMeqSolutions=try(nleqslv(c(alphaT,betaT), EMequations, xData=xData,gammaT=gammaT,cT=cT,
                            method = "Newton", control=list(allowSingular=TRUE,xtol=equationError,maxit=maxEquationIterations)),silent = T)
  if("try-error" %in% class(EMeqSolutions)){
    alphaT1=-1
    betaT1=-1
  }else{
    alphaT1=EMeqSolutions$x[1]
    betaT1=EMeqSolutions$x[2]
  }
  return(c(alphaT1,betaT1,gammaT1))
}
############################################################################

############################ MPS estimation ############################################
MPSequations<-function(alphaBetaMPS,xData,c1T,c2T,gammaMPS){
  betaDer=pbetaDiff(x = xData,  p = alphaBetaMPS[1], q = alphaBetaMPS[2])
  betaDer1A=c(betaDer$dp,0)
  betaDer0A=c(0,betaDer$dp)
  betaDer1B=c(betaDer$dq,0)
  betaDer0B=c(0,betaDer$dq)
  xData1=c(xData,1)
  xData0=c(0,xData)
  betaDiff=pbeta(xData1,alphaBetaMPS[1],alphaBetaMPS[2])-pbeta(xData0,alphaBetaMPS[1],alphaBetaMPS[2])
  eqAlpha=sum(((betaDer1A-betaDer0A)/betaDiff)*c2T/(c1T+c2T))
  eqBeta=sum(((betaDer1B-betaDer0B)/betaDiff)*c2T/(c1T+c2T))
  return(c(eqAlpha,eqBeta))
}

MPSupdateEstimates<-function(xData,alphaTmps,betaTmps,gammaTmps,equationError,maxEquationIterations){
  xData1=c(xData,1)
  xData0=c(0,xData)
  c1T=gammaTmps*(xData1-xData0)
  c2T=(1-gammaTmps)*(pbeta(xData1,alphaTmps,betaTmps)-pbeta(xData0,alphaTmps,betaTmps))
  gammaT1mps=mean(c1T/(c1T+c2T))

  MPSeqSolutions=try(nleqslv(c(alphaTmps,betaTmps), MPSequations, xData=xData,c1T=c1T,c2T=c2T,
                             method = "Newton", control=list(allowSingular=TRUE,xtol=equationError,maxit=maxEquationIterations)),silent = T)
  if("try-error" %in% class(MPSeqSolutions)){
    alphaT1mps=-1
    betaT1mps=-1
  }else{
  alphaT1mps=MPSeqSolutions$x[1]
  betaT1mps=MPSeqSolutions$x[2]
  }
  return(c(alphaT1mps,betaT1mps,gammaT1mps))
}
############################################################################

############################ Method of moments ############################################
MMequations<-function(parMM,xData){
  mu1=mean(xData)
  mu2=mean(xData^2)
  mu3=mean(xData^3)
  eq1=-mu1+parMM[3]/2+(1-parMM[3])*parMM[1]/(parMM[1]+parMM[2])
  eq2=-mu2+parMM[3]/3+(1-parMM[3])*parMM[1]*(parMM[1]+1)/((parMM[1]+parMM[2])*(parMM[1]+parMM[2]+1))
  eq3=-mu3+parMM[3]/4+(1-parMM[3])*parMM[1]*(parMM[1]+1)*(parMM[1]+2)/((parMM[1]+parMM[2])*(parMM[1]+parMM[2]+1)*(parMM[1]+parMM[2]+2))
  return(c(eq1,eq2,eq3))
}
############################################################################

############################################################################
############################################################################
############################ Bootstrap methods #############################
############################################################################
############################################################################

############################ MLE ############################################
bootstrapMLE<-function(alpha,beta,gamma,n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations){
  alphaT=rep(alpha,trialsB)
  betaT=rep(beta,trialsB)
  gammaT=rep(gamma,trialsB)
  iterations=rep(0,trialsB)
  i=1
  while (i<=trialsB) {
    simBU=rBetaUniform(n,alpha,beta,gamma)
    xData=sort(simBU[1:n])
    j=0
    estError=1+estPrecision
    alphaT[i]=alpha
    betaT[i]=beta
    gammaT[i]=gamma
    while (estError>estPrecision && j<maxIterations) {
      updatedEst=updateEstimates(xData,alphaT[i],betaT[i],gammaT[i],equationError,maxEquationIterations)
      estError= abs(alphaT[i]-updatedEst[1])+abs(betaT[i]-updatedEst[2])+abs(gammaT[i]-updatedEst[3])
      alphaT[i]=updatedEst[1]
      betaT[i]=updatedEst[2]
      gammaT[i]=updatedEst[3]
      j=j+1
      if(updatedEst[1]<0 || updatedEst[2]<0 || updatedEst[3]<0 ){break}
    }
    iterations[i]=j
    i=i+1
  }
  return(list("alpha"=alphaT,"beta"=betaT,"gamma"=gammaT,"iterations"=iterations))
}

############################ MPS ############################################
bootstrapMPS<-function(alpha,beta,gamma,n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations){
  alphaTmps=rep(alpha,trialsB)
  betaTmps=rep(beta,trialsB)
  gammaTmps=rep(gamma,trialsB)
  iterationsMPS=rep(0,trialsB)
  i=1
  while (i<=trialsB) {
    simBU=rBetaUniform(n,alpha,beta,gamma)
    xData=sort(simBU[1:n])
    j=0
    estError=1+estPrecision
    alphaTmps[i]=alpha
    betaTmps[i]=beta
    gammaTmps[i]=gamma
    while (estError>estPrecision && j<maxIterations) {
      updatedEst=MPSupdateEstimates(xData,alphaTmps[i],betaTmps[i],gammaTmps[i],equationError,maxEquationIterations)
      estError= abs(alphaTmps[i]-updatedEst[1])+abs(betaTmps[i]-updatedEst[2])+abs(gammaTmps[i]-updatedEst[3])
      alphaTmps[i]=updatedEst[1]
      betaTmps[i]=updatedEst[2]
      gammaTmps[i]=updatedEst[3]
      j=j+1
      if(updatedEst[1]<0 || updatedEst[2]<0 || updatedEst[3]<0 ){break}
    }
  iterationsMPS[i]=j
  i=i+1
  }
  return(list("alpha"=alphaTmps,"beta"=betaTmps,"gamma"=gammaTmps,"iterations"=iterationsMPS))
}

############################ MM ############################################
bootstrapMM<-function(alpha,beta,gamma,n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations){
  alphaMM=rep(0,trialsB)
  betaMM=rep(0,trialsB)
  gammaMM=rep(0,trialsB)
  codeMM=rep(0,trialsB)
  i=1
  while (i<=trialsB) {
    simBU=rBetaUniform(n,alpha,beta,gamma)
    xData=sort(simBU[1:n])
    MMeqSolutions=try(nleqslv(c(alpha,beta,gamma), MMequations, xData=xData,
                              method = "Newton", control=list(allowSingular=TRUE,xtol=equationError,maxit=maxEquationIterations)),silent = T)
    if("try-error" %in% class(MMeqSolutions)){break
      }else{
      alphaMM[i]=MMeqSolutions$x[1]
      betaMM[i]=MMeqSolutions$x[2]
      gammaMM[i]=MMeqSolutions$x[3]
      codeMM[i]=MMeqSolutions$termcd
    }
    i=i+1
  }
  return(list("alpha"=alphaMM,"beta"=betaMM,"gamma"=gammaMM,"code"=codeMM))
}

