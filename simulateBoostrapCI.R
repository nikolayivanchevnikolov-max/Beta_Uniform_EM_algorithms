rm(list=ls())
library(nleqslv)
library(mnorm)
###############################################
############ Input parameters #################
###############################################
args=c(2024,100,0.8,0.3,1.2,1020)

seedNum=as.numeric(args[1])
n=as.numeric(args[2]) # Try also n=500, n=1000, n=5000
gamma=as.numeric(args[3]) # Try also gamma=0.70, gamma=0.90, gamma=0.95
alpha=as.numeric(args[4]) # Try also alpha=0.10, alpha=0.50
beta=as.numeric(args[5]) # Try also beta=1.5, beta=2.00

trials=as.numeric(args[6])
trialsB=100
set.seed(seedNum)
###############################################
alphaSig=0.05

estPrecision=0.0001
maxIterations=1000
equationError=10^{-6}
maxEquationIterations=1000

alpha0=alpha
beta0=beta
gamma0=gamma

####################################################################################
############################### FUNCTIONS ##########################################
####################################################################################

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


####################################################################################
####################################################################################
####################################################################################

alphaT=rep(alpha0,trials)
betaT=rep(beta0,trials)
gammaT=rep(gamma0,trials)
alphaT_BCA=matrix(0,trials,2)
betaT_BCA=matrix(0,trials,2)
gammaT_BCA=matrix(0,trials,2)
alphaT_CI=matrix(0,trials,2)
betaT_CI=matrix(0,trials,2)
gammaT_CI=matrix(0,trials,2)

alphaTmps=rep(alpha0,trials)
betaTmps=rep(beta0,trials)
gammaTmps=rep(gamma0,trials)
alphaTmps_BCA=matrix(0,trials,2)
betaTmps_BCA=matrix(0,trials,2)
gammaTmps_BCA=matrix(0,trials,2)
alphaTmps_CI=matrix(0,trials,2)
betaTmps_CI=matrix(0,trials,2)
gammaTmps_CI=matrix(0,trials,2)

alphaMM=rep(0,trials)
betaMM=rep(0,trials)
gammaMM=rep(0,trials)
alphaMM_BCA=matrix(0,trials,2)
betaMM_BCA=matrix(0,trials,2)
gammaMM_BCA=matrix(0,trials,2)
alphaMM_CI=matrix(0,trials,2)
betaMM_CI=matrix(0,trials,2)
gammaMM_CI=matrix(0,trials,2)

codeMM=rep(0,trials)
iterations=rep(0,trials)
iterationsMPS=rep(0,trials)
compTime=rep(0,trials)
startTime=Sys.time()

za2=qnorm(1-alphaSig/2)
i=1
countWhile=0
while (i<=trials) {
  countWhile=countWhile+1
  simBU=rBetaUniformFixedGamma(n,alpha,beta,gamma) ## Alternative: rBetaUniform(n,alpha,beta,gamma)
  xData=sort(simBU[1:n])
  
  ################### MLE #############################################
  j=0
  estError=1+estPrecision
  alphaT[i]=alpha0
  betaT[i]=beta0
  gammaT[i]=gamma0
  while (estError>estPrecision && j<maxIterations) {
    updatedEst=updateEstimates(xData,alphaT[i],betaT[i],gammaT[i],equationError,maxEquationIterations)
    estError= abs(alphaT[i]-updatedEst[1])+abs(betaT[i]-updatedEst[2])+abs(gammaT[i]-updatedEst[3])
    alphaT[i]=updatedEst[1]
    betaT[i]=updatedEst[2]
    gammaT[i]=updatedEst[3]
    if(updatedEst[1]<0 || updatedEst[2]<0 || updatedEst[3]<0 ){break}
    j=j+1
  }
  iterations[i]=j
  #####################################################################
  
  ################### MPS #############################################
  j=0
  estError=1+estPrecision
  alphaTmps[i]=alpha0
  betaTmps[i]=beta0
  gammaTmps[i]=gamma0
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
  #####################################################################
  
  #################### MM #############################################
  MMeqSolutions=try(nleqslv(c(alpha0,beta0,gamma0), MMequations, xData=xData,
                            method = "Newton", control=list(allowSingular=TRUE,xtol=equationError,maxit=maxEquationIterations)),silent = T)
  if("try-error" %in% class(MMeqSolutions)){next
  }else{
    alphaMM[i]=MMeqSolutions$x[1]
    betaMM[i]=MMeqSolutions$x[2]
    gammaMM[i]=MMeqSolutions$x[3]
    codeMM[i]=MMeqSolutions$termcd
  }
  #####################################################################
  if (!(is.na(alphaT[i]) || is.na(betaT[i]) || is.na(gammaT[i]) ||
        is.na(alphaTmps[i]) || is.na(betaTmps[i]) || is.na(gammaTmps[i]) ||
        is.na(alphaMM[i]) || is.na(betaMM[i]) || is.na(gammaMM[i]) ||
        is.nan(alphaT[i]) || is.nan(betaT[i]) || is.nan(gammaT[i]) ||
        is.nan(alphaTmps[i]) || is.nan(betaTmps[i]) || is.nan(gammaTmps[i]) ||
        is.nan(alphaMM[i]) || is.nan(betaMM[i]) || is.nan(gammaMM[i]) ||
        is.infinite(alphaT[i]) || is.infinite(betaT[i]) || is.infinite(gammaT[i]) ||
        is.infinite(alphaTmps[i]) || is.infinite(betaTmps[i]) || is.infinite(gammaTmps[i]) ||
        is.infinite(alphaMM[i]) || is.infinite(betaMM[i]) || is.infinite(gammaMM[i]) ||
        alphaT[i]<0 || betaT[i]<0 || gammaT[i]<0 ||
        alphaTmps[i]<0 || betaTmps[i]<0 || gammaTmps[i]<0 ||
        alphaMM[i]<0 || betaMM[i]<0 || gammaMM[i]<0 ||
        alphaT[i]>alpha0*50 || betaT[i]>beta0*50 || gammaT[i]>1 ||
        alphaTmps[i]>alpha0*50 || betaTmps[i]>beta0*50 || gammaTmps[i]>1 ||
        alphaMM[i]>alpha0*50 || betaMM[i]>beta0*50 || gammaMM[i]>1)){
    #################################################################################    
    # Bootstrap MLE
    bsMLE=bootstrapMLE(alphaT[i],betaT[i],gammaT[i],n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations)
    realBSindex=which(!is.na(bsMLE$alpha) & !is.na(bsMLE$beta) & !is.na(bsMLE$gamma) &
                        !is.nan(bsMLE$alpha) & !is.nan(bsMLE$beta) & !is.nan(bsMLE$gamma) &
                        !is.infinite(bsMLE$alpha) & !is.infinite(bsMLE$beta) & !is.infinite(bsMLE$gamma) &
                        bsMLE$alpha>0 & bsMLE$beta>0 & bsMLE$gamma>0 &
                        bsMLE$alpha<alpha0*50 & bsMLE$beta<beta0*50 & bsMLE$gamma<1)
    trialsBB=length(realBSindex)
    bsMLE$alpha=bsMLE$alpha[realBSindex]
    bsMLE$beta=bsMLE$beta[realBSindex]
    bsMLE$gamma=bsMLE$gamma[realBSindex]
    
    z0MLE=rep(0,3)
    z0MLE[1]=qnorm(sum(bsMLE$alpha<alphaT[i])/trialsBB)
    z0MLE[2]=qnorm(sum(bsMLE$beta<betaT[i])/trialsBB)
    z0MLE[3]=qnorm(sum(bsMLE$gamma<gammaT[i])/trialsBB)
    a1MLE=pnorm(2*z0MLE-za2)
    a2MLE=pnorm(2*z0MLE+za2)
    
    alphaTjk=rep(0,n)
    betaTjk=rep(0,n)
    gammaTjk=rep(0,n)
    for (k in 1:n) {
      j=0
      estError=1+estPrecision
      alphaTjk[k]=alphaT[i]
      betaTjk[k]=betaT[i]
      gammaTjk[k]=gammaT[i]
      while (estError>estPrecision && j<maxIterations) {
        updatedEst=updateEstimates(xData[-k],alphaTjk[k],betaTjk[k],gammaTjk[k],equationError,maxEquationIterations)
        estError= abs(alphaTjk[k]-updatedEst[1])+abs(betaTjk[k]-updatedEst[2])+abs(gammaTjk[k]-updatedEst[3])
        alphaTjk[k]=updatedEst[1]
        betaTjk[k]=updatedEst[2]
        gammaTjk[k]=updatedEst[3]
        if(updatedEst[1]<0 || updatedEst[2]<0 || updatedEst[3]<0 ){break}
        j=j+1
      }
    }
    aiMLE=rep(0,3)
    aiMLE[1]=((sum((mean(alphaTjk)-alphaTjk)^2))^(-3/2))*sum((mean(alphaTjk)-alphaTjk)^3)/6
    aiMLE[2]=((sum((mean(betaTjk)-betaTjk)^2))^(-3/2))*sum((mean(betaTjk)-betaTjk)^3)/6
    aiMLE[3]=((sum((mean(gammaTjk)-gammaTjk)^2))^(-3/2))*sum((mean(gammaTjk)-gammaTjk)^3)/6
    
    if (is.na(aiMLE[1]) || is.na(aiMLE[2]) || is.na(aiMLE[3]) ||
        is.nan(aiMLE[1]) || is.nan(aiMLE[2]) || is.nan(aiMLE[3]) ||
        is.infinite(aiMLE[1]) || is.infinite(aiMLE[2]) || is.infinite(aiMLE[3])){
      aiMLE=0
    }
    ai1MLE=pnorm(z0MLE+(z0MLE-za2)/(1-aiMLE*(z0MLE-za2)))
    ai2MLE=pnorm(z0MLE+(z0MLE+za2)/(1-aiMLE*(z0MLE+za2)))
    
    alphaT_BCA[i,]=quantile(bsMLE$alpha,probs = c(ai1MLE[1],ai2MLE[1]))
    betaT_BCA[i,]=quantile(bsMLE$beta,probs = c(ai1MLE[2],ai2MLE[2]))
    gammaT_BCA[i,]=quantile(bsMLE$gamma,probs = c(ai1MLE[3],ai2MLE[3]))
    
    alphaT_CI[i,]=quantile(bsMLE$alpha,probs = c(alphaSig/2,1-alphaSig/2))
    betaT_CI[i,]=quantile(bsMLE$beta,probs = c(alphaSig/2,1-alphaSig/2))
    gammaT_CI[i,]=quantile(bsMLE$gamma,probs = c(alphaSig/2,1-alphaSig/2))
    
    #################################################################################    
    # Bootstrap MPS
    bsMPS=bootstrapMPS(alphaTmps[i],betaTmps[i],gammaTmps[i],n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations)
    realBSindex=which(!is.na(bsMPS$alpha) & !is.na(bsMPS$beta) & !is.na(bsMPS$gamma) &
                        !is.nan(bsMPS$alpha) & !is.nan(bsMPS$beta) & !is.nan(bsMPS$gamma) &
                        !is.infinite(bsMPS$alpha) & !is.infinite(bsMPS$beta) & !is.infinite(bsMPS$gamma) &
                        bsMPS$alpha>0 & bsMPS$beta>0 & bsMPS$gamma>0 &
                        bsMPS$alpha<alpha0*50 & bsMPS$beta<beta0*50 & bsMPS$gamma<1)
    trialsBB=length(realBSindex)
    bsMPS$alpha=bsMPS$alpha[realBSindex]
    bsMPS$beta=bsMPS$beta[realBSindex]
    bsMPS$gamma=bsMPS$gamma[realBSindex]
    
    z0MPS=rep(0,3)
    z0MPS[1]=qnorm(sum(bsMPS$alpha<alphaTmps[i])/trialsBB)
    z0MPS[2]=qnorm(sum(bsMPS$beta<betaTmps[i])/trialsBB)
    z0MPS[3]=qnorm(sum(bsMPS$gamma<gammaTmps[i])/trialsBB)
    a1MPS=pnorm(2*z0MPS-za2)
    a2MPS=pnorm(2*z0MPS+za2)
    
    alphaTmpsjk=rep(0,n)
    betaTmpsjk=rep(0,n)
    gammaTmpsjk=rep(0,n)
    for (k in 1:n) {
      j=0
      estError=1+estPrecision
      alphaTmpsjk[k]=alphaTmps[i]
      betaTmpsjk[k]=betaTmps[i]
      gammaTmpsjk[k]=gammaTmps[i]
      while (estError>estPrecision && j<maxIterations) {
        updatedEst=MPSupdateEstimates(xData[-k],alphaTmpsjk[k],betaTmpsjk[k],gammaTmpsjk[k],equationError,maxEquationIterations)
        estError= abs(alphaTmpsjk[k]-updatedEst[1])+abs(betaTmpsjk[k]-updatedEst[2])+abs(gammaTmpsjk[k]-updatedEst[3])
        alphaTmpsjk[k]=updatedEst[1]
        betaTmpsjk[k]=updatedEst[2]
        gammaTmpsjk[k]=updatedEst[3]
        j=j+1
        if(updatedEst[1]<0 || updatedEst[2]<0 || updatedEst[3]<0 ){break}
      }
    }
    aiMPS=rep(0,3)
    aiMPS[1]=((sum((mean(alphaTmpsjk)-alphaTmpsjk)^2))^(-3/2))*sum((mean(alphaTmpsjk)-alphaTmpsjk)^3)/6
    aiMPS[2]=((sum((mean(betaTmpsjk)-betaTmpsjk)^2))^(-3/2))*sum((mean(betaTmpsjk)-betaTmpsjk)^3)/6
    aiMPS[3]=((sum((mean(gammaTmpsjk)-gammaTmpsjk)^2))^(-3/2))*sum((mean(gammaTmpsjk)-gammaTmpsjk)^3)/6
    
    if (is.na(aiMPS[1]) || is.na(aiMPS[2]) || is.na(aiMPS[3]) ||
        is.nan(aiMPS[1]) || is.nan(aiMPS[2]) || is.nan(aiMPS[3]) ||
        is.infinite(aiMPS[1]) || is.infinite(aiMPS[2]) || is.infinite(aiMPS[3])){
      aiMPS=0
    }
    ai1MPS=pnorm(z0MPS+(z0MPS-za2)/(1-aiMLE*(z0MPS-za2)))
    ai2MPS=pnorm(z0MPS+(z0MPS+za2)/(1-aiMLE*(z0MPS+za2)))
    alphaTmps_BCA[i,]=quantile(bsMPS$alpha,probs = c(ai1MPS[1],ai2MPS[1]))
    betaTmps_BCA[i,]=quantile(bsMPS$beta,probs = c(ai1MPS[2],ai2MPS[2]))
    gammaTmps_BCA[i,]=quantile(bsMPS$gamma,probs = c(ai1MPS[3],ai2MPS[3]))
    
    alphaTmps_CI[i,]=quantile(bsMPS$alpha,probs = c(alphaSig/2,1-alphaSig/2))
    betaTmps_CI[i,]=quantile(bsMPS$beta,probs = c(alphaSig/2,1-alphaSig/2))
    gammaTmps_CI[i,]=quantile(bsMPS$gamma,probs = c(alphaSig/2,1-alphaSig/2))
    
    
    #################################################################################    
    # Bootstrap MM
    bsMM=bootstrapMM(alphaMM[i],betaMM[i],gammaMM[i],n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations)
    realBSindex=which(!is.na(bsMM$alpha) & !is.na(bsMM$beta) & !is.na(bsMM$gamma) &
                        !is.nan(bsMM$alpha) & !is.nan(bsMM$beta) & !is.nan(bsMM$gamma) &
                        !is.infinite(bsMM$alpha) & !is.infinite(bsMM$beta) & !is.infinite(bsMM$gamma) &
                        bsMM$alpha>0 & bsMM$beta>0 & bsMM$gamma>0 &
                        bsMM$alpha<alpha0*50 & bsMM$beta<beta0*50 & bsMM$gamma<1)
    trialsBB=length(realBSindex)
    bsMM$alpha=bsMM$alpha[realBSindex]
    bsMM$beta=bsMM$beta[realBSindex]
    bsMM$gamma=bsMM$gamma[realBSindex]
    
    z0MM=rep(0,3)
    z0MM[1]=qnorm(sum(bsMM$alpha<alphaMM[i])/trialsBB)
    z0MM[2]=qnorm(sum(bsMM$beta<betaMM[i])/trialsBB)
    z0MM[3]=qnorm(sum(bsMM$gamma<gammaMM[i])/trialsBB)
    a1MM=pnorm(2*z0MM-za2)
    a2MM=pnorm(2*z0MM+za2)
    
    alphaMMjk=rep(0,n)
    betaMMjk=rep(0,n)
    gammaMMjk=rep(0,n)
    for (k in 1:n) {
      MMeqSolutions=try(nleqslv(c(alphaMM[i],betaMM[i],gammaMM[i]), MMequations, xData=xData[-k],
                                method = "Newton", control=list(allowSingular=TRUE,xtol=equationError,maxit=maxEquationIterations)),silent = T)
      if("try-error" %in% class(MMeqSolutions)){next
      }else{
        alphaMMjk[k]=MMeqSolutions$x[1]
        betaMMjk[k]=MMeqSolutions$x[2]
        gammaMMjk[k]=MMeqSolutions$x[3]
      }
      if(MMeqSolutions$x[1]<0 || MMeqSolutions$x[2]<0 || MMeqSolutions$x[3]<0 ){next}
    }
    aiMM=rep(0,3)
    aiMM[1]=((sum((mean(alphaMMjk)-alphaMMjk)^2))^(-3/2))*sum((mean(alphaMMjk)-alphaMMjk)^3)/6
    aiMM[2]=((sum((mean(betaMMjk)-betaMMjk)^2))^(-3/2))*sum((mean(betaMMjk)-betaMMjk)^3)/6
    aiMM[3]=((sum((mean(gammaMMjk)-gammaMMjk)^2))^(-3/2))*sum((mean(gammaMMjk)-gammaMMjk)^3)/6
    
    if (is.na(aiMM[1]) || is.na(aiMM[2]) || is.na(aiMM[3]) ||
        is.nan(aiMM[1]) || is.nan(aiMM[2]) || is.nan(aiMM[3]) ||
        is.infinite(aiMM[1]) || is.infinite(aiMM[2]) || is.infinite(aiMM[3])){
      aiMM=0
    }
    ai1MM=pnorm(z0MM+(z0MM-za2)/(1-aiMM*(z0MM-za2)))
    ai2MM=pnorm(z0MM+(z0MM+za2)/(1-aiMM*(z0MM+za2)))
    alphaMM_BCA[i,]=quantile(bsMM$alpha,probs = c(ai1MM[1],ai2MM[1]))
    betaMM_BCA[i,]=quantile(bsMM$beta,probs = c(ai1MM[2],ai2MM[2]))
    gammaMM_BCA[i,]=quantile(bsMM$gamma,probs = c(ai1MM[3],ai2MM[3]))
    
    alphaMM_CI[i,]=quantile(bsMM$alpha,probs = c(alphaSig/2,1-alphaSig/2))
    betaMM_CI[i,]=quantile(bsMM$beta,probs = c(alphaSig/2,1-alphaSig/2))
    gammaMM_CI[i,]=quantile(bsMM$gamma,probs = c(alphaSig/2,1-alphaSig/2))
    
    
    endTime=Sys.time()
    compTime[i]=difftime(endTime,startTime,units = 'mins')
    startTime=Sys.time()
    i=i+1
  }
}

##################################################################
####################### RESULTS ##################################
##################################################################
bcaIndex=which(!is.nan(alphaT_BCA[,1]) & !is.nan(alphaT_BCA[,2]) &
                 !is.nan(betaT_BCA[,1]) & !is.nan(betaT_BCA[,2]) &
                 !is.nan(gammaT_BCA[,1]) & !is.nan(gammaT_BCA[,2]))

ciIndex=which(!is.nan(alphaT_CI[,1]) & !is.nan(alphaT_CI[,2]) &
                !is.nan(betaT_CI[,1]) & !is.nan(betaT_CI[,2]) &
                !is.nan(gammaT_CI[,1]) & !is.nan(gammaT_CI[,2]))
estIndex=1:trials

###################### MLE ######################################
resMLE=matrix(0,6,3)

resMLE[4,1]=mean(alphaT_BCA[bcaIndex,2]-alphaT_BCA[bcaIndex,1])
resMLE[4,2]=mean(betaT_BCA[bcaIndex,2]-betaT_BCA[bcaIndex,1])
resMLE[4,3]=mean(gammaT_BCA[bcaIndex,2]-gammaT_BCA[bcaIndex,1])

resMLE[3,1]=mean((alphaT_BCA[bcaIndex,2]>=alpha)*(alphaT_BCA[bcaIndex,1]<=alpha))
resMLE[3,2]=mean((betaT_BCA[bcaIndex,2]>=beta)*(betaT_BCA[bcaIndex,1]<=beta))
resMLE[3,3]=mean((gammaT_BCA[bcaIndex,2]>=gamma)*(gammaT_BCA[bcaIndex,1]<=gamma))

resMLE[6,1]=mean(alphaT_CI[ciIndex,2]-alphaT_CI[ciIndex,1])
resMLE[6,2]=mean(betaT_CI[ciIndex,2]-betaT_CI[ciIndex,1])
resMLE[6,3]=mean(gammaT_CI[ciIndex,2]-gammaT_CI[ciIndex,1])

resMLE[5,1]=mean((alphaT_CI[ciIndex,2]>=alpha)*(alphaT_CI[ciIndex,1]<=alpha))
resMLE[5,2]=mean((betaT_CI[ciIndex,2]>=beta)*(betaT_CI[ciIndex,1]<=beta))
resMLE[5,3]=mean((gammaT_CI[ciIndex,2]>=gamma)*(gammaT_CI[ciIndex,1]<=gamma))

resMLE[1,1]=mean((alphaT[estIndex]-alpha)^2)
resMLE[1,2]=mean((betaT[estIndex]-beta)^2)
resMLE[1,3]=mean((gammaT[estIndex]-gamma)^2)

resMLE[2,1]=(mean(alphaT[estIndex])-alpha)
resMLE[2,2]=(mean(betaT[estIndex])-beta)
resMLE[2,3]=(mean(gammaT[estIndex])-gamma)


###################### MPS ######################################
bcaIndex=which(!is.nan(alphaTmps_BCA[,1]) & !is.nan(alphaTmps_BCA[,2]) &
                 !is.nan(betaTmps_BCA[,1]) & !is.nan(betaTmps_BCA[,2]) &
                 !is.nan(gammaTmps_BCA[,1]) & !is.nan(gammaTmps_BCA[,2]))

ciIndex=which(!is.nan(alphaTmps_CI[,1]) & !is.nan(alphaTmps_CI[,2]) &
                !is.nan(betaTmps_CI[,1]) & !is.nan(betaTmps_CI[,2]) &
                !is.nan(gammaTmps_CI[,1]) & !is.nan(gammaTmps_CI[,2]))
estIndex=1:trials

resMPS=matrix(0,6,3)

resMPS[4,1]=mean(alphaTmps_BCA[bcaIndex,2]-alphaTmps_BCA[bcaIndex,1])
resMPS[4,2]=mean(betaTmps_BCA[bcaIndex,2]-betaTmps_BCA[bcaIndex,1])
resMPS[4,3]=mean(gammaTmps_BCA[bcaIndex,2]-gammaTmps_BCA[bcaIndex,1])

resMPS[3,1]=mean((alphaTmps_BCA[bcaIndex,2]>=alpha)*(alphaTmps_BCA[bcaIndex,1]<=alpha))
resMPS[3,2]=mean((betaTmps_BCA[bcaIndex,2]>=beta)*(betaTmps_BCA[bcaIndex,1]<=beta))
resMPS[3,3]=mean((gammaTmps_BCA[bcaIndex,2]>=gamma)*(gammaTmps_BCA[bcaIndex,1]<=gamma))

resMPS[6,1]=mean(alphaTmps_CI[ciIndex,2]-alphaTmps_CI[ciIndex,1])
resMPS[6,2]=mean(betaTmps_CI[ciIndex,2]-betaTmps_CI[ciIndex,1])
resMPS[6,3]=mean(gammaTmps_CI[ciIndex,2]-gammaTmps_CI[ciIndex,1])

resMPS[5,1]=mean((alphaTmps_CI[ciIndex,2]>=alpha)*(alphaTmps_CI[ciIndex,1]<=alpha))
resMPS[5,2]=mean((betaTmps_CI[ciIndex,2]>=beta)*(betaTmps_CI[ciIndex,1]<=beta))
resMPS[5,3]=mean((gammaTmps_CI[ciIndex,2]>=gamma)*(gammaTmps_CI[ciIndex,1]<=gamma))

resMPS[1,1]=mean((alphaTmps[estIndex]-alpha)^2)
resMPS[1,2]=mean((betaTmps[estIndex]-beta)^2)
resMPS[1,3]=mean((gammaTmps[estIndex]-gamma)^2)

resMPS[2,1]=(mean(alphaTmps[estIndex])-alpha)
resMPS[2,2]=(mean(betaTmps[estIndex])-beta)
resMPS[2,3]=(mean(gammaTmps[estIndex])-gamma)

###################### MM ######################################
bcaIndex=which(!is.nan(alphaMM_BCA[,1]) & !is.nan(alphaMM_BCA[,2]) &
                 !is.nan(betaMM_BCA[,1]) & !is.nan(betaMM_BCA[,2]) &
                 !is.nan(gammaMM_BCA[,1]) & !is.nan(gammaMM_BCA[,2]))

ciIndex=which(!is.nan(alphaMM_CI[,1]) & !is.nan(alphaMM_CI[,2]) &
                !is.nan(betaMM_CI[,1]) & !is.nan(betaMM_CI[,2]) &
                !is.nan(gammaMM_CI[,1]) & !is.nan(gammaMM_CI[,2]))
estIndex=1:trials

resMM=matrix(0,6,3)

resMM[4,1]=mean(alphaMM_BCA[bcaIndex,2]-alphaMM_BCA[bcaIndex,1])
resMM[4,2]=mean(betaMM_BCA[bcaIndex,2]-betaMM_BCA[bcaIndex,1])
resMM[4,3]=mean(gammaMM_BCA[bcaIndex,2]-gammaMM_BCA[bcaIndex,1])

resMM[3,1]=mean((alphaMM_BCA[bcaIndex,2]>=alpha)*(alphaMM_BCA[bcaIndex,1]<=alpha))
resMM[3,2]=mean((betaMM_BCA[bcaIndex,2]>=beta)*(betaMM_BCA[bcaIndex,1]<=beta))
resMM[3,3]=mean((gammaMM_BCA[bcaIndex,2]>=gamma)*(gammaMM_BCA[bcaIndex,1]<=gamma))

resMM[6,1]=mean(alphaMM_CI[ciIndex,2]-alphaMM_CI[ciIndex,1])
resMM[6,2]=mean(betaMM_CI[ciIndex,2]-betaMM_CI[ciIndex,1])
resMM[6,3]=mean(gammaMM_CI[ciIndex,2]-gammaMM_CI[ciIndex,1])

resMM[5,1]=mean((alphaMM_CI[ciIndex,2]>=alpha)*(alphaMM_CI[ciIndex,1]<=alpha))
resMM[5,2]=mean((betaMM_CI[ciIndex,2]>=beta)*(betaMM_CI[ciIndex,1]<=beta))
resMM[5,3]=mean((gammaMM_CI[ciIndex,2]>=gamma)*(gammaMM_CI[ciIndex,1]<=gamma))

resMM[1,1]=mean((alphaMM[estIndex]-alpha)^2)
resMM[1,2]=mean((betaMM[estIndex]-beta)^2)
resMM[1,3]=mean((gammaMM[estIndex]-gamma)^2)

resMM[2,1]=(mean(alphaMM[estIndex])-alpha)
resMM[2,2]=(mean(betaMM[estIndex])-beta)
resMM[2,3]=(mean(gammaMM[estIndex])-gamma)

###################### SAVE RESULTS ######################################
write.csv(resMLE, paste("resMLE_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gamma,"_trials_",trials,"_trialsB_",trialsB,".csv",sep=""),row.names = FALSE)
write.csv(resMPS, paste("resMPS_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gamma,"_trials_",trials,"_trialsB_",trialsB,".csv",sep=""),row.names = FALSE)
write.csv(resMM, paste("resMM_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gamma,"_trials_",trials,"_trialsB_",trialsB,".csv",sep=""),row.names = FALSE)

###################### SAVE ALL CIs ######################################
resDF=data.frame(alphaT, alphaTmps, alphaMM, 
                 alphaT_BCA[,1], alphaT_BCA[,2], alphaTmps_BCA[,1], alphaTmps_BCA[,2], alphaMM_BCA[,1], alphaMM_BCA[,2],
                 alphaT_CI[,1], alphaT_CI[,2], alphaTmps_CI[,1], alphaTmps_CI[,2], alphaMM_CI[,1], alphaMM_CI[,2],
                 betaT, betaTmps, betaMM, 
                 betaT_BCA[,1], betaT_BCA[,2], betaTmps_BCA[,1], betaTmps_BCA[,2], betaMM_BCA[,1], betaMM_BCA[,2],
                 betaT_CI[,1], betaT_CI[,2], betaTmps_CI[,1], betaTmps_CI[,2], betaMM_CI[,1], betaMM_CI[,2],
                 gammaT, gammaTmps, gammaMM, 
                 gammaT_BCA[,1], gammaT_BCA[,2], gammaTmps_BCA[,1], gammaTmps_BCA[,2], gammaMM_BCA[,1], gammaMM_BCA[,2],
                 gammaT_CI[,1], gammaT_CI[,2], gammaTmps_CI[,1], gammaTmps_CI[,2], gammaMM_CI[,1], gammaMM_CI[,2],
                 compTime)
write.csv(resDF, paste("resAllRaw_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gamma,"_trials_",trials,"_trialsB_",trialsB,".csv",sep=""),row.names = FALSE)