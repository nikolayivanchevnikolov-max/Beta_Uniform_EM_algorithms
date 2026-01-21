rm(list=ls())
library(nleqslv)
library(mnorm)
source('simulateBetaUniform.R')

options(scipen = 999)

data <- read.table("data\\bottomly_unadjusted_p-values.txt",sep="\t", header=TRUE)
genes_data <- read.csv("data\\genes_seqid.csv")
genes_data2 <- read.csv("data\\genes_seqname.csv")
allData=merge(x=data,y=genes_data, by="gene", all.x=TRUE)
allData2=merge(x=allData,y=genes_data2, by="gene", all.x=TRUE)

allData2$seqid2=gsub("chr", "", allData2$seqname)
allData2$seqid2=gsub("_random", "", allData2$seqid2)

allData2$seqDiff=(allData2$seqid != allData2$seqid2 | (is.na(allData2$seqid) & !is.na(allData2$seqid2)) | (is.na(allData2$seqid2) & !is.na(allData2$seqid))) & !(is.na(allData2$seqid) & is.na(allData2$seqid2))

table(allData2$seqid,useNA = "ifany")
table(allData2$seqid2,useNA = "ifany")
dat=allData2
chrName="12"
dat=allData2[which(allData2$seqid2==chrName),]

alpha=0.05
hist(dat$edgeR.pval.unadjusted)
abline(h=round(alpha*length(dat$edgeR.pval.unadjusted)),col=2,lty=2)
hist(dat$limma.pval.unadjusted)
abline(h=round(alpha*length(dat$limma.pval.unadjusted)),col=2,lty=2)
hist(dat$DESeq2.pval.unadjusted)
abline(h=round(alpha*length(dat$DESeq2.pval.unadjusted)),col=2,lty=2)

xData=dat$DESeq2.pval.unadjusted 

max(table(xData))
xData=sort(xData)
table(xData[which(duplicated(xData))])

i=0
while (sum(duplicated(xData))>0) {
  xData[which(duplicated(xData))]=xData[which(duplicated(xData))]+which(duplicated(xData))*10^(-8)
  xData=sort(xData)
  i=i+1
}

pdf(file = paste("plots_Chromosome/", chrName, ".pdf", sep=""), width = 8, height = 6)
par(mar=c(5,5,3.5,2))
h=hist(xData, probability = TRUE, main= paste("Chromosome",chrName, sep=" "), xlab="p-values", 
       cex.axis=1.75, cex.lab=1.75, cex.main=1.75,col="gray95")
abline(h=1,col="red4",lty=2, lwd=2)

estPrecision=0.0001
maxIterations=10000
equationError=10^{-6}
maxEquationIterations=1000

startTime=Sys.time()

#################### MM #############################################
MMeqSolutions=nleqslv(c(0.5,0.8,0.05), MMequations, xData=xData,
                      method = "Newton", control=list(allowSingular=TRUE,xtol=10^(-6),maxit=1000))
alphaMM=MMeqSolutions$x[1]
betaMM=MMeqSolutions$x[2]
gammaMM=MMeqSolutions$x[3]
codeMM=MMeqSolutions$termcd
#####################################################################
curve(dBetaUniform(x,alphaMM,betaMM,gammaMM), from=0, to=0.99, col="darkblue", lwd=2, add=TRUE)
#####################################################################

alphaT=alphaMM
betaT=betaMM
gammaT=gammaMM
################## MLE #############################################
j=0
estError=1+estPrecision
while (estError>estPrecision && j<maxIterations) {
  updatedEst=updateEstimates(xData,alphaT,betaT,gammaT,equationError,maxEquationIterations)
  estError= abs(alphaT-updatedEst[1])+abs(betaT-updatedEst[2])+abs(gammaT-updatedEst[3])
  alphaT=updatedEst[1]
  betaT=updatedEst[2]
  gammaT=updatedEst[3]
  j=j+1
}
iterationsMLE=j
#####################################################################
curve(dBetaUniform(x,alphaT,betaT,gammaT), from=0, to=0.99, col="red2", lwd=2, add=TRUE)
#####################################################################

alphaTmps=alphaT
betaTmps=betaT
gammaTmps=gammaT

# alphaTmps=alphaMM
# betaTmps=betaMM
# gammaTmps=gammaMM

# alphaTmps=0.5
# betaTmps=1.5
# gammaTmps=0.8

################### MPS #############################################
j=0
estError=1+estPrecision
while (estError>estPrecision && j<maxIterations) {
  updatedEst=MPSupdateEstimates(xData,alphaTmps,betaTmps,gammaTmps,equationError,maxEquationIterations)
  estError= abs(alphaTmps-updatedEst[1])+abs(betaTmps-updatedEst[2])+abs(gammaTmps-updatedEst[3])
  alphaTmps=updatedEst[1]
  betaTmps=updatedEst[2]
  gammaTmps=updatedEst[3]
  j=j+1
}
iterationsMPS=j
#####################################################################
#####################################################################
curve(dBetaUniform(x,alphaTmps,betaTmps,gammaTmps), from=0, to=0.99, col="darkgreen", lwd=2, lty=2,add=TRUE)
#####################################################################

legend("topright", legend=c("ML", "MPS", "MM"),
       col=c("red2","darkgreen","darkblue"), lty=c(1,2,1), cex=1.4, lwd=2)

dev.off()

endTime=Sys.time()
totalTime=endTime-startTime
print(totalTime)

library(DescTools)
AndersonDarlingTest(xData,pBetaUniform,alpha=alphaT,beta=betaT,gamma=gammaT)
AndersonDarlingTest(xData,punif)
AndersonDarlingTest(xData,pBetaUniform,alpha=alphaTmps,beta=betaTmps,gamma=gammaTmps)
AndersonDarlingTest(xData,pBetaUniform,alpha=alphaMM,beta=betaMM,gamma=gammaMM)
ks.test(xData,pBetaUniform,alpha=alphaT,beta=betaT,gamma=gammaT)
ks.test(xData,pBetaUniform,alpha=alphaTmps,beta=betaTmps,gamma=gammaTmps)
ks.test(xData,pBetaUniform,alpha=alphaMM,beta=betaMM,gamma=gammaMM)

log(prod(dBetaUniform(xData,alpha=alphaT,beta=betaT,gamma=gammaT)))
log(prod(dBetaUniform(xData,alpha=alphaTmps,beta=betaTmps,gamma=gammaTmps)))
log(prod(dBetaUniform(xData,alpha=alphaMM,beta=betaMM,gamma=gammaMM)))

qqplot(xData,rBetaUniform(n=400,alpha=alphaMM,beta=betaMM,gamma=gammaMM)[1:400])
abline(0,1,col=2)

################### classification MLE #############################################
classTh=0.5#0.739 chr12 #0.845 chr13
BetaUnifClass=as.integer((1-gammaT)*dbeta(xData, alphaT, betaT)/(gammaT+(1-gammaT)*dbeta(xData, alphaT, betaT))>classTh)
cData=as.data.frame(cbind(xData,BetaUnifClass))
sigAlpha=0.05
n=length(xData)
R=max(which(xData<=sigAlpha*(1:n)/n))
cData$BHclass=rep(0,n)
cData$BHclass[1:R]=1
cData[1:150,]
table(cData[,2:3])
################### classification MM #############################################
classTh=0.9
BetaUnifClass=as.integer((1-gammaMM)*dbeta(xData, alphaMM, betaMM)/(gammaMM+(1-gammaMM)*dbeta(xData, alphaMM, betaMM))>classTh)
cData=as.data.frame(cbind(xData,BetaUnifClass))
sigAlpha=0.05
n=length(xData)
R=max(which(xData<=sigAlpha*(1:n)/n))
cData$BHclass=rep(0,n)
cData$BHclass[1:R]=1
cData[1:55,2:3]


#############################################################################################################################
library(nleqslv)
library(mnorm)
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

##########################################################################################################################
trialsB=1000

alphaSig=0.05
za2=qnorm(1-alphaSig/2)
estPrecision=0.0001
maxIterations=1000
equationError=10^{-6}
maxEquationIterations=1000

bsMLE=bootstrapMLE(alphaT,betaT,gammaT,n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations)
alpha0=alphaT
beta0=betaT
gamma0=gammaT
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
z0MLE[1]=qnorm(sum(bsMLE$alpha<alphaT)/trialsBB)
z0MLE[2]=qnorm(sum(bsMLE$beta<betaT)/trialsBB)
z0MLE[3]=qnorm(sum(bsMLE$gamma<gammaT)/trialsBB)
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

alphaT_BCA=quantile(bsMLE$alpha,probs = c(ai1MLE[1],ai2MLE[1]))
betaT_BCA=quantile(bsMLE$beta,probs = c(ai1MLE[2],ai2MLE[2]))
gammaT_BCA=quantile(bsMLE$gamma,probs = c(ai1MLE[3],ai2MLE[3]))

alphaT_CI=quantile(bsMLE$alpha,probs = c(alphaSig/2,1-alphaSig/2))
betaT_CI=quantile(bsMLE$beta,probs = c(alphaSig/2,1-alphaSig/2))
gammaT_CI=quantile(bsMLE$gamma,probs = c(alphaSig/2,1-alphaSig/2))


#################################################################################    
# Bootstrap MPS
bsMPS=bootstrapMPS(alphaTmps,betaTmps,gammaTmps,n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations)
alpha0=alphaTmps
beta0=betaTmps
gamma0=gammaTmps

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
alphaTmps_BCA=quantile(bsMPS$alpha,probs = c(ai1MPS[1],ai2MPS[1]))
betaTmps_BCA=quantile(bsMPS$beta,probs = c(ai1MPS[2],ai2MPS[2]))
gammaTmps_BCA=quantile(bsMPS$gamma,probs = c(ai1MPS[3],ai2MPS[3]))

alphaTmps_CI=quantile(bsMPS$alpha,probs = c(alphaSig/2,1-alphaSig/2))
betaTmps_CI=quantile(bsMPS$beta,probs = c(alphaSig/2,1-alphaSig/2))
gammaTmps_CI=quantile(bsMPS$gamma,probs = c(alphaSig/2,1-alphaSig/2))


#################################################################################    
# Bootstrap MM
bsMM=bootstrapMM(alphaMM,betaMM,gammaMM,n,trialsB,estPrecision,maxIterations,equationError,maxEquationIterations)
alpha0=alphaMM
beta0=betaMM
gamma0=gammaMM
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
alphaMM_BCA=quantile(bsMM$alpha,probs = c(ai1MM[1],ai2MM[1]))
betaMM_BCA=quantile(bsMM$beta,probs = c(ai1MM[2],ai2MM[2]))
gammaMM_BCA=quantile(bsMM$gamma,probs = c(ai1MM[3],ai2MM[3]))

alphaMM_CI=quantile(bsMM$alpha,probs = c(alphaSig/2,1-alphaSig/2))
betaMM_CI=quantile(bsMM$beta,probs = c(alphaSig/2,1-alphaSig/2))
gammaMM_CI=quantile(bsMM$gamma,probs = c(alphaSig/2,1-alphaSig/2))


############################### Table ###########################
rbind(c(round(alphaT,4),round(alphaTmps,4),round(alphaMM,4)),
      c(paste("(",round(alphaT_BCA[1],2),",",round(alphaT_BCA[2],2),")",sep=""),
        paste("(",round(alphaTmps_BCA[1],2),",",round(alphaTmps_BCA[2],2),")",sep=""),
        paste("(",round(alphaMM_BCA[1],2),",",round(alphaMM_BCA[2],2),")",sep="")),
      
      c(round(betaT,4),round(betaTmps,4),round(betaMM,4)),
      c(paste("(",round(betaT_BCA[1],2),",",round(betaT_BCA[2],2),")",sep=""),
        paste("(",round(betaTmps_BCA[1],2),",",round(betaTmps_BCA[2],2),")",sep=""),
        paste("(",round(betaMM_BCA[1],2),",",round(betaMM_BCA[2],2),")",sep="")),
      
      c(round(gammaT,4),round(gammaTmps,4),round(gammaMM,4)),
      c(paste("(",round(gammaT_BCA[1],2),",",round(gammaT_BCA[2],2),")",sep=""),
        paste("(",round(gammaTmps_BCA[1],2),",",round(gammaTmps_BCA[2],2),")",sep=""),
        paste("(",round(gammaMM_BCA[1],2),",",round(gammaMM_BCA[2],2),")",sep="")) )

