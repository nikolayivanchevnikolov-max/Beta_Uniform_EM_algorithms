rm(list=ls())
seedNum=2024

gammaVec=c(0.8,0.9,0.95)
alphaVec=c(0.1,0.3,0.5)
betaVec=c(1.2,1.5,2)
n=100

alpha=0.3
beta=1.2
gamma=0.9
trialsAll=1020
trialsB=100


##################################################################################################
############################## Table a = 0.1, 0.3, 0.5 ###########################################
##################################################################################################
resMLE=array(0, dim = c(6,3,3))
resMPS=array(0, dim = c(6,3,3))
resMM=array(0, dim = c(6,3,3))
for (i in 1:length(gammaVec)) {
  resMLE[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMLE_n_",n,"_seed_",seedNum,"_a_",alphaVec[i],"_b_",beta,"_gamma_",gamma,"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMPS[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMPS_n_",n,"_seed_",seedNum,"_a_",alphaVec[i],"_b_",beta,"_gamma_",gamma,"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMM[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMM_n_",n,"_seed_",seedNum,"_a_",alphaVec[i],"_b_",beta,"_gamma_",gamma,"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
}

tableA=cbind(resMLE[,1,1], resMPS[,1,1], resMM[,1,1],
                 resMLE[,1,2], resMPS[,1,2], resMM[,1,2],
                 resMLE[,1,3], resMPS[,1,3], resMM[,1,3])
rownames(tableA)=c("MSE","Bias","CP_BCa","CI_Length_BCa","CP","CI_Length")
colnames(tableA)=rep(c("ML","MPS","MM"),3)

##################################################################################################
############################## Table b = 1.2, 1.5, 2 #############################################
##################################################################################################
resMLE=array(0, dim = c(6,3,3))
resMPS=array(0, dim = c(6,3,3))
resMM=array(0, dim = c(6,3,3))
for (i in 1:length(gammaVec)) {
  resMLE[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMLE_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",betaVec[i],"_gamma_",gamma,"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMPS[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMPS_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",betaVec[i],"_gamma_",gamma,"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMM[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMM_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",betaVec[i],"_gamma_",gamma,"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
}

tableB=cbind(resMLE[,2,1], resMPS[,2,1], resMM[,2,1],
             resMLE[,2,2], resMPS[,2,2], resMM[,2,2],
             resMLE[,2,3], resMPS[,2,3], resMM[,2,3])
rownames(tableB)=c("MSE","Bias","CP_BCa","CI_Length_BCa","CP","CI_Length")
colnames(tableB)=rep(c("ML","MPS","MM"),3)

##################################################################################################
############################ Table gamma=0.8, 0.9, 0.95 ##########################################
##################################################################################################
resMLE=array(0, dim = c(6,3,3))
resMPS=array(0, dim = c(6,3,3))
resMM=array(0, dim = c(6,3,3))
for (i in 1:length(gammaVec)) {
  resMLE[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMLE_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMPS[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMPS_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMM[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMM_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
}

tableGamma=cbind(resMLE[,3,1], resMPS[,3,1], resMM[,3,1],
                  resMLE[,3,2], resMPS[,3,2], resMM[,3,2],
                  resMLE[,3,3], resMPS[,3,3], resMM[,3,3])
rownames(tableGamma)=c("MSE","Bias","CP_BCa","CI_Length_BCa","CP","CI_Length")
colnames(tableGamma)=rep(c("ML","MPS","MM"),3)
