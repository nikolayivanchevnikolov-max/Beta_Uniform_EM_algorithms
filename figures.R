rm(list=ls())
library(ggplot2)
library(latex2exp)
seedNum=2024

gammaVec=c(0.7,0.8,0.9,0.95)#
subInd=c(1,2,3,4)
n=100
alpha=0.3
beta=1.2
trialsAll=1020
trialsB=100

resMLE=array(0, dim = c(6,3,4))
resMPS=array(0, dim = c(6,3,4))
resMM=array(0, dim = c(6,3,4))
for (i in 1:length(gammaVec)) {
  resMLE[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMLE_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMPS[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMPS_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMM[,,i]=as.matrix(read.csv(paste("superSim\\res100\\resMM_n_",n,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
}

resMLE500=array(NA, dim = c(6,3,4))
resMPS500=array(NA, dim = c(6,3,4))
resMM500=array(NA, dim = c(6,3,4))
for (i in subInd) {
  resMLE500[,,i]=as.matrix(read.csv(paste("superSim\\res500\\resMLE_n_",500,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMPS500[,,i]=as.matrix(read.csv(paste("superSim\\res500\\resMPS_n_",500,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMM500[,,i]=as.matrix(read.csv(paste("superSim\\res500\\resMM_n_",500,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
}

resMLE1000=array(NA, dim = c(6,3,4))
resMPS1000=array(NA, dim = c(6,3,4))
resMM1000=array(NA, dim = c(6,3,4))
for (i in subInd) {
  resMLE1000[,,i]=as.matrix(read.csv(paste("superSim\\res1000\\resMLE_n_",1000,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMPS1000[,,i]=as.matrix(read.csv(paste("superSim\\res1000\\resMPS_n_",1000,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
  resMM1000[,,i]=as.matrix(read.csv(paste("superSim\\res1000\\resMM_n_",1000,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll,"_trialsB_",trialsB,".csv",sep="")))
}

resMLE5000=array(NA, dim = c(6,3,4))
resMPS5000=array(NA, dim = c(6,3,4))
resMM5000=array(NA, dim = c(6,3,4))
trialsAll=c(986,996,1086,1020)
for (i in subInd) {
  resMLE5000[,,i]=as.matrix(read.csv(paste("superSim\\res5000\\resMLE_n_",5000,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll[i],"_trialsB_",trialsB,".csv",sep="")))
  resMPS5000[,,i]=as.matrix(read.csv(paste("superSim\\res5000\\resMPS_n_",5000,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll[i],"_trialsB_",trialsB,".csv",sep="")))
  resMM5000[,,i]=as.matrix(read.csv(paste("superSim\\res5000\\resMM_n_",5000,"_seed_",seedNum,"_a_",alpha,"_b_",beta,"_gamma_",gammaVec[i],"_trials_",trialsAll[i],"_trialsB_",trialsB,".csv",sep="")))
}

#############################################################################################
####################################### LOG SCALE MSE ####################################### 
#############################################################################################

######################################## Parameter a ####################################### 
plotData=cbind.data.frame(x=1-gammaVec,y1=log(resMLE[1,1,]),y2=log(resMPS[1,1,]),y3=log(resMM[1,1,]),
                          y1_500=log(resMLE500[1,1,]),y2_500=log(resMPS500[1,1,]),y3_500=log(resMM500[1,1,]),
                          y1_1000=log(resMLE1000[1,1,]),y2_1000=log(resMPS1000[1,1,]),y3_1000=log(resMM1000[1,1,]),
                          y1_5000=log(resMLE5000[1,1,]),y2_5000=log(resMPS5000[1,1,]),y3_5000=log(resMM5000[1,1,]))
msePlotMLE <- ggplot(plotData, aes(x,y1)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=100$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMLE
ggsave(msePlotMLE,filename=paste("plots_All\\log_alphaMsePlot_All_n_100_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMPS<- ggplot(plotData, aes(x,y1_500)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_500, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_500, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=500$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMPS
ggsave(msePlotMPS,filename=paste("plots_All\\log_alphaMsePlot_All_n_500_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMM<- ggplot(plotData, aes(x,y1_1000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=1000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMM
ggsave(msePlotMM,filename=paste("plots_All\\log_alphaMsePlot_All_n_1000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMPS<- ggplot(plotData, aes(x,y1_5000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=5000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
# scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
#ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMPS
ggsave(msePlotMPS,filename=paste("plots_All\\log_alphaMsePlot_All_n_5000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)


######################################## Parameter b ####################################### 
plotData=cbind.data.frame(x=1-gammaVec,y1=log(resMLE[1,2,]),y2=log(resMPS[1,2,]),y3=log(resMM[1,2,]),
                          y1_500=log(resMLE500[1,2,]),y2_500=log(resMPS500[1,2,]),y3_500=log(resMM500[1,2,]),
                          y1_1000=log(resMLE1000[1,2,]),y2_1000=log(resMPS1000[1,2,]),y3_1000=log(resMM1000[1,2,]),
                          y1_5000=log(resMLE5000[1,2,]),y2_5000=log(resMPS5000[1,2,]),y3_5000=log(resMM5000[1,2,]))
msePlotMLE <- ggplot(plotData, aes(x,y1)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=100$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-3.5,4.25), breaks = seq(-3, 4, by = 1))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMLE
ggsave(msePlotMLE,filename=paste("plots_All\\log_betaMsePlot_All_n_100_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMPS<- ggplot(plotData, aes(x,y1_500)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_500, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_500, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=500$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-3.5,4.25), breaks = seq(-3, 4, by = 1))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMPS
ggsave(msePlotMPS,filename=paste("plots_All\\log_betaMsePlot_All_n_500_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMM<- ggplot(plotData, aes(x,y1_1000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=1000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-3.5,4.25), breaks = seq(-3, 4, by = 1))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMM
ggsave(msePlotMM,filename=paste("plots_All\\log_betaMsePlot_All_n_1000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMPS<- ggplot(plotData, aes(x,y1_5000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=5000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
# scale_y_continuous(limits=c(-3.5,4.25), breaks = seq(-3, 4, by = 1))+
#ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMPS
ggsave(msePlotMPS,filename=paste("plots_All\\log_betaMsePlot_All_n_5000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

######################################## Parameter gamma ####################################### 
plotData=cbind.data.frame(x=1-gammaVec,y1=log(resMLE[1,3,]),y2=log(resMPS[1,3,]),y3=log(resMM[1,3,]),
                          y1_500=log(resMLE500[1,3,]),y2_500=log(resMPS500[1,3,]),y3_500=log(resMM500[1,3,]),
                          y1_1000=log(resMLE1000[1,3,]),y2_1000=log(resMPS1000[1,3,]),y3_1000=log(resMM1000[1,3,]),
                          y1_5000=log(resMLE5000[1,3,]),y2_5000=log(resMPS5000[1,3,]),y3_5000=log(resMM5000[1,3,]))
msePlotMLE <- ggplot(plotData, aes(x,y1)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=100$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMLE
ggsave(msePlotMLE,filename=paste("plots_All\\log_gammaMsePlot_All_n_100_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMPS<- ggplot(plotData, aes(x,y1_500)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_500, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_500, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=500$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMPS
ggsave(msePlotMPS,filename=paste("plots_All\\log_gammaMsePlot_All_n_500_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMM<- ggplot(plotData, aes(x,y1_1000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=1000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMM
ggsave(msePlotMM,filename=paste("plots_All\\log_gammaMsePlot_All_n_1000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

msePlotMPS<- ggplot(plotData, aes(x,y1_5000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=5000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of Mean Squared Error for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
# scale_y_continuous(limits=c(-8.85,0), breaks = seq(0, -7.5, by = -2.5))+
#ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
msePlotMPS
ggsave(msePlotMPS,filename=paste("plots_All\\log_gammaMsePlot_All_n_5000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)
###############################################################################################

#############################################################################################
####################################### Relative Bias ####################################### 
#############################################################################################

######################################## Parameter a ####################################### 
plotData=cbind.data.frame(x=1-gammaVec,y1=log(abs(resMLE[2,1,]/alpha)),y2=log(abs(resMPS[2,1,]/alpha)),y3=log(abs(resMM[2,1,]/alpha)),
                          y1_500=log(abs(resMLE500[2,1,]/alpha)),y2_500=log(abs(resMPS500[2,1,]/alpha)),y3_500=log(abs(resMM500[2,1,]/alpha)),
                          y1_1000=log(abs(resMLE1000[2,1,]/alpha)),y2_1000=log(abs(resMPS1000[2,1,]/alpha)),y3_1000=log(abs(resMM1000[2,1,]/alpha)),
                          y1_5000=log(abs(resMLE5000[2,1,]/alpha)),y2_5000=log(abs(resMPS5000[2,1,]/alpha)),y3_5000=log(abs(resMM5000[2,1,]/alpha)))
biasPlot <- ggplot(plotData, aes(x,y1)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=100$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\alphaBiasPlot_All_n_100_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_500)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_500, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_500, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=500$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\alphaBiasPlot_All_n_500_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_1000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=1000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\alphaBiasPlot_All_n_1000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_5000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=5000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $a=", alpha, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
# scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
#ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\alphaBiasPlot_All_n_5000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

######################################## Parameter b ####################################### 
plotData=cbind.data.frame(x=1-gammaVec,y1=log(abs(resMLE[2,2,]/beta)),y2=log(abs(resMPS[2,2,]/beta)),y3=log(abs(resMM[2,2,]/beta)),
                          y1_500=log(abs(resMLE500[2,2,]/beta)),y2_500=log(abs(resMPS500[2,2,]/beta)),y3_500=log(abs(resMM500[2,2,]/beta)),
                          y1_1000=log(abs(resMLE1000[2,2,]/beta)),y2_1000=log(abs(resMPS1000[2,2,]/beta)),y3_1000=log(abs(resMM1000[2,2,]/beta)),
                          y1_5000=log(abs(resMLE5000[2,2,]/beta)),y2_5000=log(abs(resMPS5000[2,2,]/beta)),y3_5000=log(abs(resMM5000[2,2,]/beta)))
biasPlot <- ggplot(plotData, aes(x,y1)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=100$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\betaBiasPlot_All_n_100_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_500)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_500, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_500, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=500$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\betaBiasPlot_All_n_500_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_1000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=1000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\betaBiasPlot_All_n_1000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_5000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=5000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $b=", beta, "$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
# scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
#ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\betaBiasPlot_All_n_5000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

######################################## Parameter gamma ####################################### 
plotData=cbind.data.frame(x=1-gammaVec,y1=log(abs(resMLE[2,3,]/gammaVec)),y2=log(abs(resMPS[2,3,]/gammaVec)),y3=log(abs(resMM[2,3,]/gammaVec)),
                          y1_500=log(abs(resMLE500[2,3,]/gammaVec)),y2_500=log(abs(resMPS500[2,3,]/gammaVec)),y3_500=log(abs(resMM500[2,3,]/gammaVec)),
                          y1_1000=log(abs(resMLE1000[2,3,]/gammaVec)),y2_1000=log(abs(resMPS1000[2,3,]/gammaVec)),y3_1000=log(abs(resMM1000[2,3,]/gammaVec)),
                          y1_5000=log(abs(resMLE5000[2,3,]/gammaVec)),y2_5000=log(abs(resMPS5000[2,3,]/gammaVec)),y3_5000=log(abs(resMM5000[2,3,]/gammaVec)))
biasPlot <- ggplot(plotData, aes(x,y1)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=100$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\gammaBiasPlot_All_n_100_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_500)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_500, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_500, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_500, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_500, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=500$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  #ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\gammaBiasPlot_All_n_500_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_1000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_1000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_1000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_1000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=1000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
  # scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
  ##ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\gammaBiasPlot_All_n_1000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

biasPlot<- ggplot(plotData, aes(x,y1_5000)) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=2)+
  geom_point(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=8, shape=1,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=10, shape=0,stroke = 2.5)+
  geom_point(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=8, shape=5,stroke = 2.5)+
  geom_line(data = plotData, aes(x = x, y = y1_5000, color="blue2"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y2_5000, color="seagreen4"),size=1.5,linetype="dashed")+
  geom_line(data = plotData, aes(x = x, y = y3_5000, color="red2"),size=1.5,linetype="dashed")+
  theme_bw(base_size = 30)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 1),
        legend.position='top',
        legend.text = element_text(size = 35),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.size =  unit(0.625, "in"))+
  scale_color_identity(name = TeX("$n=5000$,     Method:"),
                       breaks = c( "blue2","seagreen4","red2"),
                       labels = c( "ML","MPS","MM"),
                       guide = "legend")+
  xlab(TeX("$1-\\gamma$"))+ylab(TeX(paste("Log of MARB for $\\gamma$",sep="")))+
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))#+
# scale_y_continuous(limits=c(-8.85,0.75), breaks = seq(0, -7.5, by = -2.5))+
#ggtitle(TeX(paste("Parameters: $a=",alpha, "$, $b=", beta, "$",sep="")))
biasPlot
ggsave(biasPlot,filename=paste("plots_All\\gammaBiasPlot_All_n_5000_seed_",seedNum,"_a_",alpha,"_b_",beta,"_trials_",trialsAll[1],"_trialsB_",trialsB,".pdf",sep=""), 
       width = 12, height = 9)

