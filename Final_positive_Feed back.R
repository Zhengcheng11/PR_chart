#----Phase-II: the residual method is same as the existing Shewhart
rm(list=ls())
library(Rlab)
library(MASS)
library(iZID)
library(foreach)
library(parallel)
library(doParallel)

kernal=detectCores(logical = F)
cl=makeCluster(kernal)
registerDoParallel(cl)
#################################################################
#------------------Parameter setting----------------------------#
#################################################################
N=10000;p=0.005;ARL0=1/p;
c_gamma=c(0.05,0.1,0.2)  #the smooth parameter of the EWMA chart
lambda0=Y0=4;lambda0=Y0=4;

beta0=0.266;beta1=0.3;alpha1=0.3;

#M:the sample size of Phase-I; M.Mon:the sample size of Phase-II; B:the Monte Carlo simulation times
M=1000;M.Mon=1000;B=1000
L.low=h.low=0
L.up=h.up=20

#################################################################
#------------------Function definition----------------------------#
#################################################################
data.gen=function(cY0,clambda0,beta0,beta1,alpha1,M,N){
  
  M.y=M.lambda=matrix(NA,N,M+1)
  M.y[,1]=cY0
  M.lambda[,1]=clambda0
  
  for(i in 2:(M+1)){
    M.lambda[,i]=exp(beta0+beta1*log(M.y[,i-1]+1)+alpha1*log(M.lambda[,i-1]))
    M.y[,i]=rpois(N,M.lambda[,i])
  }
  
  M_y=M.y[,2:(M+1)]
  M_lambda=M.lambda[,2:(M+1)]
  
  return(list("M.Y"=M_y,"M.lambda"=M_lambda))
}

PR_Shewhart.RL=function(cy,clam,beta0,beta1,alpha1,L1){
  
  PR_cy=(cy-clam)/sqrt(clam)
  u.PR=mean(PR_cy)
  sd.PR=sd(PR_cy)
  
  UCL=u.PR+L1*sd.PR #is a value
  
  f.PR_cy=c(PR_cy,max(UCL)+1)
  
  t=which(f.PR_cy>UCL)[1]
  
  return(list("t"=t,"u.PR"=u.PR,"sd.PR"=sd.PR)) 
}

PR_EWMA_RL=function(cy,clam,gamma,h){
  
  PR_CY=(cy-clam)/sqrt(clam)
  
  plot.stat=t=0
  Z=c()
  zi=0
  
  while(plot.stat<h & t<length(cy)){
    t=t+1
    zi=max(0,(1-gamma)*zi+gamma*PR_CY[t])
    Z[t]=zi
    plot.stat=Z[t]
  }
  return(t)
}

Pro.e.Mlambda=function(mon.My,IC.cy,IC.clam,beta0,beta1,alpha1){
  
  Obser.y=cbind(IC.cy,mon.My)
  CM=dim(Obser.y)[2]
  e.Mlambda=matrix(NA,dim(Obser.y)[1],CM)
  e.Mlambda[,1]=IC.clam
  
  for(i in 2:CM){
    e.Mlambda[,i]=exp(beta0+beta1*log(Obser.y[,i-1]+1)+alpha1*log(e.Mlambda[,i-1]))
  }
  
  f.Mlambda=e.Mlambda[,-1]   #the estimation of lambda
  
  return(f.Mlambda)
}

Pro.Shewhart.RL1=function(mon.cy,beta0,beta1,alpha1,e.clam,p,M.Mon,B){
  
  Sim.data=apply(matrix(e.clam,length(e.clam),1),1, function(xx) rpois(B,xx))
  
  Sort.y=apply(Sim.data,2, function(x) sort(x))  #dim(Sim.data)=[B,length(e.clam)]
  ch=Sort.y[(1-p)*B,]  #the UCL of the prob Shewhart chart,is a vector
  
  f.ch=c(ch,0)
  f.mon.cy=c(mon.cy,max(ch))
  
  RL=which(f.mon.cy>f.ch)[1]
  
  return(RL)
}

PR_Shewhart.RL1=function(cy,clam,beta0,beta1,alpha1,UCL){
  
  PR_cy=(cy-clam)/sqrt(clam)
  
  f.PR_cy=c(PR_cy,max(UCL)+1)
  
  t=which(f.PR_cy>UCL)[1]
  
  return(t) 
}

###------------------h.find-------------------------###
#PR_Shewhart
PR_Shewhart.L=function(data,ARL0,L.low,L.up){
  
  for(i in 1:100){
    
    L=(L.up+L.low)/2
    
    clusterEvalQ(cl, library(VGAM))
    clusterExport(cl,varlist=c("data","PR_Shewhart.RL","L"), envir=environment())
    #PR_Shewhart.RL(cy,clam,beta0,beta1,alpha1,L1)
    Result=parApply(cl,data,1,function(xx) PR_Shewhart.RL(xx[,1],xx[,2],beta0,beta1,alpha1,L))
    
    RL=t(sapply(Result, getElement, "t"))
    u.PR=t(sapply(Result, getElement, "u.PR"))
    sd.PR=t(sapply(Result, getElement, "sd.PR"))
    
    if(abs(mean(RL)-ARL0)<0.2 | L.up-L.low<0.00002){
      return(list("L"=L,"u.PR"=u.PR,"sd.PR"=sd.PR))
      break
    }else if(mean(RL)-ARL0>0.2){
      L.up=L
    }else{L.low=L}
    print(c(L,mean(RL)))
  }
  
}

PR_EWMA.h=function(data,gamma,ARL0,h.low,h.up){
  
  for(i in 1:100){
    
    h=(h.up+h.low)/2
    
    clusterEvalQ(cl, library(VGAM))
    clusterExport(cl,varlist=c("data","PR_EWMA_RL","h","gamma"), envir=environment())
    #PR_EWMA_RL(cy,clam,gamma,h)
    RL=parApply(cl,data,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],gamma,h))  
    
    if(abs(mean(RL)-ARL0)<0.2 | h.up-h.low<0.00002){
      return(h)
      break
    }else if(mean(RL)-ARL0>0.2){
      h.up=h
    }else{h.low=h}
    print(c(h,mean(RL)))
  }
  return(h)
}
#################################################################
###-----------------------Implement------------------------------###
#################################################################
# IC data generation
cY0=rep(Y0,N);clambda0=rep(lambda0,N)

IC.data=data.gen(cY0,clambda0,beta0,beta1,alpha1,M,N)
M_y=IC.data$M.Y
M_lambda=IC.data$M.lambda

IC.cy=M_y[,M]
IC.clam=M_lambda[,M]

F.data=array(NA,dim=c(N,M,2))
F.data[,,1]=M_y
F.data[,,2]=M_lambda

PR_result=PR_Shewhart.L(F.data,ARL0,L.low,L.up)
L1=PR_result$L  
u.PR=mean(PR_result$u.PR)  
sd.PR=mean(PR_result$sd.PR) 
UCL=u.PR+L1*sd.PR
UCL

c_h.ew=c()
for(i in 1:length(c_gamma)){
  gamma=c_gamma[i]
  c_h.ew[i]=PR_EWMA.h(F.data,gamma,ARL0,h.low,h.up)
}

#####################################################################
###--------------------Monitoring---------------------------------###
#####################################################################
delta_beta0=c(0,0.02,0.05,0.1,0.15,0.2,0.3,0.5,0.7,1)
delta_beta1=delta_alpha1=c(0,0.02,0.05,0.07,0.1,0.15,0.2,0.3)

b0_OC.data=e.b0_lambda=array(NA,dim=c(N,M.Mon,length(delta_beta0)))
b1_OC.data=a1_OC.data=e.b1_lambda=e.a1_lambda=array(NA,dim=c(N,M.Mon,length(delta_beta1)))

#OC data by beta0
for(i in 1:length(delta_beta0)){
  d_beta0=delta_beta0[i]
  #data.gen(cY0,clambda0,beta0,beta1,alpha1,M,N)
  OC.data=data.gen(IC.cy,IC.clam,beta0+d_beta0,beta1,alpha1,M.Mon,N)
  b0_OC.data[,,i]=OC.data$M.Y
  #Pro.e.Mlambda(mon.My,IC.cy,IC.clam,beta0,beta1,alpha1)
  e.b0_lambda[,,i]=Pro.e.Mlambda(b0_OC.data[,,i],IC.cy,IC.clam,beta0,beta1,alpha1)
}

#OC data by beta1
for(i in 1:length(delta_beta1)){
  d_beta1=delta_beta1[i]
  OC.data_b=data.gen(IC.cy,IC.clam,beta0,beta1+d_beta1,alpha1,M.Mon,N)
  b1_OC.data[,,i]=OC.data_b$M.Y
  e.b1_lambda[,,i]=Pro.e.Mlambda(b1_OC.data[,,i],IC.cy,IC.clam,beta0,beta1,alpha1)
  
  d_alpha1=delta_alpha1[i]
  OC.data_a=data.gen(IC.cy,IC.clam,beta0,beta1,alpha1+d_alpha1,M.Mon,N)
  a1_OC.data[,,i]= OC.data_a$M.Y
  e.a1_lambda[,,i]=Pro.e.Mlambda(a1_OC.data[,,i],IC.cy,IC.clam,beta0,beta1,alpha1)
}


T_beta0= T_beta0_2=matrix(NA,length(delta_beta0),5)
T_beta1=T_beta1_2=T_alpha1=T_alpha1_2=matrix(NA,length(delta_alpha1),5)

###--------- Beta0 shift-----------------###
for(i in 1:length(delta_beta0)){
  mon.My=b0_OC.data[,,i]
  e.Mlambda=e.b0_lambda[,,i]
  
  Pro.M=array(NA,dim=c(dim(e.Mlambda)[1],dim(e.Mlambda)[2],2))
  Pro.M[,,1]=mon.My
  Pro.M[,,2]=e.Mlambda
  
  clusterEvalQ(cl, library(VGAM))
  clusterExport(cl,varlist=c("Pro.Shewhart.RL1","PR_Shewhart.RL1","PR_EWMA_RL","data.gen","Pro.M","c_gamma","beta0","beta1","alpha1","p","B","M.Mon","UCL","c_h.ew"), envir=environment())
  #Pro.Shewhart.RL1(mon.cy,beta0,beta1,alpha1,e.clam,p,M.Mon,B)
  Pro.RL=parApply(cl,Pro.M,1,function(xx) Pro.Shewhart.RL1(xx[,1],beta0,beta1,alpha1,xx[,2],p,M.Mon,B))
  
  #PR_Shewhart.RL1(cy,clam,beta0,beta1,alpha1,UCL)
  PR.RL=parApply(cl,Pro.M,1,function(xx) PR_Shewhart.RL1(xx[,1],xx[,2],beta0,beta1,alpha1,UCL))
  
  #PR_EWMA_RL(cy,clam,gamma,h)
  EWMA.RL_1=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[1],c_h.ew[1]))
  EWMA.RL_2=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[2],c_h.ew[2]))
  EWMA.RL_3=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[3],c_h.ew[3]))
  
  T_beta0[i,]=c(mean(Pro.RL),mean(PR.RL),mean(EWMA.RL_1),mean(EWMA.RL_2),mean(EWMA.RL_3))
  T_beta0_2[i,]=c(sd(Pro.RL),sd(PR.RL),sd(EWMA.RL_1),sd(EWMA.RL_2),sd(EWMA.RL_3))
  
  print(T_beta0[i,])
  #print(T_beta0_2[i,])
}

###-----------Beta 1 shift--------------###
for(i in 1:length(delta_beta1)){
  mon.My=b1_OC.data[,,i]
  e.Mlambda=e.b1_lambda[,,i]
  
  Pro.M=array(NA,dim=c(dim(e.Mlambda)[1],dim(e.Mlambda)[2],2))
  Pro.M[,,1]=mon.My
  Pro.M[,,2]=e.Mlambda
  
  clusterEvalQ(cl, library(VGAM))
  clusterExport(cl,varlist=c("Pro.Shewhart.RL1","PR_Shewhart.RL1","PR_EWMA_RL","data.gen","Pro.M","c_gamma","beta0","beta1","alpha1","p","B","M.Mon","UCL","c_h.ew"), envir=environment())
  #Pro.Shewhart.RL1(mon.cy,beta0,beta1,alpha1,e.clam,p,M.Mon,B)
  Pro.RL=parApply(cl,Pro.M,1,function(xx) Pro.Shewhart.RL1(xx[,1],beta0,beta1,alpha1,xx[,2],p,M.Mon,B))
  
  #PR_Shewhart.RL1(cy,clam,beta0,beta1,alpha1,UCL)
  PR.RL=parApply(cl,Pro.M,1,function(xx) PR_Shewhart.RL1(xx[,1],xx[,2],beta0,beta1,alpha1,UCL))
  
  #PR_EWMA_RL(cy,clam,gamma,h)
  EWMA.RL_1=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[1],c_h.ew[1]))
  EWMA.RL_2=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[2],c_h.ew[2]))
  EWMA.RL_3=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[3],c_h.ew[3]))
  
  T_beta1[i,]=c(mean(Pro.RL),mean(PR.RL),mean(EWMA.RL_1),mean(EWMA.RL_2),mean(EWMA.RL_3))
  T_beta1_2[i,]=c(sd(Pro.RL),sd(PR.RL),sd(EWMA.RL_1),sd(EWMA.RL_2),sd(EWMA.RL_3))
  
  print(T_beta1[i,])
  #print(T_beta1_2[i,])
}

###----------alpha 1 shift-------------###
for(i in 1:length(delta_alpha1)){
  mon.My=a1_OC.data[,,i]
  e.Mlambda=e.a1_lambda[,,i]
  
  Pro.M=array(NA,dim=c(dim(e.Mlambda)[1],dim(e.Mlambda)[2],2))
  Pro.M[,,1]=mon.My
  Pro.M[,,2]=e.Mlambda
  
  clusterEvalQ(cl, library(VGAM))
  clusterExport(cl,varlist=c("Pro.Shewhart.RL1","PR_Shewhart.RL","PR_EWMA_RL","data.gen","Pro.M","c_gamma","beta0","beta1","alpha1","p","B","M.Mon","UCL","c_h.ew","L1"), envir=environment())
  #Pro.Shewhart.RL1(mon.cy,beta0,beta1,alpha1,e.clam,p,M.Mon,B)
  Pro.RL=parApply(cl,Pro.M,1,function(xx) Pro.Shewhart.RL1(xx[,1],beta0,beta1,alpha1,xx[,2],p,M.Mon,B))
  
  #PR_Shewhart.RL1(cy,clam,beta0,beta1,alpha1,UCL)
  PR.RL=parApply(cl,Pro.M,1,function(xx) PR_Shewhart.RL1(xx[,1],xx[,2],beta0,beta1,alpha1,UCL))
  
  #PR_EWMA_RL(cy,clam,gamma,h)
  EWMA.RL_1=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[1],c_h.ew[1]))
  EWMA.RL_2=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[2],c_h.ew[2]))
  EWMA.RL_3=parApply(cl,Pro.M,1,function(xx) PR_EWMA_RL(xx[,1],xx[,2],c_gamma[3],c_h.ew[3]))
  
  T_alpha1[i,]=c(mean(Pro.RL),mean(PR.RL),mean(EWMA.RL_1),mean(EWMA.RL_2),mean(EWMA.RL_3))
  T_alpha1_2[i,]=c(sd(Pro.RL),sd(PR.RL),sd(EWMA.RL_1),sd(EWMA.RL_2),sd(EWMA.RL_3))
  
  print(T_alpha1[i,])
  #print(T_alpha1_2[i,])
}

Feed.neg_beta0_b1_a1=rbind(T_beta0,T_beta1,T_alpha1)
Feed.neg_beta0_b1_a1_2=rbind(T_beta0_2,T_beta1_2,T_alpha1_2)

rownames(Feed.neg_beta0_b1_a1)=c(as.character(delta_beta0),as.character(delta_beta1),as.character(delta_alpha1))
colnames(Feed.neg_beta0_b1_a1)=c("Shewhart","PR_Shewhart","EWMA_0.05","EWMA_0.1","EWMA_0.2")

rownames(Feed.neg_beta0_b1_a1_2)=c(as.character(delta_beta0),as.character(delta_beta1),as.character(delta_alpha1))
colnames(Feed.neg_beta0_b1_a1_2)=c("Shewhart","PR_Shewhart","EWMA_0.05","EWMA_0.1","EWMA_0.2")

write.csv(Feed.neg_beta0_b1_a1,file="C:/Users/25037/Desktop/MZC/TW Paper 2/Simulation Result N=10000/Feed.poi_Pro.csv")
write.csv(Feed.neg_beta0_b1_a1_2,file="C:/Users/25037/Desktop/MZC/TW Paper 2/Simulation Result N=10000/Feed.poi_Pro_2.csv")

stopCluster(cl)
