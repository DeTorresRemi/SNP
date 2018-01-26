################
###LIBRAIRIES###
################

library("MASS")
library("DEoptim")
library("np")

source("C:/Users/Remi/Documents/ENSAE/Semi and non parametric econometric/SNP/AMSTLE_parametric.R") #à modifier

###############
###FONCTIONS###
###############

#estimations de l'estimateur à noyau
kernel_estimates<-function(data, y.col, beta){
  G_hat<-npreg(data = data, txdat = index_linear(data, y.col, beta), tydat = data[,y.col])
  return(G_hat$mean)
}

#calcul de la contribution de vraisemblance pour chaque individu en semi-param
indiv_llh_sp<-function(data, y.col, beta, vector_kernel_estimates){
  n<-nrow(data)
  result<-matrix(NA, nrow = n, ncol = 1)
  for (i in (1:n)){
    if (data[i,y.col]==0) result[i,1]<-log(1-vector_kernel_estimates[i])
    else if (data[i,y.col]==1) result[i,1]<-log(vector_kernel_estimates[i])
  }
  return(result)
}

#calcul sp de la contribution minimale de vraisemblance pour chaque individu pour y=0 ou y=1
min_indiv_llh_sp<-function(data, y.col, beta, vector_kernel_estimates){
  indiv_llh0<-log(1-vector_kernel_estimates)
  indiv_llh1<-log(vector_kernel_estimates)
  indiv_llh01<-cbind(indiv_llh0,indiv_llh1)
  return(apply(indiv_llh01, FUN=min, MARGIN=1))
}

#calcul de la vraisemblance MSTLE
MSTLE_llh_sp<-function(data, y.col, beta, lambda){
  n<-dim(data)[1]
  hn<-floor(lambda*n)
  vector_kernel_estimates<-kernel_estimates(data, y.col, beta)
  if (hn==n){
    llh<--sum(indiv_llh_sp(data, y.col, beta, vector_kernel_estimates))
    observations_removed<-"none"
  } 
  else{
    #extraire indices des hn contributions de vraisemblance maximales
    index.max<-vector(mode="numeric",length = hn)
    min_indiv_llh0<-min_indiv_llh_sp(data, y.col, beta, vector_kernel_estimates)
    for (i in (1:hn)){
      index.max[i]<-which.max(min_indiv_llh0)
      min_indiv_llh0[index.max[i]]<--Inf
    }
    llh<--sum(indiv_llh_sp(data[index.max,], y.col, beta, vector_kernel_estimates))
    observations_removed<-(1:n)[-index.max]
  }
  output<-list(llh, observations_removed)
    names(output)<-c("MSTLE_llh", "index of obs. removed from llh")
    return(output)
}

#optimisation (minimisation par DEoptim de: -llh) de beta avec les paramètres préconisés par CIZEK
optim_MSTLE_sp<-function(data, y.col, lambda){
  P<-30*(ncol(data)-1)
  lower<<-rep(-3,length=ncol(data))
  upper<<--lower
  par1<<-DEoptim.control(CR=1, F=1, strategy=1, itermax=100, parallelType = 1, packages = list("np"), parVar = list("data","y.col","lambda","type","lower","upper","par1","MSTLE_llh_sp", "min_indiv_llh_sp", "kernel_estimates", "indiv_llh_sp", "index_linear", "npreg"))
  par2<<-DEoptim.control(CR=1, F=0.75, strategy=1, itermax=100, parallelType = 1, packages = list("np"), parVar = list("data","y.col","lambda","type","lower","upper","par1","MSTLE_llh_sp", "min_indiv_llh_sp", "kernel_estimates", "indiv_llh_sp", "index_linear", "npreg"))
  par3<<-DEoptim.control(CR=1, F=0.5, strategy=1, itermax=100, parallelType = 1, packages = list("np"), parVar = list("data","y.col","lambda","type","lower","upper","par1","MSTLE_llh_sp", "min_indiv_llh_sp", "kernel_estimates", "indiv_llh_sp", "index_linear", "npreg"))
  results<-vector("list", length=length(lambda))
  set.seed(1234)
  for (i in 1:length(lambda)){
    if (P<=100) result<-DEoptim(function(beta) MSTLE_llh_sp(data, y.col, beta, lambda[i])[[1]],lower, upper, control = par1)
    else if (P>100 & P<=250) result<-DEoptim(function(beta) MSTLE_llh_sp(data, y.col, beta, lambda[i])[[1]],lower, upper, control = par2)
    else if (P>250) result<-DEoptim(function(beta) MSTLE_llh_sp(data, y.col, beta, lambda[i])[[1]], lower, upper, control = par3)
    
    results[[i]]<-result$optim$bestmem
    names(results[[i]])<-c("intercept", colnames(data)[-y.col])
  }
  return(results)
}

#AMSTLE
AMSTLE_sp<-function(data, y.col, lambda){
  T1<-Sys.time()
  D<-length(lambda)
  n<-nrow(data)
  P<-30*(ncol(data)-1)
  hn=floor(lambda*n)
  beta_MSTLE<-optim_MSTLE_sp(data, y.col, lambda)
  avrg.llh<-matrix(NA, ncol=1, nrow=D)
  for (d in (1:D)){
    avrg.llh[d]<--MSTLE_llh_sp(data, y.col, beta_MSTLE[[d]], lambda[d])[[1]]/hn[d]
  }
  d.star<-which.max(avrg.llh)
  T2<-Sys.time()
  Tdiff<-difftime(T2, T1)
  names(beta_MSTLE)<-paste("with", hn, "obs. kept")
  output<-list(beta_MSTLE[[d.star]], 
               lambda[d.star], 
               hn[d.star],
               n-hn[d.star],
               MSTLE_llh_sp(data, y.col, beta_MSTLE[[d.star]], lambda[d.star])[[2]], 
               avrg.llh[d.star], 
               avrg.llh, 
               beta_MSTLE,
               Tdiff)
  names(output)<-c("coeff", 
                   "optimum lambda",
                   "nb. of obs. kept at optimum lambda",
                   "nb. of obs. removed at optimum lambda", 
                   "obs. removed at optimum lambda", 
                   "max. average llh", 
                   "average llh", 
                   "MSTLE beta estimates",
                   "time of execution")
  return(output)
}

###########
###TESTS###
###########

#estimations de l'estimateur à noyau
#kernel_estimates(data, y.col, beta)

#calcul de la contribution de vraisemblance pour chaque individu en semi-param
#indiv_llh_sp(data, y.col,beta)

#calcul sp de la contribution minimale de vraisemblance pour chaque individu pour y=0 ou y=1
#min_indiv_llh_sp(data, y.col, beta)

#calcul de la vraisemblance MSTLE
#MSTLE_llh_sp(data, y.col, beta, lambda)

#optimisation (minimisation par DEoptim de: -llh) de beta avec les paramètres préconisés par CIZEK
#optim_MSTLE_sp(data, y.col, lambda)

#AMSTLE
#AMSTLE_sp(data, y.col, lambda)
