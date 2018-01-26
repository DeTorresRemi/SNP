################
###LIBRAIRIES###
################

library("MASS")
library("DEoptim")

###############
###FONCTIONS###
###############

#calcul de X'beta
index_linear<-function (data, y.col, beta){
  nb.param<-length(beta)
  nb.obs<-nrow(data)
  X<-as.matrix(cbind(rep(1,nb.obs), data[,-y.col]))
  beta<-matrix(beta, nrow=nb.param, ncol=1)
  Xbeta<-X %*% beta
  return(Xbeta)
}

#calcul de logF
logF<-function(x, type){
  if (type=="probit") return(log(pnorm(x)))
  else if (type=="logit") return(log(1/(1+exp(-x))))
}

#calcul de log(1-F)
logcF<-function(x, type){
  if (type=="probit") return(log(1-pnorm(x)))
  else if (type=="logit") return(log(1-1/(1+exp(-x))))
}

#calcul de la contribution de vraisemblance pour chaque individu
indiv_llh<-function(data, y.col, beta, type){
  n<-nrow(data)
  result<-matrix(NA, nrow = n, ncol = 1)
  for (i in (1:n)){
    if (data[i,y.col]==0) result[i,1]<-logcF(index_linear(data, y.col, beta)[i], type)
    else if (data[i,y.col]==1) result[i,1]<-logF(index_linear(data, y.col, beta)[i], type)
  }
  return(result)
}

#calcul de la contribution minimale de vraisemblance pour chaque individu pour y=0 ou y=1
min_indiv_llh<-function(data, y.col, beta, type){
  indiv_llh0<-logcF(index_linear(data, y.col,beta), type)
  indiv_llh1<-logF(index_linear(data, y.col,beta), type)
  indiv_llh01<-cbind(indiv_llh0,indiv_llh1)
  return(apply(indiv_llh01, FUN=min, MARGIN=1))
}

#calcul de la vraisemblance MSTLE
MSTLE_llh<-function(data, y.col, beta, lambda, type){
  n<-dim(data)[1]
  hn<-floor(lambda*n)
  if (hn==n){
    llh<--sum(indiv_llh(data, y.col, beta, type))
    observations_removed<-"none"
  }
    else{
      #extraire indices des hn contributions de vraisemblance maximales
      index.max<-vector(mode="numeric",length = hn)
      min_indiv_llh0<-min_indiv_llh(data, y.col, beta, type)
      for (i in (1:hn)){
        index.max[i]<-which.max(min_indiv_llh0)
        min_indiv_llh0[index.max[i]]<--Inf
      }
      llh<--sum(indiv_llh(data[index.max,], y.col, beta, type))
      observations_removed<-(1:n)[-index.max]
    }
  output<-list(llh, observations_removed)
  names(output)<-c("MSTLE_llh", "index of obs. removed from llh")
  return(output)
}

#optimisation (minimisation par DEoptim de: -llh) de beta avec les paramètres préconisés par CIZEK
optim_MSTLE<-function(data, y.col, lambda, type){
    data<<-data
    y.col<<-y.col
    lambda<<-lambda
    type<<-type
    P<-30*(ncol(data)-1)
    lower<<-rep(-3,length=ncol(data))
    upper<<--lower
    par1<<-DEoptim.control(CR=1,F=1,strategy=1,itermax=100,parallelType = 1, parVar = list("data","y.col","lambda","type","lower","upper","par1","MSTLE_llh", "min_indiv_llh", "logF", "logcF", "indiv_llh", "index_linear"))
    par2<<-DEoptim.control(CR=1,F=0.75, strategy=1,itermax=100,parallelType = 1, parVar = list("data","y.col","lambda","type","lower","upper","par2","MSTLE_llh", "min_indiv_llh", "logF", "logcF", "indiv_llh", "index_linear"))
    par3<<-DEoptim.control(CR=1,F=0.5, strategy=1,itermax=100,parallelType = 1, parVar = list("data","y.col","lambda","type","lower","upper","par3","MSTLE_llh", "min_indiv_llh", "logF", "logcF", "indiv_llh", "index_linear"))
    results<-vector("list", length=length(lambda))
    set.seed(1234)
  for (i in 1:length(lambda)){
    if (P<=100) result<-DEoptim(function(beta) MSTLE_llh(data, y.col, beta, lambda[i], type)[[1]],lower, upper, control = par1)
      else if (P>100 & P<=250) result<-DEoptim(function(beta) MSTLE_llh(data, y.col, beta, lambda[i], type)[[1]],lower, upper, control = par2)
      else if (P>250) result<-DEoptim(function(beta) MSTLE_llh(data, y.col, beta, lambda[i], type)[[1]], lower, upper, control = par3)
  
  results[[i]]<-result$optim$bestmem
  names(results[[i]])<-c("intercept", colnames(data)[-y.col])
  }
  return(results)
}

#AMSTLE
AMSTLE<-function(data, y.col, lambda, type){
  T1<-Sys.time()
  D<-length(lambda)
  n<-nrow(data)
  hn=floor(lambda*n)
  beta_MSTLE<-optim_MSTLE(data, y.col, lambda, type)
  avrg.llh<-vector("numeric", length = D)
  for (d in (1:D)){
    avrg.llh[d]<--MSTLE_llh(data, y.col, beta_MSTLE[[d]], lambda[d], type)[[1]]/hn[d]
    }
  d.star<-which.max(avrg.llh)
  T2<-Sys.time()
  Tdiff<-difftime(T2, T1)
  names(beta_MSTLE)<-paste("with", hn, "obs. kept")
  output<-list(beta_MSTLE[[d.star]], 
               lambda[d.star], 
               hn[d.star],
               n-hn[d.star],
               MSTLE_llh(data, y.col, beta_MSTLE[[d.star]], lambda[d.star], type)[[2]], 
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
                   "temps d'exécution")
  return(output)
}

###########
###TESTS###
###########

#calcul de X'beta
#index_linear(data,y.col,beta)

#calcul de la contribution de vraisemblance pour chaque individu en semi-param
#indiv_llh(data, y.col,beta, type)

#calcul sp de la contribution minimale de vraisemblance pour chaque individu pour y=0 ou y=1
#min_indiv_llh(data, y.col, beta, type)

#calcul de la vraisemblance MSTLE
#MSTLE_llh(data, y.col, beta, lambda, type)

#optimisation (minimisation par DEoptim de: -llh) de beta avec les paramètres préconisés par CIZEK
#optim_MSTLE(data, y.col, lambda, type)

#AMSTLE
#AMSTLE(data, y.col, lambda, type)
