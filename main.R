library("xtable")

source("C:/Users/Remi/Documents/ENSAE/Semi and non parametric econometric/SNP/AMSTLE_parametric.R") #à modifier
source("C:/Users/Remi/Documents/ENSAE/Semi and non parametric econometric/SNP/AMSTLE_semiparametric.R") #à modifier

#DONNEES
data<-leuk

#TRAITEMENT
data[,1]<-data[,1]/1000
data[,2]<-as.numeric(data[,2])-1
data[which(data[,3]<=52),3]<-0
data[which(data[,3]>52),3]<-1
data

#PARAMETRES
y.col=3 #colonne de la variable expliquée binaire Y
lambda<-c(2/3,0.7,0.75,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1)
type<-"logit" #ou "probit"

#RESULTATS
beta_AMSTLE<-AMSTLE(data, y.col, lambda, type) #environ 3-4min
beta_AMSTLE_sp<-AMSTLE_sp(data, y.col, lambda) #très long à faire tourner

#TABLEAUX LATEX

#tableaux des coeff
Beta<-t(data.frame(beta_AMSTLE$coeff))
xtable(Coeff_parametrique, caption = "Beta estimé par l'AMSTLE paramétrique sur les données leuk")

Beta_sp<-t(data.frame(beta_AMSTLE_sp$coeff))
xtable(Beta_sp, caption = "Beta estimé par l'AMSTLE semi-paramétrique sur les données leuk")

#obs. retirées
xtable(data[beta_AMSTLE$`obs. removed at optimum lambda`,], digits = 0, caption = "Observations returées du calcul de la vraisemblance par l'AMSTLE paramétrique")

#graphe des average llh
av_llh<-cbind(lambda,beta_AMSTLE$`average llh`)
plot(av_llh, type = "b", pch = 2, xlab = "lambda", ylab = "log-vrais. moy.", sub = "Log-vrais. moyennes du MSTLE paramétrique par lambda")

av_llh<-cbind(lambda,beta_AMSTLE_sp$`average llh`)
plot(av_llh, type = "b", pch = 2, xlab = "lambda", ylab = "log-vrais. moy.", sub = "Log-vrais. moyennes du MSTLE paramétrique par lambda")


#####################
###RESULTATS PARAM###
#####################

# > beta_AMSTLE

# $coeff
# intercept        wbc         ag 
# 0.2119169 -0.2354463  2.5580640 
# 
# $`optimum lambda`
# [1] 0.925
# 
# $`nb. of obs. kept at optimum lambda`
# [1] 30
# 
# $`nb. of obs. removed at optimum lambda`
# [1] 3
# 
# $`obs. removed at optimum lambda`
# [1] 17 32 33
# 
# $`max. average llh`
# [1] -0.3185038
# 
# $`average llh`
# [1] -0.4339707 -0.3872124 -0.3980195 -0.3675014 -0.3538931 -0.3412541 -0.3412541
# [8] -0.3294867 -0.3185038 -0.3961827 -0.3838020 -0.4706355
# 
# $`MSTLE beta estimates`
# $`MSTLE beta estimates`$`with 22 obs. kept`
# intercept        wbc         ag 
# 0.2017635 -0.2321370  2.5409009 
# 
# $`MSTLE beta estimates`$`with 23 obs. kept`
# intercept           wbc            ag 
# -1.897032e+00  1.061708e-07 -5.824026e-06 
# 
# $`MSTLE beta estimates`$`with 24 obs. kept`
# intercept        wbc         ag 
# 0.2079765 -0.2342337  2.5520259 
# 
# $`MSTLE beta estimates`$`with 26 obs. kept`
# intercept        wbc         ag 
# 0.2116177 -0.2353845  2.5578957 
# 
# $`MSTLE beta estimates`$`with 27 obs. kept`
# intercept        wbc         ag 
# 0.2119169 -0.2354463  2.5580640 
# 
# $`MSTLE beta estimates`$`with 28 obs. kept`
# intercept        wbc         ag 
# 0.2119170 -0.2354463  2.5580637 
# 
# $`MSTLE beta estimates`$`with 28 obs. kept`
# intercept        wbc         ag 
# 0.2119170 -0.2354463  2.5580640 
# 
# $`MSTLE beta estimates`$`with 29 obs. kept`
# intercept        wbc         ag 
# 0.2119172 -0.2354463  2.5580643 
# 
# $`MSTLE beta estimates`$`with 30 obs. kept`
# intercept        wbc         ag 
# 0.2119169 -0.2354463  2.5580640 
# 
# $`MSTLE beta estimates`$`with 31 obs. kept`
# intercept           wbc            ag 
# 1.377066e+00 -1.938804e-01 -5.587935e-09 
# 
# $`MSTLE beta estimates`$`with 32 obs. kept`
# intercept           wbc            ag 
# 1.377099e+00 -1.938842e-01 -4.190952e-09 
# 
# $`MSTLE beta estimates`$`with 33 obs. kept`
# intercept         wbc          ag 
# -1.30734464 -0.03177171  2.26106578 
# 
# $`temps d'exécution`
# Time difference of 10.12492 mins

##########################
###RESULTATS SEMI-PARAM###
##########################

# > beta_AMSTLE_sp

# $coeff
# intercept        wbc         ag 
# -2.9438045 -0.1774942  2.8994133 
# 
# $`optimum lambda`
# [1] 1
# 
# $`nb. of obs. kept at optimum lambda`
# [1] 33
# 
# $`nb. of obs. removed at optimum lambda`
# [1] 0
# 
# $`obs. removed at optimum lambda`
# [1] "none"
# 
# $`max. average llh`
# [1] -0.2497761
# 
# $`average llh`
# [,1]
# [1,] -0.5777819
# [2,] -0.5635166
# [3,] -0.5487560
# [4,] -0.5127154
# [5,] -0.4962245
# [6,] -0.4812040
# [7,] -0.4812282
# [8,] -0.4667979
# [9,] -0.4540327
# [10,] -0.4396012
# [11,] -0.4258900
# [12,] -0.2497761
# 
# $`MSTLE beta estimates`
# $`MSTLE beta estimates`$`with 22 obs. kept`
# intercept         wbc          ag 
# -0.01199432  0.01812401 -0.99049193 
# 
# $`MSTLE beta estimates`$`with 23 obs. kept`
# intercept         wbc          ag 
# 2.64110297 -0.05408776  2.95573720 
# 
# $`MSTLE beta estimates`$`with 24 obs. kept`
# intercept         wbc          ag 
# -2.25073157  0.04752211 -2.59197412 
# 
# $`MSTLE beta estimates`$`with 26 obs. kept`
# intercept         wbc          ag 
# 1.69762458 -0.03176747  1.73629818 
# 
# $`MSTLE beta estimates`$`with 27 obs. kept`
# intercept         wbc          ag 
# 1.07117349  0.04363521 -2.38178607 
# 
# $`MSTLE beta estimates`$`with 28 obs. kept`
# intercept         wbc          ag 
# -2.51778848  0.04204614 -2.29816565 
# 
# $`MSTLE beta estimates`$`with 28 obs. kept`
# intercept         wbc          ag 
# 1.71943787 -0.02518891  1.37627758 
# 
# $`MSTLE beta estimates`$`with 29 obs. kept`
# intercept         wbc          ag 
# -1.01386292  0.02578659 -1.40924136 
# 
# $`MSTLE beta estimates`$`with 30 obs. kept`
# intercept         wbc          ag 
# -0.64662620  0.04671087 -2.55316645 
# 
# $`MSTLE beta estimates`$`with 31 obs. kept`
# intercept         wbc          ag 
# 1.76179438  0.05211763 -2.84854000 
# 
# $`MSTLE beta estimates`$`with 32 obs. kept`
# intercept         wbc          ag 
# -1.91058415 -0.04079748  2.22896358 
# 
# $`MSTLE beta estimates`$`with 33 obs. kept`
# intercept        wbc         ag 
# -2.9438045 -0.1774942  2.8994133 
# 
# 
# $`time of execution`
# Time difference of 40.25887 mins
