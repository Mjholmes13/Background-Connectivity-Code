library(Rfast2)
library(colf)
library(invgamma)
library(mvtnorm)
library(matlib)
library(bayesplot)
library(rstanarm)
library(ggplot2)
library(spate)
library(MASS)
library(Matrix)
library(tidyverse)
library(XML)
library(methods)
library(optimbase)

# WaveletInfo <- xmlParse(file = "~/Graduate Documents/Research/HCP/Data/wavespecs_YG.xml")
# WaveletInfoDf <- xmlToDataFrame("~/Graduate Documents/Research/HCP/Data/wavespecs_YG.xml")
# 
# rownames(WaveletInfoDf) <- c("ndim", "boundary", "wavelet", "compress", "nlevels", "K1", "K2", "K3", "K4", "K5", "K6", "K7",
#                              "K8", "K9", "Tw", "J", "T")
# #s = dim(Ystar)[2]
# 
# Nite = 50
# #tw = dim(W)[1]
# nlev = as.numeric(WaveletInfoDf["nlevels",])
# J = nlev + 1
# m_seq = seq(1,nlev,1)
# #psi_hat = matrix(NA, nrow = Nite, ncol = s) 
# #alpha_hat = matrix(NA, nrow = Nite, ncol = s)
# 
# Kj = c()
# for(i in 1:J){Kj[i] = WaveletInfoDf[5+i,]}
# Kj = as.numeric(Kj) #Size J x 1 - contains the # of compressed wavelet coefficients at each level
# 
# 
# mvector = c()
# for(i in 1:J){
#   aux0 = i*ones(Kj[i],1)
#   mvector = append(mvector,aux0)
# }

cumsumav = cumsum(Kj)

for(k in 1:s){
  Ehat_m[,m] = (I - Xstar%*%inv(t(Xstar)%*%inv(SigmaWs)%*%Xstar)%*%t(Xstar)%*%inv(SigmaWs))%*%Ystar[,s]
}

# Ehat_m = matrix(NA, nrow = tw , ncol = nlev )
# Theta_m = vector(mode = "numeric", length = nlev)
# 
# Ehat_s = matrix(NA, nrow = tw, ncol = 4)
# Theta_s = vector(mode = "numeric", length = 4)

alpha_s = 0
phi_s = 1
I = diag(tw)
#Fix s
s = 1

for(m in 1:nlev){
  SigmaWs = diag(phi_s*(2^(alpha_s))^(-m_seq[m]), nrow = tw, ncol = tw)
  Ehat_m[,m] = (I - Xstar%*%inv(t(Xstar)%*%inv(SigmaWs)%*%Xstar)%*%t(Xstar)%*%inv(SigmaWs))%*%Ystar[,s] #tw x 1
  Theta_m[m] = var(Ehat_m[,m])
}

Log_Theta = log(Theta_s)


