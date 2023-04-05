library(ciftiTools)
library(Rcmdr)
library(factoextra)
library(matconv)
library(rpubs)
library(lattice)
library(corrplot)
library(matlib)
library(wavelets)
library(Matrix)
library(R.matlab)
library(MASS)
library(corrplot)
library(Rfast2)
library(colf)
library(invgamma)
library(mvtnorm)
library(bayesplot)
library(rstanarm)
library(ggplot2)
library(spate)
library(tidyverse)
library(XML)
library(methods)
library(optimbase)
library(progress)
library(LaplacesDemon)
library(neuRosim)
library(pheatmap)
library(RColorBrewer)
library(colorRamps)
library(gtools)


##############################################     RUN IF RestingStateAnalysis.R HASN"T BEEN RUN ############################
ciftiTools.setOption("wb_path", "C:/Users/Owner/Documents/Graduate Documents/Research/HCP")


task <- read_xifti("~/Graduate Documents/Research/HCP/Data/100307/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii",
                   brainstructures =  "all")

parc <- load_parc("Yeo_17")


parc_vec <- c(as.matrix(parc)) 
parc_n <- 114

task <- move_from_mwall(task, NA) #replace the medial wall vertices with NA so that the length of the parcellation vector will match the # of cortical data vertices

sub_keys <- as.numeric(task$meta$subcort$labels) - 2
sub_keys <- parc_n + sub_keys 
brain_vec <- c(parc_vec)
brain_vec <- c(brain_vec, sub_keys)


#Need to center data - scale = FALSE causes only the mean to be subtracted
task_mat <- t(as.matrix(task))
task_mat_C <- scale(task_mat, scale = FALSE)
task_mat_C[is.na(task_mat_C)] = 0

t = 284

##############################################################################################################################
#####################################            TASK STATE SVD - Michelle Way              ##################################
##############################################################################################################################


F_PC_task <- vector(mode = "list", length = parc_n + 19)
SigmaFF_task <- vector(mode = "list", length = parc_n + 19)


V_right <- vector(mode = "list", length = parc_n + 19)
V_right_sized <- vector(mode = "list", length = parc_n + 19)
dimROI <- c()
for (p in 1:(parc_n + 19)){  
  
  taskdata_p <- task_mat_C[,brain_vec == p] #for each parcel p, we get a logical vector which essentially takes all of the corresponding T index rows of xii_mat and outputs them, matrix will be of size #Trueforparcelp x 1200 (tseries)
  dimROI[p] <-  dim(taskdata_p)[2]
  nr = dim(taskdata_p)[2]
 
  full_svd <- svd(taskdata_p/sqrt(t), nv = min(t,nr), nu = min(t,nr)) 
  sing_vals <- full_svd$d  #DIM: vector of length min(1200,nr) -> here will be nr
  V_right[[p]] <- full_svd$v   #DIM p x nv  -> here will be nr x nr
  sing_vals_matrix <- matrix(diag(sing_vals), nrow = min(nr,t), ncol = min(nr,t))  #DIM p x p
  var <- cumsum(sing_vals^2)/sum(sing_vals^2)
  index <- seq(from = 1, to = length(var), by = 1)
  num_comp <- min(index[var >= 0.5])
  V_right_sized[[p]] <- V_right[[p]][,1:num_comp]
  F_PC_task[[p]] <- taskdata_p %*% V_right[[p]][,1:num_comp] #DIM: t x nr %*% nr x num_comp -> t x num_comp 
  SigmaFF_task[[p]] <- t(F_PC_task[[p]])%*%F_PC_task[[p]]/t 
  
}

########################## Formulating Y_L and PHI #################################

phi <- bdiag(V_right_sized)

phi <- matrix(phi, nrow = 91146, ncol = 1836)

phi_t = t(phi)

Y_L <- c()

for(i in 1:length(F_PC_task)){
  Y_L <- cbind(Y_L, F_PC_task[[i]])
}

######################### #PC on final matrix ##############
nr = dim(Y_L)[2] #t x 1836 sum_pk

full_svd <- svd(Y_L/sqrt(t), nv = min(t,nr), nu = min(t,nr)) 
sing_vals <- full_svd$d 

sing_vals_matrix <- matrix(diag(sing_vals), nrow = t, ncol = t)
var <- cumsum(sing_vals^2)/sum(sing_vals^2)
index <- seq(from = 1, to = length(var), by = 1)
num_comp <- min(index[var >= 0.95])  

Y_G <- (Y_L %*% full_svd$v[,1:num_comp]) 
Y_G_forDWT <- t(Y_G) #remember to note that it is transposed for W function in matlab
s = dim(Y_G)[2]

psi <- full_svd$v[,1:num_comp]
psi_t <- t(psi)

write.matrix(Y_G_forDWT, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/Y_G.csv") 

#############################################################################################################################
#####################################                    DESIGN MATRIX                     ##################################
#############################################################################################################################


EVFunction <- function(onset, duration, TR){
  temp_indices = vector(mode = "list", length = length(onset))
  end_onset <- c()
  
  for(j in 1:length(onset)){
    end_onset[j] <- onset[j] + duration[j]
    
    shift_seq <- seq(from = TR, to = 214, by = TR) 
    zero_seq <- seq(from = 0, to = 214-TR, by = TR)
    
    ninterv <- length(shift_seq)
    interval = matrix(c(zero_seq, shift_seq), nrow = ninterv, ncol = 2, byrow = FALSE)
    
    temp_interval <- matrix(,nrow = ninterv, ncol = 2)
    
    for(i in 1:ninterv){
      temp_interval[i,1] = (onset[j] <= interval[i,2]) && (onset[j] > interval[i,1])
      temp_interval[i,2] = (end_onset[j] <= interval[i,2]) && (end_onset[j] > interval[i,1])
    }
    
    
    onset_indices1 <- which(temp_interval[,1] == TRUE)
    onset_indices2 <- which(temp_interval[,2] == TRUE)
    
    mid_onset1 <- interval[onset_indices1,2] - TR/2
    mid_onset2 <- interval[onset_indices2,2] - TR/2
    
    vol_indicator <- onset[j] > mid_onset1
    vol_indicator_end <- end_onset[j] > mid_onset2
    
    temp_indices[[j]] = seq(from = onset_indices1 + 1 - vol_indicator, to = onset_indices2 - 1 + vol_indicator_end, by = 1)
    
  }
  
  return(temp_indices)
  
}

#HRF
param1 <- 6
param2 <- 16
scaling_param <- 1/6

hrf <- function(t){
  hrf = dgamma(t, shape = param1) - scaling_param * dgamma(t, param2)
  return(hrf)
}

hrf_deriv <- function(x){
  hrf_deriv = (((x^15)*exp(-x))/7846046208000) - (((x^14)*exp(-x))/523069747200) - (((x^5)*exp(-x))/120) + (((x^4)*exp(-x))/24)
  return(hrf_deriv)
}



filename = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/100307/MNINonLinear/Results/tfMRI_MOTOR_LR/EVs"

t = 284
p = 5 
dt_aux <- seq(from = 0.36, to = 214, by = 0.72) 


hrf_deriv_eval <- c()
for(i in 1:t){
  hrf_deriv_eval[i] = hrf_deriv(dt_aux[i])
}


BoxCar_Func <- matrix(rep(0,times = t*p), nrow = t, ncol = p)

EV_list <- vector(mode = "list") 
EV_list[[1]] <- "lf.txt" #Left Foot
EV_list[[2]] <- "lh.txt" #Left Hand
EV_list[[3]] <- "rf.txt" #Right Foot
EV_list[[4]] <- "rh.txt" #Right Hand
EV_list[[5]] <- "t.txt" #Tongue

for(i in 1:p){
  EV = EV_list[[i]]
  EVPath = paste(filename,EV, sep = "/") 
  
  EV_Scan = matrix(scan(EVPath), nrow = 2, ncol = 3, byrow = TRUE)
  
  EV_Output = EVFunction(onset = EV_Scan[,1], duration = EV_Scan[,2], TR = 0.72)
  Boxcar_Indices = unlist(EV_Output)
  
  BoxCar_Func[Boxcar_Indices,i] = 1
}


#Convolution

Beta_HRF = matrix(rep(0, times = t*p), nrow = t, ncol = p)
Beta_dHRF = matrix(rep(0, times = t*p), nrow = t, ncol = p)




for(i in 1:p){
  hrf_temp = convolve(hrf(dt_aux[1:284]), BoxCar_Func[,i])  
  dhrf_temp = convolve(hrf_deriv(dt_aux[1:284]), BoxCar_Func[,i])
  
  Beta_HRF[,i] = hrf_temp[1:284]
  Beta_dHRF[,i] = dhrf_temp[1:284]
}

DesignMatrix = matrix(,nrow = t, ncol = 2*p)

for(j in 1:p){
  DesignMatrix[, (j-1)*2 + 1] = Beta_HRF[,j]
  DesignMatrix[, (j-1)*2 + 2] = Beta_dHRF[,j]
}
DesignMatrix = scale(DesignMatrix, scale = FALSE)

write.table(DesignMatrix, "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/DesignMat.csv", col.names = FALSE, row.names = FALSE)



#############################################################################################################################
#####################################                    GET DATA FROM MATLAB              ##################################
#############################################################################################################################

W <- read.csv("C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/W.csv", header = FALSE)
tw = dim(W)[1]
W <- matrix(unlist(W), nrow = tw, ncol = t)
W_Inv <- t(W)

psi_hat <- read.csv("C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/psi_hat50.csv", header = FALSE)
psi_hat <- matrix(unlist(psi_hat), nrow = s, ncol = 50) #of dim T^w x 50
psi_hat <- t(psi_hat)


alpha_hat <- read.csv("C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/alpha_hat50.csv", header = FALSE)
alpha_hat <- matrix(unlist(alpha_hat), nrow = s, ncol = 50) #of dim T^w x 50
alpha_hat <- t(alpha_hat)

WaveletInfo <- xmlParse(file = "~/Graduate Documents/Research/HCP/Data/wavespecs_YG.xml")
WaveletInfoDf <- xmlToDataFrame("~/Graduate Documents/Research/HCP/Data/wavespecs_YG.xml")
rownames(WaveletInfoDf) <- c("ndim", "boundary", "wavelet", "compress", "nlevels", "K1", "K2", "K3", "K4", "K5", "K6", "K7",
                             "K8", "K9", "Tw", "J", "T")


############################################################################################################################
###############################              MCMC BETA POSTERIOR SAMPLING                ###################################
############################################################################################################################

Nite = 50

nlev = as.numeric(WaveletInfoDf["nlevels",])
J = nlev + 1
m_seq = seq(1,nlev,1)

Kj = c()
for(i in 1:J){Kj[i] = WaveletInfoDf[5+i,]}
Kj = as.numeric(Kj) #Size J x 1 - contains the # of compressed wavelet coefficients at each level


mvector = c()
for(i in 1:J){
  aux0 = i*ones(Kj[i],1)
  mvector = append(mvector,aux0)
}

a0 = 2 
b0 = 2

kv = 100
p = dim(DesignMatrix)[2]
I_p = diag(p)
I_e = diag(tw)

DesignMatrix<- scale(DesignMatrix) 
Xstar <- W%*%DesignMatrix
Ystar = W%*%Y_G

iter = 8000

betaIV_full <- matrix(NA, nrow = p, ncol = s)
betaMCMC <- array(NA, dim = c(iter, p,s))

alpha_s = matrix(NA, nrow = tw, ncol = s)
psi_s = matrix(NA, nrow = tw, ncol = s)


for(s in 1:167){
  print(s)
  alpha_s[,s] = rep(alpha_hat[50,s], times = tw)   #These are created and worked on in the MethodOfMoments.R file
  psi_s[,s] = rep(psi_hat[50,s], times = tw)
  
  SigmaWsInv = (psi_s[,s]^-1)*(2^(mvector*alpha_s[,s])) 
  
  SigmaWsMInv = diag(SigmaWsInv, nrow = tw, ncol = tw)
  
  
  
  # STEP 2 - Initialize B*_s 
  Bstar_0s <- inv(t(Xstar)%*%SigmaWsMInv%*%Xstar)%*%t(Xstar)%*%SigmaWsMInv%*%Ystar[,s]  #p x 1
  betaIV_full[,s] <- Bstar_0s
  
  #STEP 3 - Update Psi_s from full conditional 
  
  for(i in 1:iter){
    
    #print(i)
    
    
    M = t(Ystar[,s])%*%SigmaWsMInv%*%Ystar[,s] - 2*(t(Ystar[,s])%*%SigmaWsMInv%*%Xstar%*%betaIV_full[,s]) + t(betaIV_full[,s])%*%t(Xstar)%*%SigmaWsMInv%*%Xstar%*%betaIV_full[,s]
    aux_psi = rinvgamma(1, shape = a0 + ((tw/2))/2, b0 + ((M/2)/2))
    psi_s[,s] = rep(aux_psi, times = tw)
    SigmaWsInv = (psi_s[,s]^-1)*(2^(mvector*alpha_s[,s])) 
    
    SigmaWsMInv = diag(SigmaWsInv, nrow = tw, ncol = tw)
    
    #STEP 4 - Update Beta from full conditional
    
    mu = inv(t(Xstar)%*%SigmaWsMInv%*%Xstar%*%(I_p + kv^-1 * I_p))%*%t(Xstar)%*%SigmaWsMInv%*%Ystar[,s]
    sigma = t(Xstar)%*%SigmaWsMInv%*%Xstar%*%(I_p + kv^-1 * I_p)
    betaIV_full[,s] <- rmvnorm(1,mean = mu, sigma = sigma)
    
    
    betaMCMC[i,,s] <- betaIV_full[,s]
    betaIV_full[,s] <- t(betaIV_full[,s])
    
    
  }
}

beta_mean <- colMeans(betaMCMC_burnin) 

write.table(beta_mean, "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/beta_means.csv", col.names = FALSE, row.names = FALSE)



###########################################################################################################################
#####################################        THINNING/TOLERANCE/LAG on BETAs       ########################################
###########################################################################################################################


betaMCMC_burnin <- betaMCMC[1001:8000,,]
betaMCMC_thinned <- array(NA, dim = c(700,p,s))
burn_n <- length(betaMCMC_burnin[,1,1])
for(s in 1:167){
  
  for(p in 1:10){
    betaMCMC_thinned[,p,s] <- Thin(betaMCMC_burnin[,p,s], By = 10)
    
  }
}

beta_mean <- colMeans(betaMCMC_thinned)



write.table(beta_mean, "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/beta_means.csv", col.names = FALSE, row.names = FALSE)





############################################################################################################################
###############################         CALCULATING RESIDUALS & PROJECTING BACK          ###################################
############################################################################################################################


#Don't centralize the residual matrix
index_NA <- is.na(task_mat)
ind_NA <- which(index_NA, arr.ind = T)
ind_NA <- as.matrix(ind_NA, nrow = 1582448, ncol = 2)
NA_cols <- unique(ind_NA[,2])
NA_cols <- as.vector(NA_cols)
task_mat_noNA <- task_mat[,-NA_cols]


zeros_mat <- brain_vec == 0
ind_zero <- which(zeros_mat, arr.ind = TRUE)

taskdata_nozero <- task_mat_C[,-ind_zero]
taskdata_nozero_orig <- task_mat[,-ind_zero]


ColMeans <- colMeans(task_mat)
ColMeans_nozero <- ColMeans[-ind_zero]




W_inv = t(W) 


Y_G = W_inv%*%Ystar

X = W_inv%*%Xstar


Y_GGlobSpace <- Y_G%*%psi_t #This gives us Y_L space which is t x sum(pk) Y_L = Y~ * phi



betaGlobSpace <- beta_mean%*%psi_t #(p x s) x (s x sum(pk)) = 10 x 1836
betaLocSpace <- betaGlobSpace%*%phi_t
Xbeta_original <- DesignMatrix_Scaled%*%betaLocSpace

XBeta_GlobSpace <- X%*%betaGlobSpace #t x sum(pk) = 284 x 1836 = (t x p) x  (p x sum(pk))

Y_LocalSpace <- Y_L%*%phi_t

XBeta_LocalSpace <- XBeta_GlobSpace%*%phi_t

Resid <- taskdata_nozero - XBeta_LocalSpace


# XstarBetastar <- Xstar%*%beta_mean
# 
# ResidStar <- Ystar - XstarBetastar
# 
# ResidGlobSpace <- W_Inv%*%ResidStar%*%psi_t #t x sum(pk)
# 
# phi_t <- t(phi)
# 
# ResidLocalSpace <- ResidGlobSpace%*%phi_t




write.matrix(ColMeans_nozero , file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/ColMeans.csv")
write.matrix(taskdata_nozero , file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/Y_cent.csv")
write.matrix(taskdata_nozero_orig, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/Y_original.csv") 
write.matrix(brain_vec_no0, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/ROI_vec.csv") 
write.matrix(DesignMatrix_Scaled, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/X.csv") 
###########################################################################################################################
########################    APPLY PCA same way as RESTING STATE TO RESIDUALS ##############################################
###########################################################################################################################
F_PC_BC <- vector(mode = "list", length = parc_n + 19)
SigmaFF_BC <- vector(mode = "list", length = parc_n + 19)

brain_vec_no0 <- brain_vec[brain_vec != 0]


for (p in 1:(parc_n + 19)){  
  
  data_p <- Resid[,brain_vec_no0 == p] 
  nr = dim(data_p)[2]
  
  
  if(nr <= t){
    full_svd <- svd(data_p/sqrt(t), nv = min(t,nr), nu = min(t,nr)) 
    sing_vals <- full_svd$d  #DIM: vector of length min(t,nr) -> here will be nr
    V_right <- full_svd$v   #DIM p x nv  -> here will be nr x nr
    sing_vals_matrix <- matrix(diag(sing_vals), nrow = nr, ncol = nr)  #DIM p x p
    var <- cumsum(sing_vals^2)/sum(sing_vals^2)
    index <- seq(from = 1, to = length(var), by = 1)
    num_comp <- min(index[var >= 0.5])
    F_PC_BC[[p]] <- (data_p %*% V_right[,1:num_comp]) #DIM: t x nr %*% nr x num_comp 
    SigmaFF_BC[[p]] <- t(F_PC_BC[[p]])%*%F_PC_BC[[p]]/t 
  }
  else{
    full_svd <- svd(data_p/sqrt(t), nv = min(t,nr), nu = min(t,nr))  
    sing_vals <- full_svd$d #DIM: vector of length min(1200,nr) 
    U_left <- full_svd$u #DIM n x nu  
    sing_vals_matrix <- matrix(diag(sing_vals), nrow = t, ncol = t)
    var <- cumsum(sing_vals^2)/sum(sing_vals^2)
    index <- seq(from = 1, to = length(var), by = 1)
    num_comp <- min(index[var >= 0.5])
    U_left_numcomp <- matrix(data = U_left[,1:num_comp], nrow = t, ncol = num_comp )
    F_PC_BC[[p]] <- sqrt(t) * U_left_numcomp %*% (sing_vals_matrix[1:num_comp, 1:num_comp])
    
    SigmaFF_BC[[p]] <- sing_vals_matrix[1:num_comp, 1:num_comp]%*%sing_vals_matrix[1:num_comp, 1:num_comp]
    
  }
}

dim_ROI_Resid <- c()
for(i in 1:(parc_n +19)){
  dim_ROI_Resid[i] = dim(SigmaFF_BC[[i]])[1]
}





############################################################################################################################
############## RV COEFFICIENT CALCULATION ON RESIDUAL MATRIX (pretending i have right matrix for now) ######################
############################################################################################################################
RV_BC = matrix(,nrow = parc_n + 19, ncol = parc_n + 19)

t_RV_BC = matrix(, nrow = parc_n + 19, ncol = parc_n + 19)

t = 284

for(i in 2:(parc_n + 19)) {
  for(j in 1:i){
    if(i == 17 && j != 17){
      CF_ij = SigmaFF_BC[[i]]^(-1/2)%*%(t(F_PC_BC[[i]])%*%F_PC_BC[[j]])%*%(inv(SigmaFF_BC[[j]])^(1/2))
      CF_ji = (inv(SigmaFF_BC[[j]])^(1/2))%*%(t(F_PC_BC[[j]])%*%F_PC_BC[[i]])%*%(SigmaFF_BC[[i]])^(-1/2)
      CF_ii = SigmaFF_BC[[i]]^(-1/2)*(t(F_PC_BC[[i]])%*%F_PC_BC[[i]])*(SigmaFF_BC[[i]])^(-1/2)
      CF_jj = (inv(SigmaFF_BC[[j]])^(1/2))%*%(t(F_PC_BC[[j]])%*%F_PC_BC[[j]])%*%(inv(SigmaFF_BC[[j]])^(1/2))
    }
    else if(j == 17 && i != 17){
      CF_ij = (inv(SigmaFF_BC[[i]])^(1/2))%*%(t(F_PC_BC[[i]])%*%F_PC_BC[[j]])%*%(SigmaFF_BC[[j]])^(-1/2)
      CF_ji = SigmaFF_BC[[j]]^(-1/2)%*%(t(F_PC_BC[[j]])%*%F_PC_BC[[i]])%*%(inv(SigmaFF_BC[[i]])^(1/2))
      CF_ii = (inv(SigmaFF_BC[[i]])^(1/2))%*%(t(F_PC_BC[[i]])%*%F_PC_BC[[i]])%*%(inv(SigmaFF_BC[[i]])^(1/2))
      CF_jj = SigmaFF_BC[[j]]^(-1/2)*(t(F_PC_BC[[j]])%*%F_PC_BC[[j]])*(SigmaFF_BC[[j]])^(-1/2)
    }
    else if(i == 17 && j == 17){
      CF_ij = SigmaFF_BC[[i]]^(-1/2)*(t(F_PC_BC[[i]])%*%F_PC_BC[[j]])*(SigmaFF_BC[[j]])^(-1/2)
      CF_ji = SigmaFF_BC[[j]]^(-1/2)*(t(F_PC_BC[[j]])%*%F_PC_BC[[i]])*(SigmaFF_BC[[i]])^(-1/2)
      CF_ii = SigmaFF_BC[[i]]^(-1/2)*(t(F_PC_BC[[i]])%*%F_PC_BC[[i]])*(SigmaFF_BC[[i]])^(-1/2)
      CF_jj = SigmaFF_BC[[j]]^(-1/2)*(t(F_PC_BC[[j]])%*%F_PC_BC[[j]])*(SigmaFF_BC[[j]])^(-1/2)
      
    }
    else{

      CF_ij = (inv(SigmaFF_BC[[i]])^(1/2))%*%(t(F_PC_BC[[i]])%*%F_PC_BC[[j]])%*%(inv(SigmaFF_BC[[j]])^(1/2))
      CF_ji = (inv(SigmaFF_BC[[j]])^(1/2))%*%(t(F_PC_BC[[j]])%*%F_PC_BC[[i]])%*%(inv(SigmaFF_BC[[i]])^(1/2))
      CF_ii = (inv(SigmaFF_BC[[i]])^(1/2))%*%(t(F_PC_BC[[i]])%*%F_PC_BC[[i]])%*%(inv(SigmaFF_BC[[i]])^(1/2))
      CF_jj = (inv(SigmaFF_BC[[j]])^(1/2))%*%(t(F_PC_BC[[j]])%*%F_PC_BC[[j]])%*%(inv(SigmaFF_BC[[j]])^(1/2))
      
    }
    
    RV_BC[i,j] = tr(CF_ij%*%CF_ji)/sqrt(tr(CF_ii%*%CF_ii)*tr(CF_jj%*%CF_jj)) 
    RV_BC[j,i] = RV_BC[i,j]
    
    
    
    Beta_i <- (tr(SigmaFF_BC[[i]]))^2/tr(SigmaFF_BC[[i]]%*%SigmaFF_BC[[i]]) 
    Beta_j <- (tr(SigmaFF_BC[[j]]))^2/tr(SigmaFF_BC[[j]]%*%SigmaFF_BC[[j]]) 
    
    
    ERV <- sqrt(Beta_i * Beta_j) / (t - 1)
    
    tau_isub <- sum(diag((t(F_PC_BC[[i]])%*%F_PC_BC[[i]]))^2)/(tr((t(F_PC_BC[[i]]) %*% F_PC_BC[[i]])%*%(t(F_PC_BC[[i]]) %*% F_PC_BC[[i]])))
    tau_jsub <- sum(diag((t(F_PC_BC[[j]])%*%F_PC_BC[[j]]))^2)/(tr((t(F_PC_BC[[j]]) %*% F_PC_BC[[j]])%*%(t(F_PC_BC[[j]]) %*% F_PC_BC[[j]])))
    
    tau_i <- ((t - 1)/((t - 3)*(t-1 - Beta_i))) * (t*(t+1)*tau_isub - (t-1)*(Beta_i + 2))
    tau_j <- ((t - 1)/((t - 3)*(t-1 - Beta_j))) * (t*(t+1)*tau_jsub - (t-1)*(Beta_j + 2))
    
    VRV <- ((2*(t-1-Beta_i)*(t-1-Beta_j))/((t+1)*(t-1)^2 * (t-2)))* (1 + ((t-3)/(2*t*(t-1)))*tau_i*tau_j)
    t_RV_BC[i,j] <- (RV_BC[i,j] - ERV)/sqrt(VRV)
    
  }
}

alpha = 0.05
D = 133*133 

p = (1 - (alpha)/(D))
p1 = 1 - (alpha)/(133*132/2)
threshold <- qnorm(p)
threshold1 <- qnorm(p1)

Connectivity_50 <- abs(t_RV_BC) > threshold

Connectivity_50_1 <- abs(t_RV_BC) > threshold1

diag(Connectivity_50) <- NA
significant <- which(Connectivity_50, arr.ind = T)
sig_rows <- significant[,1]
sig_cols <- significant[,2]
tRV_sig_val <- c()
RV_sig_val <- c()


for(i in 1:2286){
  tRV_sig_val[i] <- t_RV_BC[sig_rows[i],sig_cols[i]]
  RV_sig_val[i] <- RV_BC[sig_rows[i], sig_cols[i]]
}



heatmap(RV_BC, Rowv = NA, Colv = NA, revC = TRUE)
heatmap(t_RV_BC, Rowv = NA, Colv = NA, revC = TRUE)

RV_BC[1,1] <- 1

rownames(RV_BC) <- as.character(c(1:133))
colnames(RV_BC) <- c(1:133)

for(i in 1:133){
  if(odd(as.numeric(colnames(RV_BC)[i])) == FALSE){
    colnames(RV_BC)[i] <- ""
  }
}

pheatmap(RV_BC, treeheight_row = 0, treeheight_col = 0, main = "Background Connectivity RV Matrix",
         fontsize = 8, angle_col = "270", show_rownames = F, color = blue2green2red(100), show_colnames = F)
pheatmap(RV_BC, treeheight_row = 0, treeheight_col = 0, main = "Background Connectivity RV",
         fontsize = 8, angle_col = "270", show_rownames = F, cluster_cols = F, cluster_rows = F, show_colnames = T,
         color = blue2green2red(100))

pheatmap(t_RV_BC, treeheight_row = 0, treeheight_col = 0, main = "Background Connectivity RV",
         fontsize = 8, angle_col = "270", show_rownames = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                     "RdBu")))(100))




write.matrix(significant, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/SigROI_BC.csv") 
write.matrix(tRV_sig_val, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/tRV_BC.csv")
write.matrix(RV_sig_val, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/RV_BC.csv")
# #############################################################################################################################
# #################                BACKGROUND CONNECTIVITY VIA SIGMA ROI MATRIX           #####################################
# #############################################################################################################################
# ***Not fully functional
# 
# 
# psi_post_mean = colMeans(psi_s)
# 
# 
# CovarianceTime = array(NA, dim = c(t,t,s))
# newpsi = c()
# newalpha = c()
# Sigma = c()
# for(s in 1:167){
#   newpsi = rep(psi_post_mean[s],times = tw)
#   newalpha = rep(alpha_hat[50,s], times = tw)
#   Sigma = (newpsi[s])*(2^(mvector*-newalpha[s]))
#   CovarianceTime[,,s] = W_Inv%*%diag(Sigma, nrow = tw, ncol = tw)%*%W
# }
# 
# 
# find_mode <- function(x) {
#   u <- unique(x)
#   tab <- tabulate(match(x, u))
#   u[tab == max(tab)]
# }
# 
# 
# Theta_s = c()
# rangeCov = c()
# for(s in 1:167){
#   rangeCov = range(diag(CovarianceTime[,,s]))
#   print(rangeCov)
#   Theta_s[s] = find_mode(diag(CovarianceTime[,,s]))
# }
# 
# Sigma_theta = diag(Theta_s, nrow = 167, ncol = 167)
# 
# SigmaROI = psi%*%Sigma_theta%*%psi_t
# 
# test <- psi%*%psi_t
# 
# F_PC_ROI <- vector(mode = "list", length = parc_n + 19)
# SigmaFF_ROI <- vector(mode = "list", length = parc_n + 19)
# dim_ROI = vector(mode = "numeric", length = parc_n + 19)
# 
# for(i in 1:(parc_n +19)){
#   dim_ROI[i] = dim(SigmaFF_task[[i]])[1]
# }
# 
# cumulativedim_ROI = cumsum(dim_ROI)
# 
# for (i in 1:(parc_n + 19)) {
#   if(i == 1){
#     F_PC_ROI[[i]] <- SigmaROI[,1:dim_ROI[i]]
#   }
#   else{
#     upp_bound = cumulativedim_ROI[i]
#     low_bound = cumulativedim_ROI[i-1]
#     F_PC_ROI[[i]] <- SigmaROI[,(low_bound + 1):upp_bound]
#   }
#   SigmaFF_ROI[[i]] <- t(F_PC_ROI[[i]])%*%F_PC_ROI[[i]]/t 
#   
# }
# 
# 
# for(i in 2:(parc_n + 19)) {
#   for(j in 1:i){
#     
#     if(i == 94){
#       M_ij = SigmaFF_ROI[[i]]^(-1/2)%*%(t(F_PC_ROI[[i]])%*%F_PC_ROI[[j]])%*%(inv(SigmaFF_ROI[[j]])^(1/2))/t
#       M_ji = (inv(SigmaFF_ROI[[j]])^(1/2))%*%(t(F_PC_ROI[[j]])%*%F_PC_ROI[[i]])%*%(SigmaFF_ROI[[i]])^(-1/2)/t
#       M_ii = SigmaFF_ROI[[i]]^(-1/2)*(t(F_PC_ROI[[i]])%*%F_PC_ROI[[i]])*(SigmaFF_ROI[[i]])^(-1/2)/t
#       M_jj = (inv(SigmaFF_ROI[[j]])^(1/2))%*%(t(F_PC_ROI[[j]])%*%F_PC_ROI[[j]])%*%(inv(SigmaFF_ROI[[j]])^(1/2))/t
#     }
#     if(j == 94){
#       M_ij = (inv(SigmaFF_ROI[[i]])^(1/2))%*%(t(F_PC_ROI[[i]])%*%F_PC_ROI[[j]])%*%(SigmaFF_ROI[[j]]^(-1/2))/t
#       M_ji = (SigmaFF_ROI[[j]]^(-1/2))%*%(t(F_PC_ROI[[j]])%*%F_PC_ROI[[i]])%*%(inv(SigmaFF_ROI[[i]])^(1/2))/t
#       M_ii = (inv(SigmaFF_ROI[[i]])^(1/2))%*%(t(F_PC_ROI[[i]])%*%F_PC_ROI[[i]])%*%(inv(SigmaFF_ROI[[i]])^(1/2))/t
#       M_jj = SigmaFF_ROI[[j]]^(-1/2)*(t(F_PC_ROI[[j]])%*%F_PC_ROI[[j]])*(SigmaFF_ROI[[j]])^(-1/2)/t
#     }
#     
#     else{
#       M_ij = (inv(SigmaFF_ROI[[i]])^(1/2))%*%(t(F_PC_ROI[[i]])%*%F_PC_ROI[[j]])%*%(inv(SigmaFF_ROI[[j]])^(1/2))/t
#       M_ji = (inv(SigmaFF_ROI[[j]])^(1/2))%*%(t(F_PC_ROI[[j]])%*%F_PC_ROI[[i]])%*%(inv(SigmaFF_ROI[[i]])^(1/2))/t
#       M_ii = (inv(SigmaFF_ROI[[i]])^(1/2))%*%(t(F_PC_ROI[[i]])%*%F_PC_ROI[[i]])%*%(inv(SigmaFF_ROI[[i]])^(1/2))/t
#       M_jj = (inv(SigmaFF_ROI[[j]])^(1/2))%*%(t(F_PC_ROI[[j]])%*%F_PC_ROI[[j]])%*%(inv(SigmaFF_ROI[[j]])^(1/2))/t
#     }
#     
#     
#     RV_BC[i,j] = tr(M_ij%*%M_ji)/sqrt(tr(M_ii%*%M_ii)*tr(M_jj%*%M_jj)) 
#     RV_BC[j,i] = RV_BC[i,j]
#     
#     
#     
#     Beta_i <- (tr(SigmaFF_ROI[[i]]))^2/tr(SigmaFF_ROI[[i]]%*%SigmaFF_ROI[[i]]) #This is just going to be 1 for any 1x1 matrix
#     Beta_j <- (tr(SigmaFF_ROI[[j]]))^2/tr(SigmaFF_ROI[[j]]%*%SigmaFF_ROI[[j]]) #Would be trace which would be diff't for other dims
#     
#     
#     ERV <- sqrt(Beta_i * Beta_j) / (t - 1)
#     
#     tau_isub <- sum(diag((t(F_PC_ROI[[i]])%*%F_PC_ROI[[i]]))^2)/(tr((t(F_PC_ROI[[i]]) %*% F_PC_ROI[[i]])%*%(t(F_PC_ROI[[i]]) %*% F_PC_ROI[[i]])))
#     tau_jsub <- sum(diag((t(F_PC_ROI[[j]])%*%F_PC_ROI[[j]]))^2)/(tr((t(F_PC_ROI[[j]]) %*% F_PC_ROI[[j]])%*%(t(F_PC_ROI[[j]]) %*% F_PC_ROI[[j]])))
#     
#     tau_i <- ((t - 1)/((t - 3)*(t-1 - Beta_i))) * (t*(t+1)*tau_isub - (t-1)*(Beta_i + 2))
#     tau_j <- ((t - 1)/((t - 3)*(t-1 - Beta_j))) * (t*(t+1)*tau_jsub - (t-1)*(Beta_j + 2))
#     
#     VRV <- ((2*(t-1-Beta_i)*(t-1-Beta_j))/((t+1)*(t-1)^2 * (t-2)))* (1 + ((t-3)/(2*t*(t-1)))*tau_i*tau_j)
#     t_RV_BC[i,j] <- (RV_BC[i,j] - ERV)/sqrt(VRV)
#     
#   }
}



