rm(list = ls())
library(ciftiTools)
library(Rcmdr)
library(factoextra)
library(matconv)
library(rpubs)
library(lattice)
library(corrplot)
library(matlib)
library(FactoMineR)
library(plot.matrix)
library(dplyr)
library(gplots)
library(colorRamps)


##############################################################################################################################
#########################################               LOAD CIFTI AND VARS                ###################################
##############################################################################################################################
ciftiTools.setOption("wb_path", "C:/Users/Owner/Documents/Graduate Documents/Research/HCP")


xii <- read_xifti("~/Graduate Documents/Research/HCP/Data/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii",
                  brainstructures =  "all")

parc <- load_parc("Yeo_17")
parc_vec <- c(as.matrix(parc))
parc_n <- 114


sub_keys <- as.numeric(xii$meta$subcort$labels) - 2
sub_keys <- parc_n + sub_keys
brain_vec <- c(parc_vec)
brain_vec <- c(brain_vec, sub_keys)

xii <- move_from_mwall(xii, NA) #replace the medial wall vertices with NA so that the length of the parcellation vector will match the # of cortical data vertices

#Center Data - only mean to be subtracted
xii_mat <- t(as.matrix(xii))
xii_mat_C <- scale(xii_mat, scale = FALSE)
xii_mat_C[is.na(xii_mat_C)] = 0

t = 1200


##############################################################################################################################
#####################################               RESTING STATE SVD                       ##################################
##############################################################################################################################
F_PC <- vector(mode = "list", length = parc_n + 19)
SigmaFF <- vector(mode = "list", length = parc_n + 19)
test <- c()

for (p in 1:(parc_n + 19)){
  
  data_p <- xii_mat_C[,brain_vec == p] #for each parcel p, we get a logical vector which essentially takes all of the corresponding T index rows of xii_mat and outputs them, matrix will be of size #Trueforparcelp x 1200 (tseries)
  
  nr = dim(data_p)[2]
  test[p] <- nr
  
  
  if(nr <= 1200){
    full_svd <- svd(data_p/sqrt(t), nv = min(t,nr), nu = min(t,nr))
    sing_vals <- full_svd$d  #DIM: vector of length min(1200,nr) -> here will be nr
    V_right <- full_svd$v   #DIM p x nv  -> here will be nr x nr
    sing_vals_matrix <- matrix(diag(sing_vals), nrow = nr, ncol = nr)  #DIM p x p
    var <- cumsum(sing_vals^2)/sum(sing_vals^2)
    index <- seq(from = 1, to = length(var), by = 1)
    num_comp <- min(index[var >= 0.5])
    F_PC[[p]] <- (data_p %*% V_right[,1:num_comp]) #DIM: 1200 x nr %*% nr x num_comp -> 1200 x num_comp
    SigmaFF[[p]] <- t(F_PC[[p]])%*%F_PC[[p]]/t
  }
  else{
    full_svd <- svd(data_p/sqrt(t), nv = min(t,nr), nu = min(t,nr))
    sing_vals <- full_svd$d #DIM: vector of length min(1200,nr) -> here will be 1200
    U_left <- full_svd$u #DIM n x nu  -> here will be 1200 x 1200
    sing_vals_matrix <- matrix(diag(sing_vals), nrow = t, ncol = t)
    var <- cumsum(sing_vals^2)/sum(sing_vals^2)
    index <- seq(from = 1, to = length(var), by = 1)
    num_comp <- min(index[var >= 0.5])
    F_PC[[p]] <- sqrt(t) * U_left[,1:num_comp] %*% (sing_vals_matrix[1:num_comp, 1:num_comp])

    SigmaFF[[p]] <- sing_vals_matrix[1:num_comp, 1:num_comp]%*%sing_vals_matrix[1:num_comp, 1:num_comp]
  
  }
}


dim_ROI_Rest <- c()
for(i in 1:(parc_n +19)){
  dim_ROI_Rest[i] = dim(SigmaFF[[i]])[1]
}



##############################################################################################################################
#####################################               RV COEFFICIENT CALC                     ##################################
##############################################################################################################################

RV = matrix(,nrow = parc_n + 19, ncol = parc_n + 19)

t_RV = matrix(, nrow = parc_n + 19, ncol = parc_n + 19)
t = 1200

check = 0
for(i in 2:(parc_n + 19)) {
  for(j in 1:i){
    
    CF_ij = (inv(SigmaFF[[i]])^(1/2))%*%(t(F_PC[[i]])%*%F_PC[[j]])%*%(inv(SigmaFF[[j]])^(1/2))
    CF_ji = (inv(SigmaFF[[j]])^(1/2))%*%(t(F_PC[[j]])%*%F_PC[[i]])%*%(inv(SigmaFF[[i]])^(1/2))
    CF_ii = (inv(SigmaFF[[i]])^(1/2))%*%(t(F_PC[[i]])%*%F_PC[[i]])%*%(inv(SigmaFF[[i]])^(1/2))
    CF_jj = (inv(SigmaFF[[j]])^(1/2))%*%(t(F_PC[[j]])%*%F_PC[[j]])%*%(inv(SigmaFF[[j]])^(1/2))
    
    
    RV[i,j] = tr(CF_ij%*%CF_ji)/sqrt(tr(CF_ii%*%CF_ii)*tr(CF_jj%*%CF_jj))  
    RV[j,i] = RV[i,j]

    
    
    
    Beta_i <- (tr(SigmaFF[[i]]))^2/tr(SigmaFF[[i]]%*%SigmaFF[[i]]) #This is just going to be 1 for any 1x1 matrix
    Beta_j <- (tr(SigmaFF[[j]]))^2/tr(SigmaFF[[j]]%*%SigmaFF[[j]]) #Would be trace which would be diff't for other dims
    
    
    ERV <- sqrt(Beta_i * Beta_j) / (t - 1)
    
    
    
    tau_isub <- sum(diag((t(F_PC[[i]])%*%F_PC[[i]]))^2)/(tr((t(F_PC[[i]]) %*% F_PC[[i]])%*%(t(F_PC[[i]]) %*% F_PC[[i]])))
    tau_jsub <- sum(diag((t(F_PC[[j]])%*%F_PC[[j]]))^2)/(tr((t(F_PC[[j]]) %*% F_PC[[j]])%*%(t(F_PC[[j]]) %*% F_PC[[j]])))
    
    
    tau_i <- ((t - 1)/((t - 3)*(t-1 - Beta_i))) * (t*(t+1)*tau_isub - (t-1)*(Beta_i + 2))
    tau_j <- ((t - 1)/((t - 3)*(t-1 - Beta_j))) * (t*(t+1)*tau_jsub - (t-1)*(Beta_j + 2))
  
    
    VRV <- ((2*(t-1-Beta_i)*(t-1-Beta_j))/((t+1)*(t-1)^2 * (t-2)))* (1 + ((t-3)/(2*t*(t-1)))*tau_i*tau_j)
    t_RV[i,j] <- (RV[i,j] - ERV)/sqrt(VRV)
    
  }
}




alpha = 0.05
D = 133*133 

p = (1 - (alpha)/(D)) 
p1 = 1 - (alpha)/(133*132/2)
threshold <- qnorm(p)
threshold1 <- qnorm(p1)

Connectivity_50 <- abs(t_RV) > threshold

Connectivity_50_1 <- abs(t_RV) > threshold1

Connectivity_3_Ting <- abs(t_RV) > 3

diag(Connectivity_50) <- NA
diag(Connectivity_50_1) <- NA
diag(Connectivity_3_Ting) <- NA

which(Connectivity_50, arr.ind = T) 

which(Connectivity_50_1, arr.ind = T)
which(Connectivity_3_Ting, arr.ind = T) 


significant <- which(Connectivity_50, arr.ind = T)
sig_rows <- significant[,1]
sig_cols <- significant[,2]
tRV_sig_val_rest <- c()
RV_sig_val_rest <- c()

for(i in 1:156){
  tRV_sig_val_rest[i] <- t_RV[sig_rows[i],sig_cols[i]]
  RV_sig_val_rest[i] <- RV[sig_rows[i], sig_cols[i]]
}


significant_thresh3 <- which(Connectivity_3_Ting, arr.ind = T)
sig_rows3 <- significant_thresh3[,1]
sig_cols3 <- significant_thresh3[,2]
tRV_sig_val_rest3 <- c()
RV_sig_val_rest3 <- c()

for(i in 1:25){
  tRV_sig_val_rest3[i] <- t_RV[sig_rows3[i],sig_cols3[i]]
  RV_sig_val_rest3[i] <- RV[sig_rows3[i], sig_cols3[i]]
}


RV_large <- ((RV > .08936768) &  (RV < 1))
RV_large_vec <- which(RV_large,arr.ind = TRUE)
large_rows <- RV_large_vec[,1]
large_cols <- RV_large_vec[,2]
RV_large_vals <- c()
tRV_large_vals <- c()


for(i in 1:338){
  tRV_large_vals[i] <- t_RV[large_rows[i],large_cols[i]]
  RV_large_vals[i] <- RV[large_rows[i], large_cols[i]]
}


require(gplots)

heatmap(RV, Rowv = NA, Colv = NA, revC = TRUE)
heatmap(abs(t_RV), Rowv = NA, Colv = NA, revC = TRUE)


RV[1,1] <- 1
rownames(RV) <- as.character(c(1:133))
colnames(RV) <- as.character(c(1:133))

for(i in 1:133){
  if(odd(as.numeric(colnames(RV)[i])) == FALSE){
    colnames(RV)[i] <- ""
  }
}


pheatmap(RV, treeheight_row = 0, treeheight_col = 0, show_rownames = F, main = "Resting State RV Connectivity",
         fontsize = 8, angle_col = "270", show_colnames = T, cluster_cols = F, cluster_rows = F,col = blue2green2red(100))

pheatmap(RV, treeheight_row = 0, treeheight_col = 0, show_rownames = F, main = "Resting State RV Connectivity",
         fontsize = 8, angle_col = "270", show_colnames = T, col = blue2green2red(100))



write.matrix(significant, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/SigROI_Rest.csv")
write.matrix(tRV_sig_val_rest, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/tRV_Rest.csv")
write.matrix(RV_sig_val_rest, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/RV_Rest.csv")


write.matrix(significant_thresh3, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/SigROI_Rest3.csv")
write.matrix(tRV_sig_val_rest3, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/tRV_Rest3.csv")
write.matrix(RV_sig_val_rest3, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/RV_Rest3.csv")

write.matrix(RV_large_vec, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/Large_RV.csv")
write.matrix(tRV_large_vals, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/tRV_Rest_Large.csv")
write.matrix(RV_large_vals, file = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/RV_Rest_Large.csv")
