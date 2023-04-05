rm(list = ls())
library(ciftiTools)
library(Rcmdr)
library(factoextra)
library(plot.matrix)

####################################### FUNCTIONS #############################
EVFunction <- function(onset, duration, TR){
  temp_indices = vector(mode = "list", length = length(onset))
  end_onset <- c()
  
  for(j in 1:length(onset)){
    end_onset[j] <- onset[j] + duration[j]
    
    shift_seq <- seq(from = TR, to = 214, by = TR) #214 is the run duration - 3min &34seconds = 214s
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

################################################################################
# ALREADY RUN THIS IN MOTOR ANALYSIS STUFF
# ciftiTools.setOption("wb_path", "C:/Users/Owner/Documents/Graduate Documents/Research/HCP")
filename = "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/100307/MNINonLinear/Results/tfMRI_MOTOR_LR/EVs"
# 
# task <- read_xifti("~/Graduate Documents/Research/HCP/Data/100307/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii",
#                   brainstructures =  "all")


t = 284
p = 5 #lh, rh, lf, rf, tongue, & cue - not sure whether to include cue...probably?
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


#Beta_HRF_1 = matrix(rep(0, times = t*p), nrow = t, ncol = p)
#Beta_dHRF_1 = matrix(rep(0, times = t*p), nrow = t, ncol = p)

for(i in 1:p){
  hrf_temp = convolve(hrf(dt_aux[1:284]), BoxCar_Func[,i])   #Ask Michelle about these different lengths??? 
  dhrf_temp = convolve(hrf_deriv(dt_aux[1:284]), BoxCar_Func[,i])
  
  Beta_HRF[,i] = hrf_temp[1:284]
  Beta_dHRF[,i] = dhrf_temp[1:284]
  # Beta_HRF_1[,i] = hrf_temp[1:284]
  # Beta_dHRF_1[,i] = dhrf_temp[1:284]
}

DesignMatrix = matrix(,nrow = t, ncol = 2*p)

for(j in 1:p){
  DesignMatrix[, (j-1)*2 + 1] = Beta_HRF[,j]
  DesignMatrix[, (j-1)*2 + 2] = Beta_dHRF[,j]
}

write.table(DesignMatrix, "C:/Users/Owner/Documents/Graduate Documents/Research/HCP/Data/DesignMat.csv", col.names = FALSE, row.names = FALSE)
