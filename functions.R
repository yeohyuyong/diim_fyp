library(openxlsx)


download_data = function(){
  # download 2019 IOT and extract total output and final demand column
  data = read.xlsx("/Users/yeohyy/Library/CloudStorage/GoogleDrive-yeoh.yuyong@gmail.com/My Drive/School/NTU/FYP/Data/mrsd_70_Res_UR_ind.xlsx", sheet = "2019IOTable")
  iot2019 = data[2:(2+15-1), 3:(3+15-1)]
  x = data[2:(2+15-1), ncol(data)]
  c = data[2:(2+15-1), ncol(data)-1]
  x = as.matrix(sapply(x, as.numeric))
  c = as.matrix(sapply(c, as.numeric))
  
  # extract technical coefficient matrix
  A = read.xlsx("/Users/yeohyy/Library/CloudStorage/GoogleDrive-yeoh.yuyong@gmail.com/My Drive/School/NTU/FYP/Data/mrsd_70_Res_UR_ind.xlsx", sheet = "A")
  A = A[2:(2+15-1), 3:(3+15-1)]
  A = as.matrix(sapply(A, as.numeric))
  
  # extract sector initial inoperability
  q0 = read.xlsx("/Users/yeohyy/Library/CloudStorage/GoogleDrive-yeoh.yuyong@gmail.com/My Drive/School/NTU/FYP/Data/mrsd_70_Res_UR_ind.xlsx", sheet = "sector inoperability",colNames = FALSE)
  q0 = q0[1:(1+15-1),2]
  
  c_star = read.xlsx("/Users/yeohyy/Library/CloudStorage/GoogleDrive-yeoh.yuyong@gmail.com/My Drive/School/NTU/FYP/Data/mrsd_70_Res_UR_ind.xlsx", sheet = "c_star")
  c_star = c_star[1:(1+15-1),7]
  
  
  
  A_star = solve(diag(as.vector(x)))%*%A%*%diag(as.vector(x))
  
 
  
  return(list(
     x = x,
     c = c,
     A = A,
     q0 = q0,
     c_star = c_star,
     A_star = A_star
  ))
}


# construct_H_matrix <- function(A) {
#   # Calculate Leontief inverse
#   I_minus_A = diag(nrow(A)) - A
#   L = solve(I_minus_A)  # Leontief inverse: (I - A)^{-1}
# 
#   # Compute H = A %*% L
#   H = A %*% L
# 
#   return(H)
# }
# 
# 
# perform_pca = function(H_matrix) {
#   # Conduct PCA using eigendecomposition
#   eigen_result <- eigen(H_matrix)
#   eigenvalues <- eigen_result$values
#   eigenvectors <- eigen_result$vectors
#   
#   # Calculate variance explained
#   variance_explained <- eigenvalues / sum(eigenvalues)
#   cumulative_variance <- cumsum(variance_explained)
#   
#   return(list(
#     eigenvalues = Re(eigenvalues),
#     eigenvectors = Re(eigenvectors),
#     variance_explained = Re(variance_explained),
#     cumulative_variance = Re(cumulative_variance)
#   ))
# }
# 
# 
# 
# 
# 
# pca_result <- perform_pca(H)
# 
# # Extract first and second principal components
# PC1 <- pca_result$eigenvectors[, 1]
# PC2 <- pca_result$eigenvectors[, 2]
# 
# # Create data frame for visualization
# pca_data <- data.frame(
#   Sector_ID = 1:nrow(A),
#   PC1 = PC1,
#   PC2 = PC2,
#   Eigenvalue_1 = pca_result$eigenvalues[1],
#   Eigenvalue_2 = pca_result$eigenvalues[2]
# )
# 
# top_3_sectors = order(pca_data[,2], decreasing = TRUE)[1:3]




DIIM = function(q0, A_star,c_star,x,lockdown_duration, total_duration,key_sectors=NULL,risk_management=1) {
  
  a_ii = diag(A_star)
  T = total_duration
  
  num_sectors = length(q0)
  inoperability_evolution = matrix(NA,nrow = num_sectors,ncol=total_duration)
  
  
  # if key_sectors is not NULL, reduce those sectors' initial q0 by 10%
  if (!is.null(key_sectors)) {
    q0[key_sectors] = q0[key_sectors] * 0.9
  }
  
  # we assume after 2 years the economic activity return to 99% of pre lock down level
  qT = q0 * 1/100
  
  k = log(q0/qT)/(T*(1-a_ii))
  K <- diag(as.vector(k))

  inoperability_evolution[,1] = q0
  
  for (t in 2:total_duration) {
    if (t <= lockdown_duration) {
      inoperability_evolution[, t] = inoperability_evolution[, t - 1] + K %*% (A_star %*% inoperability_evolution[, t - 1] + c_star - inoperability_evolution[, t - 1]) * risk_management
    }
    else {
      # after lockdown c* become 0 as no more external shock
      inoperability_evolution[, t] = inoperability_evolution[, t - 1] + K %*% (A_star %*% inoperability_evolution[, t -1] - inoperability_evolution[, t - 1]) * risk_management
    }
  }
  
  # vector x is the planned YEARLY output of all sectors
  # convert to planned DAILY output of all sectors
  # note that output is in millions
  x_daily = x/366 #366 days in 2020
  EL_evolution = t(apply(inoperability_evolution, 1, cumsum)) * as.vector(x_daily)
  EL_end = EL_evolution[,ncol(EL_evolution)]
  total_economic_loss = sum(EL_evolution[,ncol(EL_evolution)])
  
  return(
    list(
      inoperability_evolution = inoperability_evolution,
      EL_evolution = EL_evolution,
      EL_end = EL_end,
      total_economic_loss = total_economic_loss
    )
  )
}


simulation_ml_vs_diim = function(q0, A_star,c_star,x,lockdown_duration, total_duration){
  
  # run model to calculate the economic loss (without intervention)
  model = DIIM(q0, A_star,c_star,x,lockdown_duration, total_duration)
  EL_evolution = model$EL_evolution
  
  # obtain the sectors with top economic loss
  max_econ_loss <- apply(EL_evolution, 1, max)
  sorted_indices <- order(max_econ_loss, decreasing = TRUE)
  top_econ_loss_3 = sorted_indices[1:5]
  
  # rerun the model with intervention for the top economic sectors
  
  
  model_diim = DIIM(q0, A_star,c_star,x,lockdown_duration, total_duration,key_sectors = top_econ_loss_3)
  
  ml_key_sectors = c(3,2,5,8,10) # obtained from PCA
  model_ml = DIIM(q0, A_star,c_star,x,lockdown_duration, total_duration,key_sectors = ml_key_sectors)
  
  model_tot_econ_loss = model$total_economic_loss
  model_diim_tot_econ_loss = model_diim$total_economic_loss
  model_ml_tot_econ_loss = model_ml$total_economic_loss
  
  return(list(
    lockdown_duration = lockdown_duration,
    total_duration = total_duration,
    model_tot_econ_loss = model_tot_econ_loss,
    model_diim_tot_econ_loss = model_diim_tot_econ_loss,
    model_ml_tot_econ_loss = model_ml_tot_econ_loss
  ))
}



