library(openxlsx)
library(ggplot2)
library(igraph)
library(factoextra)
library(FactoMineR)
library(dendextend)
library(tidyverse)
library('corrr')
library(ggcorrplot)
library("FactoMineR")

data = read.xlsx("/Users/yeohyy/Library/CloudStorage/GoogleDrive-yeoh.yuyong@gmail.com/My Drive/School/NTU/FYP/Data/mrsd_70_Res_UR_ind.xlsx", sheet = "2019IOTable")
iot2019 = data[2:(2+15-1), 3:(3+15-1)]
x = data[2:(2+15-1), ncol(data)]

A = read.xlsx("/Users/yeohyy/Library/CloudStorage/GoogleDrive-yeoh.yuyong@gmail.com/My Drive/School/NTU/FYP/Data/mrsd_70_Res_UR_ind.xlsx", sheet = "A")
A = A[2:(2+15-1), 3:(3+15-1)]
A = as.matrix(sapply(A, as.numeric))
n_sectors = nrow(A)

I_minus_A = diag(nrow(A)) - A
L = solve(I_minus_A)


calculate_linkages <- function(L_matrix, x_vector = NULL) {
  # Leontief inverse should be passed as L
  
  # Backward linkages (column sums normalized)
  BL <- colSums(L_matrix)
  BL_normalized <- BL / mean(BL)
  
  # Forward linkages (row sums normalized)
  FL <- rowSums(L_matrix)
  FL_normalized <- FL / mean(FL)
  
  # Total linkages (average)
  total_linkages <- (BL_normalized + FL_normalized) / 2
  
  return(data.frame(
    Sector_ID = 1:length(BL),
    BL = BL,
    BL_normalized = BL_normalized,
    FL = FL,
    FL_normalized = FL_normalized,
    Total_Linkage = total_linkages
  ))
}

linkages <- calculate_linkages(L)

cat("\nRasmussen-Hirschman Linkages:\n")
print(linkages %>% arrange(desc(Total_Linkage)))



construct_H_matrix <- function(A) {
  # Calculate Leontief inverse
  I_minus_A = diag(nrow(A)) - A
  L = solve(I_minus_A)  # Leontief inverse: (I - A)^{-1}
  
  # Compute H = A %*% L
  H = A %*% L
  
  return(H)
}

H = construct_H_matrix(A)


perform_pca <- function(H_matrix) {
  # Conduct PCA using eigendecomposition
  eigen_result <- eigen(H_matrix)
  eigenvalues <- eigen_result$values
  eigenvectors <- eigen_result$vectors
  
  # Calculate variance explained
  variance_explained <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(variance_explained)
  
  return(list(
    eigenvalues = Re(eigenvalues),
    eigenvectors = Re(eigenvectors),
    variance_explained = Re(variance_explained),
    cumulative_variance = Re(cumulative_variance)
  ))
}

pca_result <- perform_pca(H)

# Extract first and second principal components
PC1 <- pca_result$eigenvectors[, 1]
PC2 <- pca_result$eigenvectors[, 2]

# Create data frame for visualization
pca_data <- data.frame(
  Sector_ID = 1:nrow(A),
  PC1 = PC1,
  PC2 = PC2,
  Eigenvalue_1 = pca_result$eigenvalues[1],
  Eigenvalue_2 = pca_result$eigenvalues[2]
)



plot_scree <- function(pca_result, sector_names = NULL) {
    n_components <- min(10, length(pca_result$variance_explained))  
    scree_data <- data.frame(
      Component = 1:n_components,
      Variance_Explained = pca_result$variance_explained[1:n_components] * 100,
      Cumulative = pca_result$cumulative_variance[1:n_components] * 100
  )
  
  p <- ggplot(scree_data, aes(x = Component)) +
    geom_line(aes(y = Variance_Explained, color = "Individual"), size = 1) +
    geom_point(aes(y = Variance_Explained, color = "Individual"), size = 3) +
    geom_line(aes(y = Cumulative, color = "Cumulative"), size = 1) +
    geom_point(aes(y = Cumulative, color = "Cumulative"), size = 3) +
    geom_hline(
      yintercept = 70,
      linetype = "dashed",
      color = "gray",
      size = 0.5
    ) +
    geom_hline(
      yintercept = 80,
      linetype = "dashed",
      color = "gray",
      size = 0.5
    ) +
    labs(
      title = "Scree Plot: Variance Explained by Principal Components",
      x = "Principal Component",
      y = "Variance Explained (%)",
      color = "Type"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p)
  
  return(scree_data)
}

scree_data <- plot_scree(pca_result)

cat("\nVariance Explained by Principal Components:\n")
print(
  data.frame(
    PC = 1:5,
    Variance = pca_result$variance_explained[1:5] * 100,
    Cumulative = pca_result$cumulative_variance[1:5] * 100
  )
)



plot_factor_map <- function(pca_data, sector_labels = NULL) {
  if (is.null(sector_labels)) {
    sector_labels <- paste("S", pca_data$Sector_ID, sep = "")
  }
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 4, color = "steelblue", alpha = 0.6) +
    geom_text(aes(label = Sector_ID), size = 3, hjust = -0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray", alpha = 0.5) +
    labs(
      title = "Factor Map: First Two Principal Components",
      x = paste("PC1 (", round(pca_result$variance_explained[1]*100, 1), "% variance)"),
      y = paste("PC2 (", round(pca_result$variance_explained[2]*100, 1), "% variance)")
    ) +
    theme_minimal() +
    coord_equal()
  
  print(p)
  
  return(p)
}

factor_map <- plot_factor_map(pca_data)





pca_matrix = pca_data[,2:3]
str(pca_matrix)

set.seed(100)

## function to compute total within-cluster sum of square 
wcss <- function(k) {
  kmeans(pca_matrix, k, nstart = 10 )$tot.withinss
}

k.values <- 1:10
# apply wcss to all k values
wcss_k<-sapply(k.values, wcss)
plot(k.values, wcss_k,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# or plot the graph using gglop
library(ggplot2)
#create a dataframe for ggplot
elbow<-data.frame(k.values, wcss_k)
ggplot(elbow, aes(x = k.values, y = wcss_k)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 30, by = 1))
#select the elbow point k=3
set.seed(100)
k3 <- kmeans(pca_matrix, centers = 3, nstart = 10)
str(k3)
k3 

## Hierarchical clustering using Complete Linkage
hc1 <- hclust(dist(pca_matrix), method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang =-.1)
#draw rectangles with different colors for 7 clusters respectively in the dendrogram, 
#where argument border specifies the colors of the rectangles, default value of border=1 for black.
rect.hclust(hc1, k=3, border = 2:8)



