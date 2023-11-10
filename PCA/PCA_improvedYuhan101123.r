#Load libraries for the script.

library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(gridExtra)
library(bio3d)
library(parallel)
library(fastcluster)

#Set the working directory for the analysis.

setwd("/home/ucbecla/Scratch")

# Function for PCA analysis
perform_pca <- function(trj, pdb, output_prefix) {
  # Read DCD trajectory and associated PDB file
  trj <- read.dcd(trj)
  pdb <- read.pdb(pdb)
  
  # Select carbon-alpha atoms
  ca.inds <- atom.select(pdb, elety = "CA")
  
  # Create a 3-dimensional matrix for coordinates
  xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)
  
  # Save matrices
  save(trj, file = paste0(output_prefix, "_TRJ.RData"))
  save(xyz, file = paste0(output_prefix, "_XYZ.RData"))
  
  # Perform PCA analysis
  pc <- pca.xyz(xyz[, ca.inds$xyz], mass = pdb)
  
  # Save PCA results
  save(pc, file = paste0(output_prefix, "_PCA.RData"))
  
  # Remove unnecessary objects
  rm(trj, xyz)
  
  return(pc)
}

# Example usage
pc <- perform_pca("54runsOneQuarter.dcd", "backbone_beOption.pdb", "run54_onequarter")

# Elbow plot for determining the number of clusters
k <- list()
for (i in 1:10) {
  k[[i]] <- kmeans(pc$z, i) 
}

betweenss_totss <- sapply(k, function(x) x$betweenss/x$totss)

# Plot Elbow plot
tiff(paste0(output_prefix, "_Elbow.tif"), width = 7, height = 7, units = "in", res = 1200)
plot(1:10, betweenss_totss, type = "b", ylab = "Between SS / Total SS", xlab = "Clusters (k)")
dev.off()

# Determine the number of clusters to use
num_clusters <- 4

# Do hierarchical clustering
hc <- hclust(dist(pc$z[, 1:2]), method = "average")
dend <- as.dendrogram(hc)

# Save clustering results
save(hc, file = paste0(output_prefix, "_HC.RData"))
save(dend, file = paste0(output_prefix, "_Dendrogram.RData"))

# Cut the dendrogram to get cluster assignments
grps <- cutree(hc, k = num_clusters)
write.table(grps, paste0(output_prefix, "_Clusters.txt"), sep = "\t")

# Plot PC plot with colored clusters
tiff(paste0(output_prefix, "_PC_Plot.tif"), width = 7, height = 7, units = "in", res = 1200)
plot(pc, col = c("black", "red", "blue", "green")[grps])
dev.off()

# Get midpoints of each cluster
get_mid <- function(z, clust) {
  mid_clust <- colMeans(z[grps == clust, 1:2])
  rel <- z[grps == clust, 1:2] - mid_clust
  frame <- which(sqrt(rel[, 1]^2 + rel[, 2]^2) == min(sqrt(rel[, 1]^2 + rel[, 2]^2)))
  frame <- which(sqrt(rel[, 1]^2 + rel[, 2]^2) %in% min(sqrt(rel[, 1]^2 + rel[, 2]^2)))[1]
  mid_rep <- z[grps == clust, 1:2][frame,]
  rep_frame <- which(z[, 1:2] == mid_rep)[1]
  rep_frame <- which(z[, 1:2] %in% mid_rep)[1]
  return(rep_frame)
}

# Get midpoints for each cluster
midpoints <- lapply(1:num_clusters, function(c) get_mid(pc$z, c))

# Print midpoints
for (i in seq_along(midpoints)) {
  print(paste("Midpoint for Cluster", i, "is at Frame", midpoints[[i]]))
}