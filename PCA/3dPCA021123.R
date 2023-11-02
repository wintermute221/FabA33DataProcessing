"""
Script Name: 3D PCA Script
Author: Yuhan Wang
Creation Date: 2nd Nov, 2023

Description:
This script aims to generate 3D plots with the input of the atom coordinates from frames generated from MD simulations

- This script is intended for educational purposes only.
- Use at your own risk.

"""


# Install library poorly

# install.packages("plotly")

# import library

library(plotly)
library("scatterplot3d")
library(viridis)
library(ggplot2)
library(grid)
library(plyr)
require(gridExtra)
library(ggsci)
library(extrafont)
#library(gridBase)
library(bio3d)

# read PDB file and trajectory 

trj <- read.dcd("4_5_338_50_r5_md_0_1_backbone_50ns.dcd")
pdb <- read.pdb("backbone_beOption.pdb")

#Select only the carbon-alpha atoms, reduces computational cost.

ca.inds <- atom.select(pdb, elety = "CA")

# select the alpha carbon atoms and save their coordinates into a csv file

ca <- pdb$atom[ca.inds$atom,]

# get the corrdinates
ca_xyz <- ca[,c("x","y","z")]

# write into a csv file

#write.csv(ca_xyz, file= "AlphaCarbon_xyz_backbone_beOption.csv")

#Create a 3-dimensional matrix which stores the coordinates of every frame,
#And split among number of cores determined previously (can be removed if not necessary).
xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)

pc <- pca.xyz(xyz[, ca.inds$xyz], mass = pdb)


#Do the clustering using hierarchical clustering. "Average" method is quoted as the best for MD frames.

hc <- hclust(dist(pc$z[, 1:2]), method = "average")

# assign the cluster index into a numeric vector

grps <- cutree(hc, k=6)

# Visualize the dendrogram

plot(hc, main="Hierarchical Clustering Dendrogram", xlab="Index in original data", sub="", cex=0.6)

# write the pac matrix into a csv file
write.csv(pca, file= "4_5_338_50_r5_md_0_1_backbone_50ns.csv")

## plot a 2D figure of PCs
pca <- plot(pc, col = grps)

## 3D plots of PC

# Read data
data = read.csv("4_5_338_50_r5_md_0_1_backbone_50ns.csv")

# Assign colours for different groups
colors <- c("#000099", "#CC00CC", "#009900", "#FF9933", "#66FFFF", "#d466ff")
colors <- colors[(grps)]

# plot the 3D figure for PCA
scatterplot3d(data[,2:4], pch = 16, color=colors, angle = 45) # try different angles 135, 225, and 315


## interactive 3D plots
# Read the AlphaCarbon_xyz_backbone_beOption csv file to get the coordinate

csv_ac <- read.csv("AlphaCarbon_xyz_backbone_beOption.csv")

# plot the backbone protein
pdb_3d <- plot_ly(csv_ac, x = csv_ac[,2],y=csv_ac[,3],z=csv_ac[,4], color = "#BF382A")

# Adjust the size of the dots
pdb_3d <- pdb_3d %>% add_markers(marker = list(size = 3))

# plot pca

colors <- c("#000099", "#CC00CC", "#009900", "#FF9933", "#66FFFF", "#d466ff")
colors <- colors[(grps)]

# Combine the grps vector into a dataframe
df <- data.frame(data, grps)

onedcd_3d <- plot_ly(df, x = df[,2],y=df[,3],z=df[,4],color = ~grps)

# Adjust the size of the dots
onedcd_3d <- onedcd_3d %>% add_markers(marker = list(size = 3))



