```
# This script is for i) firstly applying PCA on simulation trajectories so as to grasp the dominant motions from the simulations; ii) then based on the PCA results which assumably generate an albow plot including the variance and the eigenvectors of all the eigenvalues calculated, a clustering analysis then could be done
# This script was used on a total of 54 runs * 6 replicas from the Fab A33 simulations finished between April - June 2022, and Oct - Dec 2022 to further analyse the data so as to have an idea of the correlation between certain residues on the protein with aggregation kinetics

```


#Load libraries for the script.

library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(gridExtra)
library(extrafont)
library(bio3d)
library(parallel)
library(fastcluster)
library(factoextra)


#Read in your DCD trajectory and associated PDB file (must have the same number of atoms).

pdb <- read.pdb("backbone_beOption.pdb")
trj <- read.dcd("A_Half_Aligned.dcd")

#Select only the carbon-alpha atoms, reduces computational cost.

ca.inds <- atom.select(pdb, elety = "CA")

#Create a 3-dimensional matrix which stores the coordinates of every frame,
#And split among number of cores determined previously (can be removed if not necessary).
xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)

#Save this matrix and trajectory as R Object files, if you want to load them up later,
#Saves you having to run the "fit.xyz" command for every new repeat of the same trajectory.
save(trj, file = "All_HalfFrame_TRJ.RData")
save(xyz, file = "All_HalfFrame_XYZ.RData")

#Remove the "trj" object from this session, it is no longer used.
#This saves RAM for the rest of the analysis.

rm(trj)

#Perform the PCA analysis.

#load("All_HalfFrame_XYZ.RData")

pc <- pca.xyz(xyz[, ca.inds$xyz], mass = pdb)

#Remove the matrix object from memory to save RAM.

rm(xyz)

#Save the PCA object as an R Object file.

save(pc, file = "All_HalfFrame_PCA.RData")

#Get the Elbow plot of your PCA. This helps determine the number of clusters to use.

#load("All_HalfFrame_PCA.RData")

k <- list()
for(i in 1:10){      #Number of clusters e.g. 1 to 10
  k[[i]] <- kmeans(pc$z, i) 
}

betweenss_totss <- list()
for(i in 1:10){
  betweenss_totss[[i]] <- k[[i]]$betweenss/k[[i]]$totss
}

tiff("Elbow_LowpH.tif", width=7, height=7, units="in", res=1200)
plot(1:10, betweenss_totss, type = "b", ylab = "Between SS / Total SS", xlab = "Clusters (k)")

#Comment out everything after here for first run; then use determined number of clusters in this 2nd half.

#Do the clustering using hierarchical clustering. "Average" method is quoted as the best for MD frames.

hc <- hclust(dist(pc$z[, 1:2]), method = "average")

#Change type of object hc is. This is useful for plotting the dendrogram later.
#Save the dendrogram object as an R Object file to load it later.

dend = as.dendrogram(hc)
save(hc, file = "A_LowpH_350_HC.RData")
save(dend, file = "A_LowpH_350_Dendrogram.RData")

#Set the number of clusters (k) to use as determined by elbow plot.
#Save the clustering as a text file. This gives you frame number vs cluster number in 2 columns.
#Then remove unnecessary objects from memory.

#load("A_LowpH_350_HC.RData")

grps <- cutree(hc, k = 4)
write.table(grps, "4_Clusters_LowpH.txt", sep="\t")
rm(hc)
rm(dend)

#Get the PC plot coloured with the clusters identified above. PC1vsPC, PC2vsPC3, PC1vsPC3 and Scree plot

tiff("PC_Plot_4Clusters_LowpH.tif", width=7, height=7, units="in", res=1200)
plot(pc, col=c("black", "red", "blue", "green")[grps])

#Get the mid-point of each cluster. This is the closest structure to every other in that cluster.

get_mid <- function(z, clust){
    mid_clust <- colMeans(z[grps == clust,1:2])
    rel <- z[grps == clust,1:2] - mid_clust
    frame <- which(sqrt(rel[,1]**2+rel[,1]**2) ==  min(sqrt(rel[,1]**2+rel[,1]**2)))
    frame <- which(sqrt(rel[,1]**2+rel[,1]**2) %in%  min(sqrt(rel[,1]**2+rel[,1]**2)))[1]
    mid_rep <- z[grps == clust,1:2][frame,]
    rep_frame <- which(z[,1:2] == mid_rep)[1]
    rep_frame <- which(z[,1:2] %in% mid_rep)[1]
    return(rep_frame)
}

#Reduce or increase to the number of clusters you are using.

mid_c1 <- get_mid(pc$z,1)
print(mid_c1)
mid_c2 <- get_mid(pc$z,2)
print(mid_c2)
mid_c3 <- get_mid(pc$z,3)
print(mid_c3)
mid_c4 <- get_mid(pc$z,4)
print(mid_c4)
#mid_c5 <- get_mid(pc$z,5)
#print(mid_c5)
#mid_c6 <- get_mid(pc$z,6)
#print(mid_c6)
#mid_c7 <- get_mid(pc$z,7)
#print(mid_c7)

dev.off()
stopCluster(cl)