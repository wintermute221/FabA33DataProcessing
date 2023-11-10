#Load libraries for the script.

library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(gridExtra)
library(bio3d)
library(parallel)
library(fastcluster)



analyze_dcd_trajectory <- function(dcd_file, pdb_file, output_prefix) {
 
 
  # Set the working directory
  # setwd("/home/ucbecla/Scratch")

  # # Split the analysis over several cores
  # setup.ncore(ncore)
  # cl <- makeForkCluster(nnodes = ncore)

  # Read DCD trajectory and associated PDB file
  trj <- read.dcd(dcd_file)
  pdb <- read.pdb(pdb_file)

  # Select only the carbon-alpha atoms
  ca.inds <- atom.select(pdb, elety = "CA")

  # Create a 3-dimensional matrix for coordinates
  xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)

  # # Save the trajectory as an R Object file
  # save(trj, file = paste0(output_prefix, "_TRJ.RData"))

  # Remove the "trj" object
  rm(trj)

  # # Save the XYZ matrix as an R Object file
  # save(xyz, file = paste0(output_prefix, "_XYZ.RData"))

  # Calculate the dynamic cross-correlation matrix
  cij <- dccm(xyz[, ca.inds$xyz])

  # Write cross-correlation data as a table
  write.table(cij, paste0(output_prefix, "_Cross_Corr.txt"), sep = "\t")

  # plot
  tiff(paste0(output_prefix, ".tif"), width=7, height=7, units="in", res=1200)
  plot(cij)

  # Close the plotting device
  dev.off()

  # Remove specific objects to free up memory
  rm(pdb, ca.inds, xyz, cij)

  # Garbage collection to free up additional memory
  gc()

  # # Stop the parallel cluster
  # stopCluster(cl)
}

# Example usage
# analyze_dcd_trajectory("Trajectory_File.dcd", "PDB_File.pdb", "E7_cpH")

# Example usage in a for loop to loop over all the dcd files 
# Set the directory path
directory_path <- "/your/directory/path"

# Get a list of all files in the directory
all_files <- list.files(directory_path)

# Filter files with ".dcd" extension
dcd_files <- grep("\\.dcd$", all_files, value = TRUE)

# Iterate over each dcd file
for (file in dcd_files) {
  # Extract the file prefix (remove extension)
  file_prefix <- tools::file_path_sans_ext(file)
  
  analyze_dcd_trajectory(file, "Fabs.pdb", file_prefix)

}
