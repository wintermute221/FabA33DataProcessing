# load bio3d
library(bio3d)

#set directory as where the DCD files are
setwd('/Users/wintermute/Library/CloudStorage/GoogleDrive-wintermute.backup@gmail.com/My Drive/MacbookPro/GoogleDrive/UCL/PaulGroup/MolecularDynamics/Gromacs/FabA33/DataAnalysisAll200323/DCCM170423/ConvertDCD140623/54DCD140623')

# read pdb
pdb <- read.pdb("Fab.pdb")

# define working directory
filepath <- getwd()

# get all the dcd files in the directory
files <- list.files(pattern = "\\.dcd$")

for (i in files) {

    #read dcd file
    trj <- read.dcd(i)

    #select alpha carbon
    inds <- atom.select(pdb, elety='CA')

   ## lsq fit of trj on pdb

    xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

    #Calculate the dynamic cross-correlation matrix.


    cij <- dccm(xyz[,inds$xyz])

    #save as png image in specific directory with 600*350 resolution
    #png(file = file.path(filepath, paste0(toString(i), ".png")))
    png(file = file.path(filepath, paste0(tools::file_path_sans_ext(i), ".png")))

    # a DCCM we want to save
    plot(cij)

    # a function call to save the file
    dev.off()

}

# improvement for future use: maybe I can edit the png line a bit so that I could output all the plots into a specific directory