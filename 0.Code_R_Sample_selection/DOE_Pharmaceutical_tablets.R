#===========#
# Library   #
#===========#

library(AlgDesign)
library(tidyverse)
library(dplyr)
library(factoextra)
library(readxl)

#==========================#
# Calibration Dataset      #
#==========================#

x_all <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/3x_all_pharma.xlsx")

#=======================#
# Principal Components  #
#=======================#

pca <- prcomp(x_all, scale = F)
fviz_eig(pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             geom = "point"
)

# Eigenvalues
eig.val <- get_eigenvalue(pca)
eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 


#===============================================================#

#----------------------------#
#  DOE sample selection n=49 #
#----------------------------#

PCAX=pca$x[,c(0:20)]

Z_0=optFederov(data = PCAX, nTrials = 49, center=TRUE,approximate=FALSE,criterion="D")
Z=as.vector(Z_0$rows)

# Substracting 1 to adjust data to python
z_python = Z-1
dput(z_python)


#------------------------#
#  DOE samples plot      #
#------------------------#


## Plot
ind.sup_DOE <-x_all[Z,]
ind.sup.coord_DOE <- predict(pca, newdata = ind.sup_DOE)

# Plot of active individuals
p_DOE <- fviz_pca_ind(pca, repel = TRUE, geom = "point")

# Add supplementary individuals
fviz_add(p_DOE, ind.sup.coord_DOE, color ="red",addlabel = F)
