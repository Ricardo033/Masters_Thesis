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

x_all <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/1x_all_milk.xlsx")

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
#  DOE sample selection n=60 #
#----------------------------#

PCAX=pca$x[,c(0:20)]

Z_1=optFederov(data = PCAX, nTrials = 60, center=TRUE,approximate=FALSE,criterion="D")
Z1=as.vector(Z_1$rows)

# Substracting 1 to adjust data to python
z_python1 = Z1-1
dput(z_python1)



#------------------------------#
#  DOE sample selection n=90   #
#------------------------------#

PCAX=pca$x[,c(0:20)]

Z_2=optFederov(data = PCAX, nTrials = 90, center=TRUE,approximate=FALSE,criterion="D")
Z2=as.vector(Z_2$rows)

# Substracting 1 to adjust data to python
z_python2 = Z2-1
dput(z_python2)


#------------------------------#
#  DOE sample selection n=120  #
#------------------------------#

PCAX=pca$x[,c(0:20)]

Z_3=optFederov(data = PCAX, nTrials = 120, center=TRUE,approximate=FALSE,criterion="D")
Z3=as.vector(Z_3$rows)

# Substracting 1 to adjust data to python
z_python3 = Z3-1
dput(z_python3)

#------------------------------#
#  DOE sample selection n=193  #
#------------------------------#

## DOE to calculate the Z matrix
PCAX=pca$x[,c(0:20)]

Z_0=optFederov(data = PCAX, nTrials = 193, center=TRUE,approximate=FALSE,criterion="D")
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
