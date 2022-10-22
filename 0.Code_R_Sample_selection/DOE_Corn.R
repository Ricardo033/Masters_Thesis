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

x_all <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/2x_all_corn.xlsx")

#-----------------------------------------------------------------------------#

#=======================#
# Principal Components  #
#=======================#

pca <- prcomp(x_all, scale = FALSE)
fviz_eig(pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             geom="point",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             
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


#------------------------#
#  DOE sample selection  #
#------------------------#

## DOE to calculate the Z matrix
PCAX=pca$x[,c(0:20)]

Z_0=optFederov(data = PCAX, nTrials = 37, center=TRUE,approximate=FALSE,criterion="D")
Z=as.vector(Z_0$rows)

# Substracting 1 to adjust data to python
z_python1 = Z-1
dput(z_python1)

#============================================#
#  Plot of selected samples in PCA space     #
#============================================#

#------------------------#
#  Ks sample selection   #
#------------------------#

ks=c( 0,  1,  2,  3,  4,  7,  8, 10, 12, 15, 16, 18, 20, 21, 22, 24, 26,
      27, 30, 31, 32, 33, 35, 36, 39, 40, 41, 42, 43, 44, 45, 47, 48, 49,
      50, 51, 55)

# As R deos not start a df from 0 we sum 1 to each position
ks = ks+1

ind.sup_ks <-x_all[ks,]
ind.sup.coord_ks <- predict(pca, newdata = ind.sup_ks)

# Plot of active individuals
p_ks <- fviz_pca_ind(pca, repel = TRUE, geom = "point")
# Add supplementary individuals
fviz_add(p_ks, ind.sup.coord_ks, color ="blue",addlabel = F)


#------------------------#
#  dup sample selection  #
#------------------------#

dup=c( 3,  4,  5,  6,  7,  8,  9, 11, 14, 15, 16, 17, 18, 19, 20, 23, 25,
       26, 27, 28, 29, 30, 31, 33, 37, 38, 40, 41, 43, 44, 45, 46, 47, 48,
       50, 53, 55)

# As R deos not start a df from 0 we sum 1 to each position
dup = dup +1

ind.sup_dup <-x_all[dup,] 
ind.sup.coord_dup <- predict(pca, newdata = ind.sup_dup)

# Plot of active individuals
p_dup <- fviz_pca_ind(pca, repel = TRUE, geom = "point")

# Add supplementary individuals
fviz_add(p_dup, ind.sup.coord_dup, color ="orange",addlabel = F )

#------------------------#
#  puch sample selection  #
#------------------------#

puch=c( 0,  1,  2,  3,  4,  6,  7,  8, 10, 12, 15, 16, 18, 20, 21, 22, 23,
        24, 26, 27, 30, 32, 33, 35, 36, 39, 40, 41, 42, 43, 45, 47, 48, 49,
        50, 51, 54)


puch = puch +1

ind.sup_puch <-x_all[puch,] 
ind.sup.coord_puch <- predict(pca, newdata = ind.sup_puch)

# Plot of active individuals
p_puch <- fviz_pca_ind(pca, repel = TRUE, geom = "point")

# Add supplementary individuals
fviz_add(p_puch, ind.sup.coord_puch, color ="green",addlabel = F)

#------------------------#
#  clus sample selection  #
#------------------------#


clus=c( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        18, 20, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 38, 39, 40, 41,
        43, 45, 48)


clus = clus +1

ind.sup_clus <-x_all[clus,] 
ind.sup.coord_clus <- predict(pca, newdata = ind.sup_clus)

# Plot of active individuals
p_clus <- fviz_pca_ind(pca, repel = TRUE,geom = "point")

# Add supplementary individuals
fviz_add(p_clus, ind.sup.coord_clus, color ="red",addlabel = F)


#------------------------#
#  DOE sample selection  #
#------------------------#

## Plot
ind.sup_DOE <-x_all[Z,]
ind.sup.coord_DOE <- predict(pca, newdata = ind.sup_DOE)

# Plot of active individuals
p_DOE <- fviz_pca_ind(pca, repel = TRUE,geom = "point")

# Add supplementary individuals
fviz_add(p_DOE, ind.sup.coord_DOE, color ="brown",addlabel = F)

#-----------------------------#
#  diuss_max sample selection  #
#-----------------------------#

diuss_max=c( 0,  2,  3,  4,  7,  9, 10, 13, 16, 17, 18, 19, 20, 24, 25, 26, 27,
             28, 30, 31, 33, 34, 36, 37, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51,
             53, 54, 55)


diuss_max = diuss_max +1

ind.sup_diuss_max <-x_all[diuss_max,] 
ind.sup.coord_diuss_max <- predict(pca, newdata = ind.sup_diuss_max)

# Plot of active individuals
p_diuss_max <- fviz_pca_ind(pca, repel = TRUE, geom = "point")

# Add supplementary individuals
fviz_add(p_diuss_max, ind.sup.coord_diuss_max, color ="purple",addlabel = F)

#-----------------------------#
#  diuss_sum sample selection  #
#-----------------------------#

diuss_sum=c( 0,  2,  3,  4,  7,  9, 10, 11, 13, 14, 15, 17, 19, 20, 21, 23, 24,
             25, 26, 28, 29, 31, 32, 34, 36, 37, 40, 42, 43, 44, 46, 47, 50, 51,
             52, 53, 54)


diuss_sum = diuss_sum +1

ind.sup_diuss_sum <-x_all[diuss_sum,] 
ind.sup.coord_diuss_sum <- predict(pca, newdata = ind.sup_diuss_sum)

# Plot of active individuals
p_diuss_sum <- fviz_pca_ind(pca, repel = TRUE,geom = "point")

# Add supplementary individuals
fviz_add(p_diuss_sum, ind.sup.coord_diuss_sum, color ="darkgrey",addlabel = F)



