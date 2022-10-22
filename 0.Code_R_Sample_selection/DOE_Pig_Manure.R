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

x_all <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/0x_all_pig.xlsx")

#=======================#
# Principal Components  #
#=======================#
pca <- prcomp(x_all, scale = FALSE)


fviz_eig(pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping,
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


#------------------------#
#  DOE sample selection  #
#------------------------#

## DOE to calculate the Z matrix
PCAX=pca$x[,c(0:20)]

#Comparing loadings and scores with python's results 

scores = data.frame(pca$x)
loadings = data.frame(pca$rotation)


Z_0=optFederov(data = PCAX, nTrials = 133, center=TRUE,approximate=FALSE,criterion="D")
Z=as.vector(Z_0$rows)

# Substracting 1 to adjust data to python
z_python = Z-1
dput(z_python)

#------------------------#
#  Ks sample selection   #
#------------------------#

ks=c( 1,   4,   9,  10,  11,  14,  15,  18,  19,  20,  21,  25,  26,
      27,  28,  29,  30,  32,  35,  40,  41,  45,  46,  47,  48,  51,
      53,  54,  56,  59,  62,  63,  66,  67,  68,  69,  72,  75,  76,
      77,  79,  81,  82,  88,  89,  90,  93,  96,  97, 102, 104, 107,
      110, 111, 112, 113, 115, 118, 123, 124, 125, 126, 127, 128, 129,
      131, 133, 137, 138, 142, 143, 144, 145, 146, 147, 148, 154, 155,
      156, 160, 162, 164, 166, 168, 169, 171, 172, 175, 178, 180, 181,
      182, 184, 187, 189, 190, 191, 193, 196, 198, 200, 203, 204, 206,
      207, 208, 211, 212, 213, 214, 215, 216, 221, 223, 224, 226, 227,
      228, 229, 230, 232, 236, 237, 238, 239, 242, 243, 245, 247, 250,
      252, 253, 254)

# As R deos not start a df from 0 we sum 1 to each position
ks = ks+1

ind.sup_ks <-x_all[ks,]
ind.sup.coord_ks <- predict(pca, newdata = ind.sup_ks)

# Plot of active individuals
p_ks <- fviz_pca_ind(pca, repel = TRUE, geom = "point")
# Add supplementary individuals
fviz_add(p_ks, ind.sup.coord_ks, color ="cyan3",addlabel = F)


#------------------------#
#  dup sample selection  #
#------------------------#

dup=c( 0,   1,   2,   3,   8,  10,  12,  13,  14,  15,  16,  18,  19,
       22,  23,  24,  29,  30,  33,  35,  39,  40,  41,  43,  47,  49,
       51,  52,  53,  54,  56,  59,  60,  61,  62,  63,  65,  66,  67,
       68,  71,  72,  74,  78,  80,  81,  83,  84,  85,  87,  90,  92,
       93,  97,  99, 105, 106, 108, 109, 112, 113, 115, 116, 118, 122,
       123, 124, 127, 128, 129, 130, 134, 135, 136, 137, 138, 140, 141,
       143, 144, 145, 147, 148, 149, 150, 151, 152, 154, 156, 157, 158,
       159, 162, 163, 164, 165, 168, 169, 170, 172, 173, 175, 177, 180,
       183, 184, 187, 190, 191, 193, 195, 196, 197, 200, 206, 207, 211,
       215, 216, 220, 222, 228, 229, 232, 235, 236, 239, 240, 248, 251,
       252, 253, 255)

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

puch=c( 1,   8,   9,  10,  11,  12,  14,  15,  17,  19,  20,  21,  25,
        26,  27,  28,  30,  32,  35,  38,  40,  41,  45,  46,  47,  48,
        51,  53,  55,  56,  59,  63,  66,  67,  68,  69,  75,  76,  77,
        78,  79,  81,  82,  88,  89,  90,  91,  93,  96,  97, 102, 104,
        107, 108, 110, 111, 112, 113, 115, 118, 123, 124, 125, 127, 128,
        129, 131, 133, 134, 137, 138, 142, 143, 144, 145, 146, 147, 154,
        156, 160, 164, 166, 168, 169, 172, 175, 176, 177, 178, 179, 180,
        181, 182, 184, 187, 189, 191, 193, 196, 198, 200, 203, 204, 206,
        207, 208, 211, 212, 213, 214, 215, 216, 221, 223, 224, 226, 227,
        228, 229, 230, 232, 236, 237, 238, 239, 241, 242, 243, 245, 247,
        250, 252, 254)


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


clus=c( 1,   2,   6,   9,  11,  12,  13,  18,  19,  20,  22,  24,  25,
        26,  27,  28,  29,  30,  31,  32,  33,  35,  36,  37,  40,  42,
        43,  45,  46,  47,  48,  53,  55,  56,  59,  61,  63,  64,  65,
        67,  68,  70,  71,  72,  73,  74,  75,  76,  77,  79,  80,  81,
        82,  84,  85,  87,  88,  89,  91,  93,  97, 100, 101, 102, 105,
        107, 108, 110, 112, 120, 124, 125, 126, 127, 128, 129, 133, 136,
        139, 140, 143, 145, 147, 153, 156, 160, 162, 163, 166, 168, 170,
        171, 172, 175, 176, 178, 179, 180, 181, 182, 184, 185, 186, 187,
        188, 189, 191, 195, 196, 197, 198, 203, 206, 207, 209, 211, 214,
        215, 219, 221, 222, 227, 228, 229, 230, 231, 234, 236, 241, 242,
        247, 250, 251)


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

diuss_max=c( 1,   2,   3,   4,   7,  10,  11,  13,  14,  16,  17,  18,  22,
             25,  26,  27,  28,  29,  32,  33,  34,  36,  37,  38,  39,  47,
             50,  51,  52,  53,  54,  57,  59,  62,  63,  65,  66,  67,  69,
             70,  71,  73,  76,  78,  79,  80,  89,  90,  91,  93,  94,  96,
             97,  98, 100, 102, 103, 104, 105, 109, 111, 116, 118, 119, 120,
             121, 122, 124, 126, 127, 128, 130, 131, 132, 134, 135, 138, 140,
             142, 145, 149, 155, 156, 160, 161, 163, 164, 165, 166, 167, 171,
             173, 175, 176, 178, 180, 182, 183, 185, 186, 190, 191, 192, 198,
             201, 202, 206, 207, 208, 210, 211, 215, 216, 219, 220, 223, 225,
             229, 230, 231, 234, 235, 237, 238, 240, 241, 242, 245, 246, 248,
             249, 250, 254)


diuss_max = diuss_max +1

ind.sup_diuss_max <-x_all[diuss_max,] 
ind.sup.coord_diuss_max <- predict(pca, newdata = ind.sup_diuss_max)

# Plot of active individuals
p_diuss_max <- fviz_pca_ind(pca, repel = TRUE, geom = "point")

# Add supplementary individuals
fviz_add(p_diuss_max, ind.sup.coord_diuss_max, color ="deeppink1",addlabel = F)

#-----------------------------#
#  diuss_sum sample selection  #
#-----------------------------#

diuss_sum=c( 1,   3,   4,   5,   6,   7,   8,  10,  12,  13,  14,  16,  18,
             19,  22,  23,  25,  27,  29,  38,  41,  42,  43,  46,  50,  51,
             52,  53,  54,  56,  58,  59,  60,  64,  65,  66,  67,  68,  70,
             71,  73,  75,  76,  89,  90,  97,  98,  99, 100, 105, 107, 108,
             113, 114, 115, 116, 118, 119, 120, 121, 124, 126, 127, 129, 131,
             132, 136, 137, 138, 140, 141, 142, 145, 148, 149, 154, 155, 157,
             158, 161, 164, 166, 171, 173, 175, 178, 179, 180, 183, 184, 186,
             187, 188, 190, 193, 194, 199, 201, 202, 206, 208, 209, 210, 212,
             213, 214, 215, 216, 217, 219, 220, 221, 222, 223, 225, 226, 229,
             230, 233, 234, 236, 237, 238, 239, 240, 241, 242, 245, 246, 248,
             249, 251, 253)


diuss_sum = diuss_sum +1

ind.sup_diuss_sum <-x_all[diuss_sum,] 
ind.sup.coord_diuss_sum <- predict(pca, newdata = ind.sup_diuss_sum)

# Plot of active individuals
p_diuss_sum <- fviz_pca_ind(pca, repel = TRUE,geom = "point")

# Add supplementary individuals
fviz_add(p_diuss_sum, ind.sup.coord_diuss_sum, color ="darkgray",addlabel = F)



#------------------------#
#  All sample selection  #
#------------------------#

## Plot
ind.sup_all <-x_all
ind.sup.coord_all <- predict(pca, newdata = ind.sup_all)

# Plot of active individuals
p_all <- fviz_pca_ind(pca, repel = TRUE,geom = "point")

# Add supplementary individuals
fviz_add(p_all, ind.sup.coord_all, color ="purple",addlabel = F)
