#===========#
# Library   #
#===========#

library(AlgDesign)
library(tidyverse)
library(dplyr)
library(factoextra)
library(readxl)
#=====================#
# Importing Dataset   #
#=====================#

pig_data <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/0pig_data.xlsx")

## X Variable ##
x_initial = pig_data[,c(0:141)]

## Y variable
y = pig_data[,c(142:151)]
y_labels = c("pH","DM","OM","total_N","ammonium_N","P","K","Na","Ca","Mg")
colnames(y) = y_labels

## Data with labeled Y
data = cbind(x_initial,y)

#==========================#
# Calibration Dataset      #
#==========================#

#drop_obs = [48,89,53] outliers 
x_all_keep <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/0x_all_keep_pig.xlsx")

y_all <- read_excel("C:/Users/ricar/Notebook_files/Thesis/Models_for_presentation/0y_all_keep_pig.xlsx")
colnames(y_all) = "DM"

#-----------------------------------------------------------------------------#

#=======================#
# Principal Components  #
#=======================#

pca <- prcomp(x_all_keep, scale = TRUE)
fviz_eig(pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
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


#============================================#
#  Plot of selected samples in PCA space     #
#============================================#

# As R deos not start a df from 0 we sum 1 to each position


#------------------------#
#  Ks sample selection   #
#------------------------#

ks=c( 1,   3,   4,   9,  10,  12,  14,  15,  16,  19,  20,  21,  22,
      25,  27,  28,  30,  32,  35,  40,  41,  45,  46,  50,  51,  53,
      54,  55,  56,  57,  61,  66,  68,  70,  71,  73,  74,  75,  79,
      80,  84,  86,  87,  88,  93,  94,  98,  99, 101, 103, 104, 105,
      107, 109, 110, 112, 115, 116, 120, 121, 122, 125, 126, 128, 130,
      132, 134, 135, 139, 140, 141, 142, 143, 144, 153, 157, 159, 163,
      165, 166, 167, 169, 171, 172, 173, 174, 175, 176, 177, 178, 179,
      181, 184, 186, 187, 188, 190, 193, 195, 196, 197, 200, 203, 204,
      205, 207, 208, 210, 211, 212, 213, 215, 216, 218, 221, 223, 224,
      225, 226, 227, 233, 234, 235, 237, 239, 240, 241, 242, 244, 247,
      250, 251, 252)
Ks = Ks+1

ind.sup_ks <-data[ks,]
ind.sup.coord_ks <- predict(pca, newdata = ind.sup_ks)

# Plot of active individuals
p_ks <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_ks, ind.sup.coord_ks, color ="blue")


#------------------------#
#  dup sample selection  #
#------------------------#


dup=c( 0,   1,   5,   7,   9,  10,  11,  13,  16,  17,  18,  19,  22,
       23,  25,  26,  28,  29,  30,  35,  40,  42,  43,  44,  47,  51,
       52,  54,  57,  60,  61,  62,  66,  67,  68,  74,  75,  76,  78,
       79,  80,  81,  82,  83,  85,  86,  87,  90,  91,  92,  94,  95,
       96,  97,  98,  99, 100, 103, 106, 109, 110, 111, 113, 114, 117,
       118, 119, 121, 123, 124, 125, 126, 127, 130, 131, 135, 136, 137,
       140, 143, 145, 147, 152, 154, 156, 157, 158, 159, 161, 163, 164,
       166, 169, 171, 177, 178, 179, 180, 183, 184, 187, 188, 191, 192,
       195, 196, 198, 200, 202, 206, 207, 208, 209, 212, 214, 215, 217,
       221, 222, 223, 225, 226, 227, 229, 232, 234, 237, 242, 243, 245,
       246, 248, 252)

# As R deos not start a df from 0 we sum 1 to each position
dup = dup +1

ind.sup_dup <-data[dup,] 
ind.sup.coord_dup <- predict(pca, newdata = ind.sup)

# Plot of active individuals
p_dup <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_dup, ind.sup.coord_dup, color ="orange")

#------------------------#
#  puch sample selection  #
#------------------------#


puch=c( 1,   6,   9,  12,  14,  15,  19,  20,  21,  25,  27,  28,  29,
        30,  31,  32,  35,  37,  40,  45,  46,  50,  51,  54,  55,  57,
        61,  65,  66,  68,  70,  73,  74,  75,  76,  79,  80,  84,  86,
        87,  88,  90,  93,  94,  95,  99, 101, 103, 104, 105, 107, 108,
        109, 110, 112, 115, 120, 121, 122, 124, 125, 126, 128, 130, 131,
        132, 134, 139, 140, 141, 142, 143, 144, 151, 153, 157, 159, 163,
        165, 167, 169, 170, 171, 172, 173, 175, 176, 177, 178, 179, 181,
        184, 185, 186, 187, 188, 192, 193, 197, 200, 201, 203, 204, 205,
        207, 208, 209, 210, 211, 212, 213, 215, 216, 218, 221, 223, 224,
        225, 226, 227, 230, 231, 233, 234, 235, 237, 239, 240, 242, 244,
        247, 250, 252)


puch = puch +1

ind.sup_puch <-data[puch,] 
ind.sup.coord_puch <- predict(pca, newdata = ind.sup_puch)

# Plot of active individuals
p_puch <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_puch, ind.sup.coord_puch, color ="green")

#------------------------#
#  clus sample selection  #
#------------------------#


clus=c( 1,   2,   6,   9,  11,  12,  13,  15,  18,  19,  20,  22,  24,
        25,  26,  27,  28,  29,  30,  31,  32,  33,  35,  36,  37,  40,
        42,  43,  45,  46,  47,  49,  53,  54,  57,  59,  61,  62,  63,
        65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  77,  78,
        79,  80,  82,  83,  85,  86,  88,  90,  94,  97,  98,  99, 102,
        104, 105, 107, 109, 117, 121, 122, 123, 124, 125, 126, 130, 133,
        134, 136, 137, 140, 142, 144, 150, 153, 157, 159, 160, 163, 165,
        167, 169, 172, 175, 176, 177, 178, 179, 181, 182, 183, 184, 185,
        186, 188, 192, 193, 194, 195, 200, 203, 204, 206, 208, 211, 212,
        213, 216, 218, 219, 224, 225, 226, 227, 228, 231, 233, 238, 239,
        244, 247, 248)


clus = clus +1

ind.sup_clus <-data[clus,] 
ind.sup.coord_clus <- predict(pca, newdata = ind.sup_clus)

# Plot of active individuals
p_clus <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_clus, ind.sup.coord_clus, color ="red")


#------------------------#
#  DOE sample selection  #
#------------------------#


## DOE to calculate the Z matrix
PCAX=pca$x[,c(0:25)]

Z_0=optFederov(data = PCAX, nTrials = 133, center=TRUE,approximate=FALSE,criterion="D")
Z=as.vector(Z_0$rows)

# Z List of samples
Z
# Binary for of Z in X
DOE = rep(0,252)
j=0

for (i in 0:253){
  if (i %in% Z){
    DOE[i] = 1
    j=j+1} }

DOE # Binary Z  
Z # Selected samples
j # num of 1 in DOE


dput(DOE)
dput(Z)

## Plot
ind.sup_DOE <-data[Z,]
ind.sup.coord_DOE <- predict(pca, newdata = ind.sup)

# Plot of active individuals
p_DOE <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_DOE, ind.sup.coord_DOE, color ="brown")



#-----------------------------#
#  diuss_max sample selection  #
#-----------------------------#


diuss_max=c( 7,   8,  11,  12,  19,  21,  23,  26,  30,  32,  36,  37,  39,
             40,  44,  45,  46,  47,  48,  50,  51,  53,  54,  56,  57,  58,
             63,  64,  66,  70,  72,  73,  74,  75,  78,  80,  82,  84,  86,
             94,  97,  99, 100, 101, 102, 104, 105, 107, 110, 111, 112, 113,
             114, 116, 117, 120, 121, 124, 126, 132, 133, 135, 136, 138, 141,
             144, 145, 146, 148, 151, 153, 154, 156, 159, 161, 164, 165, 166,
             169, 171, 172, 174, 175, 176, 177, 178, 180, 181, 182, 184, 185,
             186, 187, 188, 189, 190, 191, 192, 194, 195, 196, 197, 198, 200,
             201, 205, 207, 210, 211, 213, 215, 217, 218, 219, 220, 222, 228,
             229, 230, 232, 233, 235, 237, 238, 239, 240, 243, 245, 246, 248,
             249, 250, 251)


diuss_max = diuss_max +1

ind.sup_diuss_max <-data[diuss_max,] 
ind.sup.coord_diuss_max <- predict(pca, newdata = ind.sup_diuss_max)

# Plot of active individuals
p_diuss_max <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_diuss_max, ind.sup.coord_diuss_max, color ="pink")

#-----------------------------#
#  diuss_sum sample selection  #
#-----------------------------#


diuss_sum=c( 6,   7,   9,  11,  12,  13,  14,  17,  18,  19,  21,  23,  26,
             28,  29,  31,  32,  35,  37,  39,  41,  42,  44,  46,  47,  48,
             50,  51,  53,  54,  55,  56,  57,  58,  60,  61,  62,  66,  70,
             72,  73,  74,  75,  76,  77,  79,  80,  81,  84,  86,  87,  88,
             89,  90,  96,  98,  99, 101, 102, 103, 104, 105, 106, 108, 110,
             111, 112, 114, 116, 120, 122, 123, 124, 132, 134, 136, 140, 141,
             142, 144, 146, 148, 153, 156, 157, 166, 167, 169, 171, 172, 173,
             176, 177, 178, 179, 182, 183, 184, 185, 189, 190, 191, 192, 194,
             195, 196, 197, 198, 200, 201, 202, 206, 210, 211, 215, 219, 222,
             223, 229, 230, 231, 232, 233, 236, 238, 239, 240, 243, 244, 245,
             248, 249, 252)


diuss_sum = diuss_sum +1

ind.sup_diuss_sum <-data[diuss_sum,] 
ind.sup.coord_diuss_sum <- predict(pca, newdata = ind.sup_diuss_sum)

# Plot of active individuals
p_diuss_sum <- fviz_pca_ind(pca, repel = TRUE)

# Add supplementary individuals
fviz_add(p_diuss_sum, ind.sup.coord_diuss_sum, color ="darkgrey")



