########################
# Define colors
color_1 <- "deepskyblue2"
color_2 <- "seagreen2"
color_3 <- "orange2"
color_4 <- "darkorchid4"
color_5 <- "firebrick2"
options(digits=4)

library(ggplot2)
library(factoextra)
library(cluster)

# Size of object
library(pryr)
########################


load('rda/clinical_trial_complete.rda')
X <- clinical_trial_complete[,-c(1,11,13,14,15,16,18)]
X <- X[,-c(3,6,7,9,10)]
X_quan <- X[,-c(1,3,6)]  #col 6 is death

# Categorical variables
Y <- X[,c(1,3,6)] # col 6 is death
X <- X_quan



#### toys examples
#X <- X[1:100,]
#Y <- Y[1:100,]

rm(clinical_trial_complete,X_quan)


##################################################################################################################
# PCA for graphical representaion
##################################################################################################################

X_pcs <- prcomp(X, scale=TRUE)
colors_death <- c(color_3,color_2)[1*(Y[,3]==0)+1]
plot(X_pcs$x[,1:2],pch=19,col=colors_death)
rm(colors_death)


##################################################################################################################
# Agglomerative hierarchical clustering analysis for the NCI60 data set
##################################################################################################################

# We use the Gower distance

C <- cbind(X,Y)
summary(C)
rm(X,Y)  # No needed more?

##################################################################################################################
# Compute the Manhattan distance matrix between the observations in the data matrix


#C_Gower <- daisy(C, metric = "gower", stand = FALSE)
load('rda/C_Gower.rda')
object_size(C_Gower)

##################################################################################################################
# Single linkage

linkages_objs <- list()

single_C <- hclust(C_Gower, method = "single")
linkages_objs[[1]] <- single_C
save(linkages_objs, file = 'rda/linkages_objs.rda')
# Plot dendogram of the solution for k=5 as in K-means

## SIMPLE VERSION
##plot(single_C,main="Single linkage",cex=0.8, xlab="", sub="")
##rect.hclust(single_C,k=5,border=color_1)


single_hcd <- as.dendrogram(single_C)

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/11_H_Single_C.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Customized plot; remove labels
plot(single_hcd, ylab = "Height", nodePar = nodePar, leaflab = "none", main="Single linkage \n 5 cluster")
rect.hclust(single_C, k=5,border=color_1)

# Step 3: Run dev.off() to create the file!
dev.off()

# See the assignment

cl_single_C <- cutree(single_C, 5)
#cl_single_C
table(cl_single_C)

# Make a plot of the first two PCs with the five clusters

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/12_H_cluster_Single.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_single_C <- c(color_1,color_2,color_3,color_4,color_5)[cl_single_C]
plot(X_pcs$x[,1:2],
     pch=19,
     col=colors_single_C,
     main = "PCA graph of individuals \n First 2 PC with 5 clusters using Single linkage",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")


# Step 3: Run dev.off() to create the file!
dev.off()

# BOTH GRAPHS

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/13_H_cluster_Single.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_single_C <- c(color_1,color_2,color_3,color_4,color_5)[cl_single_C]
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Customized plot; remove labels
par(mfrow=c(1,2))
plot(single_hcd, ylab = "Height", nodePar = nodePar, leaflab = "none", main="Single linkage \n 5 cluster")
rect.hclust(single_C, k=5,border=color_1)
plot(X_pcs$x[,1:2],
     pch=19,
     col=colors_single_C,
     main = "PCA graph of individuals \n First 2 PC with 5 clusters using Single linkage",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")


# Step 3: Run dev.off() to create the file!
dev.off()

# Have a look at the silhouette

sil_single_C <- silhouette(cl_single_C,C_Gower)
linkages_objs[[2]] <- sil_single_C
save(linkages_objs, file = 'rda/linkages_objs.rda')

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/14_Single_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# Step 2: Create the plot with R code
plot(sil_single_C, col = color_1)
# Step 3: Run dev.off() to create the file!
dev.off()

# This solution is awful

rm(sil_single_C,colors_single_C,cl_single_C,single_C, single_hcd)

##################################################################################################################
# Complete linkage

complete_C <- hclust(C_Gower, method ="complete")
linkages_objs[[3]] <- complete_C
save(linkages_objs, file = 'rda/linkages_objs.rda')

complete_hcd <- as.dendrogram(complete_C)
# See the assignment

cl_complete_C <- cutree(complete_C, 5)
#cl_complete_C
table(cl_complete_C)

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/15_cluster_complete.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_complete_C <- c(color_1,color_2,color_3,color_4,color_5)[cl_complete_C]
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Customized plot; remove labels
par(mfrow=c(1,2))
plot(complete_hcd, ylab = "Height", nodePar = nodePar, leaflab = "none", main="Complete linkage \n 5 cluster")
rect.hclust(complete_C, k=5,border=color_1)
plot(X_pcs$x[,1:2],
     pch=19,
     col = colors_complete_C,
     main = "PCA graph of individuals \n First 2 PC with 5 clusters using complete linkage",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")


# Step 3: Run dev.off() to create the file!
dev.off()

# Have a look at the silhouette

sil_complete_C <- silhouette(cl_complete_C,C_Gower)
linkages_objs[[4]] <- sil_complete_C
save(linkages_objs, file = 'rda/linkages_objs.rda')

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/16_complete_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# Step 2: Create the plot with R code
plot(sil_complete_C,col=color_1)
# Step 3: Run dev.off() to create the file!
dev.off()

# This solution is better but not as good as the ones from partition methods

rm(sil_complete_C,colors_complete_C,cl_complete_C,complete_C, complete_hcd)

##################################################################################################################
# Average linkage

average_C <- hclust(C_Gower,method="average")
linkages_objs[[5]] <- average_C
save(linkages_objs, file = 'rda/linkages_objs.rda')


average_hcd <- as.dendrogram(average_C)

# See the assignment

cl_average_C <- cutree(average_C,5)
#cl_average_C
table(cl_average_C)

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/17_cluster_average.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_average_C <- c(color_1,color_2,color_3,color_4,color_5)[cl_average_C]
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Customized plot; remove labels
par(mfrow=c(1,2))
plot(average_hcd, ylab = "Height", nodePar = nodePar, leaflab = "none", main="Average linkage \n 5 cluster")
rect.hclust(average_C, k=5,border=color_1)
plot(X_pcs$x[,1:2],
     pch=19,
     col = colors_average_C,
     main = "PCA graph of individuals \n First 2 PC with 5 clusters using average linkage",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")


# Step 3: Run dev.off() to create the file!
dev.off()


# Have a look at the silhouette

sil_average_C <- silhouette(cl_average_C,C_Gower)
linkages_objs[[6]] <- sil_average_C
save(linkages_objs, file = 'rda/linkages_objs.rda')


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/18_average_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# Step 2: Create the plot with R code
plot(sil_average_C,col=color_1)
# Step 3: Run dev.off() to create the file!
dev.off()

# This solution is not good either

rm(sil_average_C,colors_average_C,cl_average_C,average_C, average_hcd)

##################################################################################################################
# Ward linkage

ward_C <- hclust(C_Gower,method="ward")
linkages_objs[[7]] <- ward_C
save(linkages_objs, file = 'rda/linkages_objs.rda')

ward_hcd <- as.dendrogram(ward_C)

# See the assignment

cl_ward_C <- cutree(ward_C,5)
#cl_ward_C
table(cl_ward_C)

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/19_cluster_ward.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_ward_C <- c(color_1,color_2,color_3,color_4,color_5)[cl_ward_C]
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Customized plot; remove labels
par(mfrow=c(1,2))
plot(ward_hcd, ylab = "Height", nodePar = nodePar, leaflab = "none", main="Ward linkage \n 5 cluster")
rect.hclust(ward_C, k=5,border=color_1)
plot(X_pcs$x[,1:2],
     pch=19,
     col = colors_ward_C,
     main = "PCA graph of individuals \n First 2 PC with 5 clusters using ward linkage",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")


# Step 3: Run dev.off() to create the file!
dev.off()

# Have a look at the silhouette

sil_ward_C <- silhouette(cl_ward_C,C_Gower)
linkages_objs[[8]] <- sil_ward_C
save(linkages_objs, file = 'rda/linkages_objs.rda')


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/20_ward_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# Step 2: Create the plot with R code
plot(sil_ward_C,col=color_1)
# Step 3: Run dev.off() to create the file!
dev.off()

# This solution is probably the best one among the agglomerative hierarchical clustering methods

rm(sil_ward_C,colors_ward_C,cl_ward_C,ward_C, ward_hcd, nodePar)

##################################################################################################################
# Divisive hierarchical clustering analysis for the NCI60 data set
##################################################################################################################

# Delete objects
rm(C, C_Gower)

# Only works for quantitative
diana_X <- diana(X, metric = "manhattan")
linkages_objs[[9]] <- diana_X
save(linkages_objs, file = 'rda/linkages_objs.rda')

# Plot dendogram of the solution

plot(diana_X,main="DIANA")

# Hit two times Return to see the dendogram
# The heights here are the diameters of the clusters before splitting
# Take k=5

rect.hclust(diana_X, k = 5, border = color_1)

# See the assignment

cl_diana_X <- cutree(diana_X,5)
#cl_diana_X
table(cl_diana_X)

# Make a plot of the first two PCs with the five clusters

colors_diana_X <- c(color_1,color_2,color_3,color_4,color_5)[cl_diana_X]
plot(X_pcs$x[,1:2],pch=19,col=colors_diana_X,main="First two PCs for the cancer cells data set",xlab="First PC",ylab="Second PC")

# Have a look at the silhouette

sil_diana_C <- silhouette(cl_diana_C,C_Gower)
plot(sil_diana_C,col=color_1)

# This solution is probably the best one among the hierarchical clustering methods