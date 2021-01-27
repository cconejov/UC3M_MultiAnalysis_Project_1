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

########################


##################################################################################################################
# Analysis 1: All the population
##################################################################################################################

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
## Reorder X
#X <- X[,c(4,1,2,3,5:9)]


# Remove the files
rm(clinical_trial_complete,X_quan)

# SUMMARY

## ONLY QUANTITATIVE
summary(X)

## QUANTITATIVE (GOWER DISTANCE)
summary(Y)

##################################################################################################################
# PCA (Visualization)
##################################################################################################################

X_pcs <- prcomp(X, scale = TRUE)

# Same as 
#Y <- scale(X)
#Y_pcs <- prcomp(Y)
#summary(Y_pcs$x)
#head(Y_pcs$x)
#rm(Y_pcs, Y)


summary(X_pcs$x)
head(X_pcs$x)


##################################################################################################################
# PC scores

dim(X_pcs$x)
head(X_pcs$x)

##################################################################################################################
# Make a plot of the first two PCs

colors_death <- c(color_3,color_2)[1*(Y[,3]==0)+1]
plot(X_pcs$x[,1:2],pch=19,col=colors_death)

##################################################################################################################
# CLUSTER ANALYSIS
##################################################################################################################

##################################################################################################################
# a) Partitional clustering
##################################################################################################################

## a.1) Selecting k with different methods

### a.1.1) Select K with WSS
fviz_nbclust(X, kmeans, method="wss", k.max = 10) # takes time save as pdf

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/1_k_wss_all.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
fviz_nbclust(X, kmeans, method="wss", k.max = 10) # takes time save as pdf

# Step 3: Run dev.off() to create the file!
dev.off()

# The plot suggests is totally unclear k = 3 or 4 are possible options

##################################################################################################################
### a.1.2) # Select K with Silhouette

fviz_nbclust(X,kmeans,method="silhouette",k.max=10)

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/2_k_silhouette_all.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
fviz_nbclust(X,kmeans,method="silhouette",k.max=10)# takes time save as pdf

# Step 3: Run dev.off() to create the file!
dev.off()

# The plot suggests K=2 

##################################################################################################################
### a.1.3) Select K with the gap statistic

gap_stat <- clusGap(X,FUN=kmeans,K.max=5,B=10) # boostrap. Take time
fviz_gap_stat(gap_stat,
              linecolor = "steelblue",
              maxSE = list(method = "firstmax", SE.factor = 1))


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/3_k_firstmax_all.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
fviz_gap_stat(gap_stat,
              linecolor = "steelblue",
              maxSE = list(method = "firstmax", SE.factor = 1))# takes time save as pdf

# Step 3: Run dev.off() to create the file!
dev.off()


# The plot suggests K=1

## ******************************************
## a.2) run k-means with K=2 from the previous 

clustering_models <- list()
kmeans_X <- kmeans(X, centers = 2, iter.max = 1000, nstart = 100)
clustering_models[[1]] <- kmeans_X


# The cluster solution

#kmeans_X$cluster
table(kmeans_X$cluster)

# Make a plot of the first two PCs split in these four clusters

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/4_kmeans_Cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_kmeans_X <- c(color_1,color_2)[kmeans_X$cluster]
plot(X_pcs$x[,1:2],
     pch=19,
     col=colors_kmeans_X,
     #main = "First two PCs for the crash2 data set",
     main = "PCA graph of individuals \n First 2 PC with 2 cluster using kmeans",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")
abline(v=0, lty = 2)
abline(h=0, lty = 2)

# Step 3: Run dev.off() to create the file!
dev.off()



## ******************************************
## a.3) Compute the silhouette

sil_kmeans_X <- silhouette(kmeans_X$cluster, dist(X,"euclidean"))


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/5_kmeans_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
plot(sil_kmeans_X,col=color_1)

# Step 3: Run dev.off() to create the file!
dev.off()

rm(kmeans_X,colors_kmeans_X,sil_kmeans_X)

##################################################################################################################
# b) K-meloids
##################################################################################################################

# ********************
# b.1) PAM
# ********************


pam_X <- pam(X, k = 2, metric = "manhattan", stand = FALSE)
clustering_models[[2]] <- pam_X
save(clustering_models, file = 'rda/clustering_models.rda')

#kmeans_X$cluster
table(pam_X$cluster)


# Make a plot of the first two PCs with the solution


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/6_PAM_cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_pam_X <- c(color_1,color_2)[pam_X$cluster]
plot(X_pcs$x[,1:2],
     pch=19,
     col=colors_pam_X,
     #main = "First two PCs for the crash2 data set",
     main = "PCA graph of individuals \n First 2 PC with 2 clusters using PAM",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")
abline(v=0, lty = 2)
abline(h=0, lty = 2)

# Step 3: Run dev.off() to create the file!
dev.off()

# Have a look at the silhouette


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/7_PAM_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
sil_pam_X <- silhouette(pam_X$cluster,dist(X,method="manhattan"))
plot(sil_pam_X,col=color_1)

# Step 3: Run dev.off() to create the file!
dev.off()

# The solution is close to K-means but it migth be sligthly better


rm(pam_X,sil_pam_X,colors_pam_X)

# ********************
# b.1) CLARA
# ********************


clara_X <- clara(X, k = 2, metric = "manhattan", stand = FALSE)
clustering_models[[3]] <- clara_X
save(clustering_models, file = 'rda/clustering_models.rda')


# Make a plot of the first two PCs with the solution
# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/8_CLARA_cluster.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
colors_clara_X <- c(color_1,color_2)[clara_X$cluster]
plot(X_pcs$x[,1:2],
     pch=19,
     col=colors_clara_X,
     #main = "First two PCs for the crash2 data set",
     main = "PCA graph of individuals \n First 2 PC with 2 clusters using CLARA",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")
abline(v=0, lty = 2)
abline(h=0, lty = 2)

# Step 3: Run dev.off() to create the file!
dev.off()

# Have a look at the silhouette

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/9_CLARA_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
sil_clara_X <- silhouette(clara_X$cluster,dist(X,method="manhattan"))
plot(sil_clara_X,col=color_1)

# Step 3: Run dev.off() to create the file!
dev.off()

# The solution is similar to the one given by PAM but the computational cost is less

rm(clara_X,sil_clara_X,colors_clara_X)

# ********************
# b.3) PAM-GOWER DISTANCE
# ********************

## MIXED DATA:

C <- cbind(X,Y)

C_Gower <- daisy(C, metric = "gower")

# save the C_Gower
library(pryr)
object_size(C_Gower)

save(C_Gower, file = 'rda/C_Gower.rda')


## save the sumarry like an object!!!
summ_C_Gower <- summary(C_Gower)
clustering_models[[4]] <- summ_C_Gower
save(clustering_models, file = 'rda/clustering_models.rda')

# In the output, "I" means quantitative variable while "N" means qualitative variable

##################################################################################################################
# Find the closest clients with the Gower distance

C_Gower_mat <- as.matrix(C_Gower)
rm(C_Gower)
C[which(C_Gower_mat==min(C_Gower_mat[C_Gower_mat!=min(C_Gower_mat)]),arr.ind = TRUE)[1,],]

##################################################################################################################
# Find the most distant clients with the Gower distance

C[which(C_Gower_mat==max(C_Gower_mat[C_Gower_mat!=max(C_Gower_mat)]),arr.ind = TRUE)[1,],]

##################################################################################################################
# Consider K=2:20, run PAM and select K using the silhouette

ncols <- 9 

C_K <- matrix(NA,nrow=1,ncol = ncols)
for (i in 1:ncols){
    print(i)
    pam_C_Gower_mat <- pam(C_Gower_mat, k = i + 1, diss = TRUE)
    C_K[i] <- pam_C_Gower_mat$silinfo$avg.width
}

clustering_models[[5]] <- C_K
save(clustering_models, file = 'rda/clustering_models.rda')

plot(2:10,C_K,pch=19,col="deepskyblue2",xlab="Number of clusters",ylab="Average silhouette")
which.max(C_K)+1

##################################################################################################################
# Run the algorithm for K = 5  and get some information from the results

pam_C_Gower_mat <- pam(C_Gower_mat, k = 5, diss = TRUE)

clustering_models[[6]] <- pam_C_Gower_mat
save(clustering_models, file = 'rda/clustering_models.rda')


load('rda/C_Gower.rda')
C_Gower_mat <- as.matrix(C_Gower)
rm(C_Gower)


# Medoids

C[pam_C_Gower_mat$medoids,]

# Have a look at the silhouette

sil_pam_C_Gower_mat <- silhouette(pam_C_Gower_mat$cluster,C_Gower_mat)
clustering_models[[7]] <- sil_pam_C_Gower_mat
save(clustering_models, file = 'rda/clustering_models.rda')


rm(C_Gower_mat)

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/Report2/10_PAMGOWER_silhouette.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code

plot(sil_pam_C_Gower_mat,col=color_1)

# Step 3: Run dev.off() to create the file!
dev.off()


summary(sil_pam_C_Gower_mat)
