# Two possibilities: 

## 1. Obtain PCs on the whole data set 
## 2. A PC's for two (each) groups
# See what happens if we obtain PCs on the whole data set
# Use the sample correlation matrix


##################################################################################################################
# Analysis 1: All the population
##################################################################################################################

load('rda/clinical_trial_complete.rda')
X <- clinical_trial_complete[,-c(1,11,13,14,15,16,18)]
X <- X[,-c(3,6,7,9,10)]
X_quan <- X[,-c(1,3,6)]

# Categorical variables
Y <- X[,c(1,3,6)]
X <- X_quan

# Remove the file
rm(clinical_trial_complete)
rm(X_quan)

# Define colors for plots

color_1 <- "deepskyblue2"
color_2 <- "seagreen2"
color_3 <- "orange2"
color_4 <- "darkorchid4"


#Define sample size and the dimension of the data matrix

n <- nrow(X)
p <- ncol(X)


##################################################################################################################

X_pcs <- prcomp(X,scale=TRUE)

##################################################################################################################
# PC scores

dim(X_pcs$x)
head(X_pcs$x)


# Plot the PC for each group
colnames(Y)

######### SEX
colors_sex <- c(color_2,color_3)[1*(Y[,1]=="male")+1]
plot(X_pcs$x[,1:2],pch=19,col=colors_sex, main = "First 2 PC by Sex")


######### PAIN
colors_pain <- c(color_2,color_3,color_4)[1 + ifelse(Y[,2]== "blunt", 1,
                                              ifelse(Y[,2]=="penetrating",2,
                                                     0))]
plot(X_pcs$x[,1:2],pch=19,col=colors_pain)


######### death 1 (death)

colors_death <- c(color_3,color_2)[1*(Y[,3]==0)+1]
plot(X_pcs$x[1:100,1:2],
     pch=19,
     col=colors_death, 
     main = "PCA graph of individuals \n First 2 PC by death",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")
abline(v=0, lty = 2)
abline(h=0, lty = 2)


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/pca_death.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
plot(X_pcs$x[,1:2],     pch=19,
     col=colors_death, 
     main = "PCA graph of individuals \n First 2 PC by death",
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)")
abline(v=0, lty = 2)
abline(h=0, lty = 2)


# Step 3: Run dev.off() to create the file!
dev.off()






##################################################################################################################
# Interpretation of the first PC: Weights for the first PC

plot(1:p,X_pcs$rotation[,1],pch=19,col=color_1,main="Weights for the first PC",
     xlab="Variables",ylab="Score", xlim = c(1, 9), ylim = c(-1, 1))
abline(h=0)
text(1:p,X_pcs$rotation[,1],labels=colnames(X),pos=1,col=color_4,cex=0.75)


# Interpretation of the second PC: Weights for the second PC

plot(1:p,X_pcs$rotation[,2],pch=19,col=color_1,main="Weights for the second PC",
     xlab="Variables",ylab="Score",xlim = c(1, 9), ylim = c(-1, 1))
abline(h=0)
text(1:p,X_pcs$rotation[,2],labels=colnames(X),pos=1,col=color_4,cex=0.75)


##################################################################################################################
# Have a look at the important variables in the first two PCs
# Note the different groups in the data
# The radius is arbitrary

plot(X_pcs$rotation[,1:2],
     pch=19,
     col=color_1,
     main = "PCA graph of variables \n Weights for the first two PCs",
     xlim = c(-1,1),
     ylim = c(-1,1),
     xlab = "PC1 (19.46%)",
     ylab = "PC2 (17.22%)"
     )
abline(h=0,v=0)
text(X_pcs$rotation[,1:2],labels=colnames(X),pos=1,col=color_4,cex=0.75)
library(plotrix)
draw.circle(0,0,1,border=color_4,lwd=1)


##################################################################################################################
# The biplot is an alternative way to plot points and the first two PCs together 
# However, it is only useful when the data set is not too large




# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/pca_biplot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
biplot(X_pcs,col=c(color_1,color_4),cex=c(0.5,0.8))

# Step 3: Run dev.off() to create the file!
dev.off()


#######################################
# Have a look at these eigenvalues

# Screeplot

library(factoextra)
fviz_eig(X_pcs,ncp=17,addlabels=T,barfill=color_1,barcolor=`color_4`)

get_eigenvalue(X_pcs)

# 4 PC's explain 62% of the data




##################################################################################################################
# Plot the scores (the four PCs)



# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/pca_4pc2.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
pairs(X_pcs$x[,1:4],col=colors_death,pch=19,main="The first four PCs",lower.panel=NULL)

# Step 3: Run dev.off() to create the file!
dev.off()





# The PCs show that the second PC is the key to show the two groups

# The first one is a PC that creates a ranking of colleges. See the first ten colleges

sort(X_pcs$x[,1],decreasing=TRUE)[1:10]

##################################################################################################################
# Plot the correlations between the original data set and the scores
# This is useful to understand also which are the important variables in the data set in terms of variability
# First, plot all the PCs
# Second, plot only the first four PCs

library(corrplot)
corrplot(cor(X,X_pcs$x),is.corr=T)
corrplot(cor(X,X_pcs$x[,1:4]),is.corr=T)

##################################################################################################################
##################################################################################################################
# playing with factominer
##################################################################################################################
##################################################################################################################
library(FactoMineR)


res<-PCA(X, scale.unit=TRUE, ncp=4, graph = FALSE)
plot(res, axes=c(1, 2), choix="ind", col.ind="red",new.plot=TRUE)
plot(res, axes=c(1, 2), choix="var", col.var="blue",new.plot=TRUE)

##################################################################################################################
##################################################################################################################
# Analysis 2: All the population
##################################################################################################################
##################################################################################################################

# death

load('rda/clinical_trial_complete.rda')
X <- clinical_trial_complete[,-c(1,11,13,14,15,16,18)]

## Population: Death
X_death <- X[X$death == 1,]
X_death <- X_death[,-c(3,6,7,9,10)]

Y_death <- X_death[,c(1,3,6)]
X_death <- X_death[,-c(1,3,6)]

n_death <- nrow(X_death)
p_death <- ncol(X_death)

X_pcs_death <- prcomp(X_death,scale=TRUE)


## Survival

X_surv <- X[X$death == 0,]
X_surv <- X_surv[,-c(3,6,7,9,10)]

Y_surv <- X_surv[,c(1,3,6)]
X_surv <- X_surv[,-c(1,3,6)]

n_surv <- nrow(X_surv)
p_surv <- ncol(X_surv)


X_pcs_surv <- prcomp(X_surv,scale=TRUE)



######### SEX
# blue female
# gray: male

colors_sex <- c("blue","gray")[1*(Y_death[,1]=="male")+1]
plot(X_pcs_death$x[,1:2],pch=19,col=colors_sex, main = "First 2 PC by Sex death patients")
plot(X_pcs_surv$x[,1:2],pch=19,col=colors_sex, main = "First 2 PC by Sex survival patients")



# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/pca_groups.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
par(mfrow=c(1,2))
plot(X_pcs_death$x[,1:2],pch=19,col=colors_sex, 
     main = "PCA graph of death patients \n First 2 PC by sex",
     xlab = "PC1 (19.1%)",
     ylab = "PC2 (16.1%)")
plot(X_pcs_surv$x[,1:2],pch=19,col=colors_sex, 
     main = "PCA graph of survival patients \n First 2 PC by sex",
     xlab = "PC1 (19.8%)",
     ylab = "PC2 (17.5%)")
par(mfrow=c(1,1))

# Step 3: Run dev.off() to create the file!
dev.off()






##################################################################################################################
# Interpretation of the first PC: Weights for the first PC

par(mfrow=c(1,2))
plot(1:p_death,X_pcs_death$rotation[,1],pch=19,col=color_1,main="Weights for the first PC death patients",
     xlab="Variables",ylab="Score", ylim = c(-1,0.6), xlim = c(0,9))
abline(h=0)
text(1:p_death,X_pcs_death$rotation[,1],labels=colnames(X_death),pos=1,col=color_4,cex=0.75)

plot(1:p_surv,X_pcs_surv$rotation[,1],pch=19,col=color_1,main="Weights for the first PC survival patients",
     xlab="Variables",ylab="Score" , ylim = c(-0.6,0.7), xlim = c(0,9))
abline(h=0)
text(1:p_surv,X_pcs_surv$rotation[,1],labels=colnames(X_surv),pos=1,col=color_4,cex=0.75)
par(mfrow=c(1,1))



# Interpretation of the second PC: Weights for the second PC
par(mfrow=c(1,2))
plot(1:p_death,X_pcs_death$rotation[,2],pch=19,col=color_1,main="Weights for the second PC deaths patients",
     xlab="Variables",ylab="Score", ylim = c(-0.3,0.7), xlim = c(0,9))
abline(h=0)
text(1:p_death,X_pcs_death$rotation[,2],labels=colnames(X_death),pos=1,col=color_4,cex=0.75)

plot(1:p_surv,X_pcs_surv$rotation[,2],pch=19,col=color_1,main="Weights for the second PC deaths patients",
     xlab="Variables",ylab="Score", ylim = c(-0.7,0.6), xlim = c(0,9))
abline(h=0)
text(1:p_surv,X_pcs_surv$rotation[,2],labels=colnames(X_surv),pos=1,col=color_4,cex=0.75)

par(mfrow=c(1,1))


## both plots
par(mfrow=c(1,2))
plot(X_pcs_death$rotation[,1:2],
     pch=19,
     col=color_1,
     main = "PCA graph of variables (death patients) \n Weights for the first two PCs",
     xlim = c(-1,1),
     ylim = c(-1,1),
     xlab = "PC1 (19.1%)",
     ylab = "PC2 (16.1%)"
)
abline(h=0,v=0)
text(X_pcs_death$rotation[,1:2],labels=colnames(X_death),pos=1,col=color_4,cex=0.75)
draw.circle(0,0,1,border=color_4,lwd=1)
## plot 2
plot(X_pcs_surv$rotation[,1:2],
     pch=19,
     col=color_1,
     main = "PCA graph of variables (survival patients) \n Weights for the first two PCs",
     xlim = c(-1,1),
     ylim = c(-1,1),
     xlab = "PC1 (19.8%)",
     ylab = "PC2 (17.5%)"
)
abline(h=0,v=0)
text(X_pcs_surv$rotation[,1:2],labels=colnames(X_surv),pos=1,col=color_4,cex=0.75)
draw.circle(0,0,1,border=color_4,lwd=1)
par(mfrow=c(1,1))



## biplot

# Step 1: Call the pdf command to start the plot
pdf(file = "figure_output/pca_group_biplot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
par(mfrow=c(1,2))
biplot(X_pcs_death,col=c(color_1,color_4),cex=c(0.5,0.8))
biplot(X_pcs_surv,col=c(color_1,color_4),cex=c(0.5,0.8))
# Step 3: Run dev.off() to create the file!
dev.off()



## EIGENVALUES
wpc_death <-  fviz_eig(X_pcs_death,ncp=8,addlabels=T,barfill=color_1,barcolor=`color_4`)
fviz_eig(X_pcs_surv,ncp=8,addlabels=T,barfill=color_1,barcolor=`color_4`)

