################################################################################
# Read in and combine data
################################################################################
setwd("F:/Dropbox/UM/03.projects/ProstateCancer/256x256_patches/output")
setwd("C:/Users/CHENG-BANG/Dropbox/UM/03.projects/ProstateCancer/256x256_patches/output")

rwdataO <- read.csv("gs_orig.csv")
rwdataA <- read.csv("1A.csv")
rwdataD <- read.csv("1D.csv") 
rwdataH <- read.csv("1H.csv") 
rwdataV <- read.csv("1V.csv")
rmid <- read.csv('blank_patch.csv')

rwdataO$Row <- as.numeric(data.frame(t(sapply(strsplit(rwdataA$Row,split = "slide"),c)))[,2])
rwdataO[rwdataO=="NaN"]<-0
rwdataO <- rwdataO[order(rwdataO$Row),]
nameO <- names(rwdataO)
nameO <- paste("O.",nameO, sep="")
names(rwdataO) <- nameO


rwdataA$Row <- as.numeric(data.frame(t(sapply(strsplit(rwdataA$Row,split = "slide"),c)))[,2])
rwdataA[rwdataA=="NaN"]<-0
rwdataA <- rwdataA[order(rwdataA$Row),]
nameA <- names(rwdataA)
nameA <- paste("A.",nameA, sep="")
names(rwdataA) <- nameA

rwdataD$Row <- as.numeric(data.frame(t(sapply(strsplit(rwdataD$Row,split = "slide"),c)))[,2])
rwdataD[rwdataD=="NaN"]<-0
rwdataD <- rwdataD[order(rwdataD$Row),]
nameD <- names(rwdataD)
nameD <- paste("D.",nameD, sep="")
names(rwdataD) <- nameD

rwdataH$Row <- as.numeric(data.frame(t(sapply(strsplit(rwdataH$Row,split = "slide"),c)))[,2])
rwdataH[rwdataH=="NaN"]<-0
rwdataH <- rwdataH[order(rwdataH$Row),]
nameH <- names(rwdataH)
nameH <- paste("H.",nameH, sep="")
names(rwdataH) <- nameH

rwdataV$Row <- as.numeric(data.frame(t(sapply(strsplit(rwdataV$Row,split = "slide"),c)))[,2])
rwdataV[rwdataV=="NaN"]<-0
rwdataV <- rwdataV[order(rwdataV$Row),]
nameV <- names(rwdataV)
nameV <- paste("V.",nameV, sep="")
names(rwdataV) <- nameV

rwdata <- cbind(rwdataO,rwdataA[,-1],rwdataD[,-1],rwdataH[,-1],rwdataV[,-1])
names(rwdata)[1]<-"Row"
rwdata$GS <- rep(c(3,4,5),times=c(1000,1000,1000))
write.csv(rwdata, file="rwdata.csv", row.names = F)

idx_odd <- seq(3000) %% 2     # The patch id starts from 0. 
rdata <- rwdata[idx_odd == 1, ]
sdata <- rwdata[idx_odd != 1, ]
################################################################################
# 50 Replications of LASSO to identify the critical features
################################################################################
library(glmnet)
# Run cross-validation for the Lasso model
set.seed(1024)
selected_idx <- c()
kn <- 50
for(i in 1:kn){
  tic <- Sys.time()
  rm_col_idx = which(names(rdata) %in% c("Row","GS"))
  cv_model_lasso <- cv.glmnet(x = as.matrix(rdata[, -rm_col_idx]),
                              y = rdata[, "GS"],
                              type.measure = "mse",
                              nfold = 15,
                              alpha = 1)
  
  # The optimal lambda value
  optimal_lambda <- cv_model_lasso$lambda.min
  
  # Fit the model with optimal lambda
  model <- glmnet(x = as.matrix(rdata[, -1]), 
                  y = rdata[, 1], 
                  alpha = 1,
                  lambda = optimal_lambda)
  
  # Get the coefficient estimates
  coefs <- coef(model)
  idx <- which(coefs!=0)
  
  # Find non-zero coefficients (i.e., selected variables)
  selected_vars <- rownames(coefs)[idx]
  if("(Intercept)" %in% selected_vars){
    selected_vars <- selected_vars[selected_vars != "(Intercept)"]
  }
  selected_idx <- c(selected_idx,which(names(rdata) %in% selected_vars))
  print(i)
  print(Sys.time()-tic)
}
################################################################
selected_var_idx <- sort(unique(selected_idx))
Dweights <- table(sort(selected_idx))/kn
Dweights.df <- data.frame(Dweights)

epsilon <- 0.9                   ## <----Decide the threshold 
Dweights2 <- Dweights[Dweights>epsilon]
Dweights2.df <- data.frame(Dweights2)
selected_var_idx2 <- selected_var_idx[Dweights2.df$Var1]


rdata_selected <- rdata[,selected_var_idx]
rdata_selected_w <- rdata_selected*Dweights

sdata_selected <- sdata[,selected_var_idx]
sdata_selected_w <- sdata_selected*Dweights


rdata_selected2 <- rdata[,selected_var_idx2]
rdata_selected_w2 <- rdata_selected2*Dweights2.df$Freq

sdata_selected2 <- sdata[,selected_var_idx2]
sdata_selected_w2 <- sdata_selected2*Dweights2.df$Freq

################################################################################
# Visualization in UMAP
################################################################################
library(umap)
library(ggplot2)
library(rgl)
rumap <- umap(rdata_selected)
sumap <- predict(rumap, sdata_selected)
windows(width = 16, height = 8);par(mfrow=c(1,2));plot(rumap$layout, col=rdata$GS-1, pch=13+rdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Real");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1);plot(sumap, col=sdata$GS-1, pch=13+sdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Synthetic");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1)

rumap.3d <- umap(rdata_selected,n_components=3)
sumap.3d <- predict(rumap.3d, sdata_selected)
plot3d(rumap.3d$layout,col=rdata$GS-1)
plot3d(sumap.3d,col=sdata$GS-1)

rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(rumap.3d$layout, type = "s", size = .5 ,col=rdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap.3d$layout[,1]),min(rumap.3d$layout[,2]),min(rumap.3d$layout[,3]),
                max(rumap.3d$layout[,1]),min(rumap.3d$layout[,2]),min(rumap.3d$layout[,3]),
                min(rumap.3d$layout[,1]),max(rumap.3d$layout[,2]),min(rumap.3d$layout[,3]),
                min(rumap.3d$layout[,1]),min(rumap.3d$layout[,2]),max(rumap.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)


rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(sumap.3d, type = "s", size = .5 ,col=sdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap.3d$layout[,1]),min(rumap.3d$layout[,2]),min(rumap.3d$layout[,3]),
                max(rumap.3d$layout[,1]),min(rumap.3d$layout[,2]),min(rumap.3d$layout[,3]),
                min(rumap.3d$layout[,1]),max(rumap.3d$layout[,2]),min(rumap.3d$layout[,3]),
                min(rumap.3d$layout[,1]),min(rumap.3d$layout[,2]),max(rumap.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)


rumap2 <- umap(rdata_selected_w)
sumap2 <- predict(rumap2, sdata_selected_w)
windows(width = 16, height = 8);par(mfrow=c(1,2));plot(rumap2$layout, col=rdata$GS-1, pch=13+rdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Real");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1);plot(sumap2, col=sdata$GS-1, pch=13+sdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Synthetic");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1)

rumap2.3d <- umap(rdata_selected_w,n_components=3)
sumap2.3d <- predict(rumap2.3d, sdata_selected_w)
## 3D ##
rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(rumap2.3d$layout, type = "s", size = .5 ,col=rdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap2.3d$layout[,1]),min(rumap2.3d$layout[,2]),min(rumap2.3d$layout[,3]),
                max(rumap2.3d$layout[,1]),min(rumap2.3d$layout[,2]),min(rumap2.3d$layout[,3]),
                min(rumap2.3d$layout[,1]),max(rumap2.3d$layout[,2]),min(rumap2.3d$layout[,3]),
                min(rumap2.3d$layout[,1]),min(rumap2.3d$layout[,2]),max(rumap2.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)


rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(sumap2.3d, type = "s", size = .5 ,col=sdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap2.3d$layout[,1]),min(rumap2.3d$layout[,2]),min(rumap2.3d$layout[,3]),
                max(rumap2.3d$layout[,1]),min(rumap2.3d$layout[,2]),min(rumap2.3d$layout[,3]),
                min(rumap2.3d$layout[,1]),max(rumap2.3d$layout[,2]),min(rumap2.3d$layout[,3]),
                min(rumap2.3d$layout[,1]),min(rumap2.3d$layout[,2]),max(rumap2.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
######

rumap3 <- umap(rdata_selected2)
sumap3 <- predict(rumap3, sdata_selected2)
windows(width = 16, height = 8);par(mfrow=c(1,2));plot(rumap3$layout, col=rdata$GS-1, pch=13+rdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Real");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1);plot(sumap3, col=sdata$GS-1, pch=13+sdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Synthetic");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1)

rumap3.3d <- umap(rdata_selected2,n_components=3)
sumap3.3d <- predict(rumap3.3d, sdata_selected2)
## 3D ##
rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(rumap3.3d$layout, type = "s", size = .5 ,col=rdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap3.3d$layout[,1]),min(rumap3.3d$layout[,2]),min(rumap3.3d$layout[,3]),
                max(rumap3.3d$layout[,1]),min(rumap3.3d$layout[,2]),min(rumap3.3d$layout[,3]),
                min(rumap3.3d$layout[,1]),max(rumap3.3d$layout[,2]),min(rumap3.3d$layout[,3]),
                min(rumap3.3d$layout[,1]),min(rumap3.3d$layout[,2]),max(rumap3.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)


rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(sumap3.3d, type = "s", size = .5 ,col=sdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap3.3d$layout[,1]),min(rumap3.3d$layout[,2]),min(rumap3.3d$layout[,3]),
                max(rumap3.3d$layout[,1]),min(rumap3.3d$layout[,2]),min(rumap3.3d$layout[,3]),
                min(rumap3.3d$layout[,1]),max(rumap3.3d$layout[,2]),min(rumap3.3d$layout[,3]),
                min(rumap3.3d$layout[,1]),min(rumap3.3d$layout[,2]),max(rumap3.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
######

rumap4 <- umap(rdata_selected_w2)
sumap4 <- predict(rumap4, sdata_selected_w2)
windows(width = 16, height = 8);par(mfrow=c(1,2));plot(rumap4$layout, col=rdata$GS-1, pch=13+rdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Real");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1);plot(sumap4, col=sdata$GS-1, pch=13+sdata$GS, xlim=c(-10,10), ylim=c(-10,10),xlab="UMAP1", ylab="UMAP2",main="Synthetic");box(lwd=3);legend(-9, 9, legend=c("3", "4","5"),col=c(2,3,4), pch=c(16,17,18),box.lty=1)

rumap4.3d <- umap(rdata_selected_w2,n_components=3)
sumap4.3d <- predict(rumap4.3d, sdata_selected_w2)
## 3D ##
rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(rumap4.3d$layout, type = "s", size = .5 ,col=rdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap4.3d$layout[,1]),min(rumap4.3d$layout[,2]),min(rumap4.3d$layout[,3]),
                max(rumap4.3d$layout[,1]),min(rumap4.3d$layout[,2]),min(rumap4.3d$layout[,3]),
                min(rumap4.3d$layout[,1]),max(rumap4.3d$layout[,2]),min(rumap4.3d$layout[,3]),
                min(rumap4.3d$layout[,1]),min(rumap4.3d$layout[,2]),max(rumap4.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)


rgl.open()
rgl.bg(col = "white")
par3d(windowRect = 10 + c(0, 0, 400, 400))
plot3d(sumap4.3d, type = "s", size = .5 ,col=sdata$GS-1, cex = 1,axes = FALSE, xlab = "", ylab = "", zlab ="")
rgl.viewpoint(theta = -70, phi = 30, fov = 0, zoom = 1)
vec <- matrix(c(min(rumap4.3d$layout[,1]),min(rumap4.3d$layout[,2]),min(rumap4.3d$layout[,3]),
                max(rumap4.3d$layout[,1]),min(rumap4.3d$layout[,2]),min(rumap4.3d$layout[,3]),
                min(rumap4.3d$layout[,1]),max(rumap4.3d$layout[,2]),min(rumap4.3d$layout[,3]),
                min(rumap4.3d$layout[,1]),min(rumap4.3d$layout[,2]),max(rumap4.3d$layout[,3])),
              ncol=3,
              byrow = T)
arrow3d(vec[1,],vec[2,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[3,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)
arrow3d(vec[1,],vec[4,], type="rotation", barblen = 0.05, s = 0.05, col = "black",width = 1/3, thickness = 0.25*width)

################################################################################
# Utilize supervised learning to create a score system (randomForest)
################################################################################
# Install required packages if not already installed
# install.packages("randomForest")
# install.packages("foreach")
# install.packages("doParallel")

# Load the required packages
library(randomForest)
library(foreach)
library(doParallel)

# Set the number of cores to be used
num_cores <- 12

# Initialize the parallel backend using doParallel
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Generate some example data
data <- rdata_selected
rY <- rep(c(3,4,5), times = c(500,500,500))
# Split the data into training and testing sets
set.seed(1024)
tY <- c()
t_pred <- data.frame(Y=c(),Y.hat=c(),E=c())
nf <- 5
for(i in 1:nf){
  tidx3 <- sample(1:500, 0.8 * 500) 
  tidx4 <- sample(501:1000, 0.8 * 500)
  tidx5 <- sample(1001:1500, 0.8 * 500)
  train_indices <- sort(c(tidx3, tidx4, tidx5))
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  train_y <- rY[train_indices]
  test_y <- rY[-train_indices]
  
  rf <- foreach(ntree=rep(3000, 8), .combine=randomForest::combine,
                .multicombine=TRUE, .packages='randomForest') %dopar% {
                  randomForest(as.matrix(train_data), train_y, ntree=ntree)
                }
  predictions <- predict(rf, test_data)
  train_pred <- data.frame(Y=test_y,Y.hat=predictions,E=test_y-predictions)
  t_pred <- rbind(t_pred, train_pred)
}

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()  # Switch back to sequential processing

# Make predictions on the test data using the trained model
rpredictions <- predict(rf, rdata_selected)
rel_pred <- data.frame(Y=rY,Y.hat=rpredictions,E=rY-rpredictions)

spredictions <- predict(rf, sdata_selected)
syn_pred <- data.frame(Y=rY,Y.hat=spredictions,E=rY-spredictions)
################################################################################

# Load the required packages
library(randomForest)
library(foreach)
library(doParallel)

# Set the number of cores to be used
num_cores <- 18

# Initialize the parallel backend using doParallel
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Generate some example data
data <- rdata_selected
rY3<- rep(c(1,0,0), times = c(500,500,500))
rY4<- rep(c(0,1,0), times = c(500,500,500))
rY5<- rep(c(0,0,1), times = c(500,500,500))
# Split the data into training and testing sets
set.seed(1024)
tY3 <- c()
tY4 <- c()
tY5 <- c()
t3_pred <- data.frame(Y=c(),Y.hat=c(),E=c())
t4_pred <- data.frame(Y=c(),Y.hat=c(),E=c())
t5_pred <- data.frame(Y=c(),Y.hat=c(),E=c())
nf <- 5
for(i in 1:nf){
  tidx3 <- sample(1:500, 0.8 * 500) 
  tidx4 <- sample(501:1000, 0.8 * 500)
  tidx5 <- sample(1001:1500, 0.8 * 500)
  train_indices <- sort(c(tidx3, tidx4, tidx5))
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  train_y3 <- rY3[train_indices]
  test_y3 <- rY3[-train_indices]
  train_y4 <- rY4[train_indices]
  test_y4 <- rY4[-train_indices]
  train_y5 <- rY5[train_indices]
  test_y5 <- rY5[-train_indices]
  
  rf3 <- foreach(ntree=rep(3000, 8), .combine=randomForest::combine,
                .multicombine=TRUE, .packages='randomForest') %dopar% {
                  randomForest(as.matrix(train_data), train_y3, ntree=ntree)
                }
  predictions3 <- predict(rf3, test_data)
  train_pred3 <- data.frame(Y=test_y3,Y.hat=predictions3,E=test_y3-predictions3)
  t3_pred <- rbind(t3_pred, train_pred3)
  
  rf4 <- foreach(ntree=rep(3000, 8), .combine=randomForest::combine,
                 .multicombine=TRUE, .packages='randomForest') %dopar% {
                   randomForest(as.matrix(train_data), train_y4, ntree=ntree)
                 }
  predictions4 <- predict(rf4, test_data)
  train_pred4 <- data.frame(Y=test_y4,Y.hat=predictions4,E=test_y4-predictions4)
  t4_pred <- rbind(t4_pred, train_pred4)
  
  rf5 <- foreach(ntree=rep(3000, 8), .combine=randomForest::combine,
                 .multicombine=TRUE, .packages='randomForest') %dopar% {
                   randomForest(as.matrix(train_data), train_y5, ntree=ntree)
                 }
  predictions5 <- predict(rf5, test_data)
  train_pred5 <- data.frame(Y=test_y5,Y.hat=predictions5,E=test_y5-predictions5)
  t5_pred <- rbind(t5_pred, train_pred5)
}

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()  # Switch back to sequential processing

# Make predictions on the test data using the trained model
rpredictions3 <- predict(rf3, rdata_selected)
rel_pred3 <- data.frame(Y=rY3,Y.hat=rpredictions3,E=rY3-rpredictions3)

spredictions3 <- predict(rf3, sdata_selected)
syn_pred3 <- data.frame(Y=rY3,Y.hat=spredictions,E=rY-spredictions)
################################################################################

library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
# training set
windows(width = 16, height = 8)
p0 <- ggplot(data=train_pred, aes(x=Y.hat, group=factor(Y), fill=factor(Y))) +
  geom_density(adjust=1.5, alpha=.3) +
  theme_ipsum() + ggtitle('Real - all')
  
p0

windows(width = 5, height = 5)
p0v <- ggplot(data=train_pred, aes(y=Y.hat, x=Y, group=factor(Y), fill=factor(Y))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.07) +
  labs(x = "GS Group", y = "Prediction Score")+
  guides(fill = guide_legend(title = "GS"))+
  ggtitle('Real Images (Training Set)')+
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=16, face= "bold", colour= "black",hjust = 0.5 ),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=12,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )

p0v

# all real data
windows(width = 16, height = 8)
p1 <- ggplot(data=rel_pred, aes(x=Y.hat, group=factor(Y), fill=factor(Y))) +
  geom_density(adjust=1.5, alpha=.3) +
  theme_ipsum() + ggtitle('Real - all')
p1

windows(width = 16, height = 8)
p1 <- ggplot(data=rel_pred, aes(x=Y.hat, group=factor(Y), fill=factor(Y))) +
  geom_density(adjust=1.5, alpha=.3) +
  theme_ipsum() + ggtitle('Real - all')
p1

windows(width = 5, height = 5)
p1v <- ggplot(data=rel_pred, aes(y=Y.hat, x=Y, group=factor(Y), fill=factor(Y))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.07) +
  labs(x = "GS Group", y = "Prediction Score")+
  guides(fill = guide_legend(title = "GS"))+
  ggtitle('Real Images (All)')+
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=16, face= "bold", colour= "black",hjust = 0.5 ),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=12,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )

p1v

# testing set: real data
windows(width = 16, height = 8)
p2 <- ggplot(data=t_pred, aes(x=Y.hat, group=factor(Y), fill=factor(Y))) +
  geom_density(adjust=1.5, alpha=.3) +
  theme_ipsum() + ggtitle('Real - test')
p2

windows(width = 5, height = 5)
p2v <- ggplot(data=t_pred, aes(y=Y.hat, x=Y, group=factor(Y), fill=factor(Y))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.07) +
  labs(x = "GS Group", y = "Prediction Score")+
  guides(fill = guide_legend(title = "GS"))+
  ggtitle('Real Images (Testing Set)')+
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=16, face= "bold", colour= "black",hjust = 0.5 ),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=12,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )

p2v

# all syn data
windows(width = 16, height = 8)
p3 <- ggplot(data=syn_pred, aes(x=Y.hat, group=factor(Y), fill=factor(Y))) +
  geom_density(adjust=1.5, alpha=.3) +
  theme_ipsum() + ggtitle('Synthetic - all')
p3

windows(width = 5, height = 5)
p3v <- ggplot(data=syn_pred, aes(y=Y.hat, x=Y, group=factor(Y), fill=factor(Y))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.07) +
  labs(x = "GS Group", y = "Prediction Score")+
  guides(fill = guide_legend(title = "GS"))+
  ggtitle('Synthetic Image')+
  theme(
  # LABELS APPEARANCE
  plot.title = element_text(size=14, face= "bold", colour= "black",hjust = 0.5  ),
  axis.title.x = element_text(size=14, face="bold", colour = "black"),    
  axis.title.y = element_text(size=14, face="bold", colour = "black"),    
  axis.text.x = element_text(size=12, face="bold", colour = "black"), 
  # axis.text.y = element_text(size=12,  colour = "black"), # unbold
  axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
  axis.line.x = element_line(color="black", size = 0.3),
  axis.line.y = element_line(color="black", size = 0.3),
  panel.border = element_rect(colour = "black", fill=NA, size=0.3)
)
p3v
################################################################################
# Hotelling T-square test
################################################################################
library(DescTools)
library(QuantPsyc)
library(energy)
testdata<-rbind(rdata_selected,sdata_selected)
# testdata$Group <- rep(c("r","s"),times=c(1500,1500))
# varnames <- names(testdata)
# write.csv(varnames, file="varnames.csv")

library(umap)



testdatapca <- prcomp(testdata)
windows();plot(cumsum(testdatapca$sdev^2)/sum(testdatapca$sdev^2),type="l")
k <- 20
extracted <- testdatapca$x[,1:k]
windows();plot(0:k,cumsum(c(0,testdatapca$sdev[1:k])^2)/sum(testdatapca$sdev^2),
               type="l",
               lwd=3,
               ylim=c(0,1),
               xlab = "# of Principle Components",
               ylab = "% of variability explained",
               main="Principle Components",
               cex.lab = 1.5,
               cex.main = 2);abline(h=.9, lty=2, col=2,lwd=3);abline(h=.95,lty=3, col=3,lwd=3);box(lwd=3);abline(v=4,lwd=2,lty=2,col=2);abline(v=10,lwd=2,lty=3,col=3)
k <- 10
extracted <- testdatapca$x[,1:k]

rdata_pcak <- extracted[1:1500,]
sdata_pcak <- extracted[1501:3000,]
# 
# rmu<-mult.norm(rdata_pcak)
# smu<-mult.norm(sdata_pcak)
# 
# mvnorm.etest(rdata_pcak,R=3000)
# 
# rmu$mult.test
# smu$mult.test

HotellingsT2Test(rdata_pcak, sdata_pcak)
HotellingsT2Test(rdata_pcak[1:500,], sdata_pcak[1:500,])
HotellingsT2Test(rdata_pcak[501:1000,], sdata_pcak[501:1000,])
HotellingsT2Test(rdata_pcak[1001:1500,], sdata_pcak[1001:1500,])

colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
stat_tab <- data.frame(Real_avg=colMeans(rdata_pcak), Real_sd=colSd(rdata_pcak),
                       Syn_avg=colMeans(sdata_pcak), Syn_sd=colSd(sdata_pcak))



HotellingsT2Test(rdata_pcak[1:500,], rdata_pcak[501:1000,])
HotellingsT2Test(rdata_pcak[501:1000,], rdata_pcak[1001:1500,])
HotellingsT2Test(rdata_pcak[1:500,], rdata_pcak[1001:1500,])

HotellingsT2Test(sdata_pcak[1:500,], sdata_pcak[501:1000,])
HotellingsT2Test(sdata_pcak[501:1000,], sdata_pcak[1001:1500,])
HotellingsT2Test(sdata_pcak[1:500,], sdata_pcak[1001:1500,])
################################################################################
library(tidyverse)
library(hrbrthemes)
library(viridis)
# create a dataset
data <- as.data.frame(extracted)
data$Group <- rep(c("Real","Synthetic"),times=c(1500,1500))
# Plot
windows(4,4)
data %>%
  ggplot( aes(x=Group, y=PC1, fill=Group)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="gray", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Spatial Dynamics PC1: Real vs. Synthetic") +
  xlab("")

windows(4,4)
data %>%
  ggplot( aes(x=Group, y=PC2, fill=Group)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Spatial Dynamics PC2: Real vs. Synthetic") +
  xlab("")

windows(4,4)
data %>%
  ggplot( aes(x=Group, y=PC3, fill=Group)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Spatial Dynamics PC3: Real vs. Synthetic") +
  xlab("")

windows(8,6);par(mfrow=c(2,3))
par(mar=c(2.5,2.5,2.5,2.5));boxplot(PC1 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"),main="SRQA-PC1",cex.names=2.5);box(lwd=3)
boxplot(PC2 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC2",cex.names=2.5);box(lwd=3)
boxplot(PC3 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC3",cex.names=2.5);box(lwd=3)
boxplot(PC4 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC4",cex.names=2.5);box(lwd=3)
boxplot(PC5 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC5",cex.names=2.5);box(lwd=3)
boxplot(PC6 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC6",cex.names=2.5);box(lwd=3)

windows(6,6);par(mfrow=c(2,2))
par(mar=c(2.5,2.5,2.5,2.5));boxplot(PC1 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"),main="SRQA-PC1",cex.names=2.5);box(lwd=3)
boxplot(PC2 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC2",cex.names=2.5);box(lwd=3)
boxplot(PC3 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC3",cex.names=2.5);box(lwd=3)
boxplot(PC4 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC4",cex.names=2.5);box(lwd=3)
# boxplot(PC5 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC5",cex.names=2.5);box(lwd=3)
# boxplot(PC6 ~ Group, data = data,col = c("#FFE0B2", "#FFA726"), main="SRQA-PC6",cex.names=2.5);box(lwd=3)


################################################################################
library(tidyverse)
library(ggplot2)
library(scales)

# make some mock data
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
rstat_tab <- data.frame(Mean=colMeans(rdata_pcak), SD=colSd(rdata_pcak))
rstat_tab$upper <- rstat_tab$Mean + 1.96*rstat_tab$SD/sqrt(1500)
rstat_tab$lower <- rstat_tab$Mean - 1.96*rstat_tab$SD/sqrt(1500)
rstat_tab$variable <- factor(rownames(rstat_tab),levels = rownames(rstat_tab))


rstat_tab3 <- data.frame(Mean=colMeans(rdata_pcak[1:500,]), SD=colSd(rdata_pcak[1:500,]))
rstat_tab3$upper <- rstat_tab3$Mean + 1.96*rstat_tab3$SD/sqrt(500)
rstat_tab3$lower <- rstat_tab3$Mean - 1.96*rstat_tab3$SD/sqrt(500)
rstat_tab3$variable <- factor(rownames(rstat_tab3),levels = rownames(rstat_tab3))

rstat_tab4 <- data.frame(Mean=colMeans(rdata_pcak[501:1000,]), SD=colSd(rdata_pcak[501:1000,]))
rstat_tab4$upper <- rstat_tab4$Mean + 1.96*rstat_tab4$SD/sqrt(500)
rstat_tab4$lower <- rstat_tab4$Mean - 1.96*rstat_tab4$SD/sqrt(500)
rstat_tab4$variable <- factor(rownames(rstat_tab4),levels = rownames(rstat_tab4))

rstat_tab5 <- data.frame(Mean=colMeans(rdata_pcak[1001:1500,]), SD=colSd(rdata_pcak[1001:1500,]))
rstat_tab5$upper <- rstat_tab5$Mean + 1.96*rstat_tab5$SD/sqrt(500)
rstat_tab5$lower <- rstat_tab5$Mean - 1.96*rstat_tab5$SD/sqrt(500)
rstat_tab5$variable <- factor(rownames(rstat_tab5),levels = rownames(rstat_tab5))



sstat_tab <- data.frame(Mean=colMeans(sdata_pcak), SD=colSd(sdata_pcak))
sstat_tab$upper <- sstat_tab$Mean + 1.96*sstat_tab$SD/sqrt(1500)
sstat_tab$lower <- sstat_tab$Mean - 1.96*sstat_tab$SD/sqrt(1500)
sstat_tab$variable <- factor(rownames(sstat_tab),levels = rownames(sstat_tab))

sstat_tab3 <- data.frame(Mean=colMeans(sdata_pcak[1:500,]), SD=colSd(sdata_pcak[1:500,]))
sstat_tab3$upper <- sstat_tab3$Mean + 1.96*sstat_tab3$SD/sqrt(500)
sstat_tab3$lower <- sstat_tab3$Mean - 1.96*sstat_tab3$SD/sqrt(500)
sstat_tab3$variable <- factor(rownames(sstat_tab3),levels = rownames(sstat_tab3))

sstat_tab4 <- data.frame(Mean=colMeans(sdata_pcak[501:1000,]), SD=colSd(sdata_pcak[501:1000,]))
sstat_tab4$upper <- sstat_tab4$Mean + 1.96*sstat_tab4$SD/sqrt(500)
sstat_tab4$lower <- sstat_tab4$Mean - 1.96*sstat_tab4$SD/sqrt(500)
sstat_tab4$variable <- factor(rownames(sstat_tab4),levels = rownames(sstat_tab4))

sstat_tab5 <- data.frame(Mean=colMeans(sdata_pcak[1001:1500,]), SD=colSd(sdata_pcak[1001:1500,]))
sstat_tab5$upper <- sstat_tab5$Mean + 1.96*sstat_tab5$SD/sqrt(500)
sstat_tab5$lower <- sstat_tab5$Mean - 1.96*sstat_tab5$SD/sqrt(500)
sstat_tab5$variable <- factor(rownames(sstat_tab5),levels = rownames(sstat_tab5))

mydata <- rstat_tab

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-1.5, 4.5)

mydata <- rstat_tab3

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-5,17)

mydata <- rstat_tab4

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-2,5)


mydata <- rstat_tab5

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-5,3)


mydata <- sstat_tab

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-1, 4)


mydata <- sstat_tab3

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-5,17)

mydata <- sstat_tab4

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "")+
  ylim(-2,5)

mydata <- sstat_tab5

windows(width = 4, height = 4)
mydata %>%
  ggplot(aes(x = variable, y = Mean, group = 1)) +
  geom_polygon(aes(y = upper), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lower), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", linewidth = 2) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) + 
  coord_polar() +
  labs(x = "", y = "") +
  ylim(-5,3)
# library(ggradar2)
# library(ggplot2)
# rstat_tab_t <- data.frame(t(rstat_tab[,-5]))
# names(rstat_tab_t) <- rstat_tab[,5]
# v.mean <- rstat_tab_t[1,]
# row.names(v.mean) <- "ALL"
# ci.df <- rbind(rstat_tab_t[3,],rstat_tab_t[4,])
# ci.df$type = c("h","l")
# row.names(ci.df) <- c("ALL","ALL1")
# ggradar2::ggradar2(v.mean, ci= ci.df, gridline.label = seq(0, 100, 10), group.line.width=1, group.point.size=2)
