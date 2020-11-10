########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: Creating several machine learning models for Prostate cancer Metastasis
## Then Comparing the performance of these models 

###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")

## Load necessary packages
library(caret)
library(nnet)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(randomForest)
library(gbm)

#######################################################################

## Load data
load("./Objs/PeriNeuralData.rda")

## Load the RF important genes (ordered according to MeanDecreaseGini)
load("./Objs/RandomForest_Importance.rda")

## Keep only the top 500 genes
keepGns <- rownames(rf_importances)[1:300]

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
boxplot(usedTrainMat, outline = FALSE)

usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")
boxplot(usedTestMat, outline = FALSE)


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
names_train <- c(as.vector(rownames(usedTrainMat)))
colnames(Training) <- names_train


## Making sure that sample names are identical in both Training and usedTrainGroup
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
names_Test <- c(as.vector(rownames(usedTestMat)))
colnames(Testing) <- names_Test
all(rownames(Testing) == names(usedTestGroup))

Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

######
control <- trainControl(method="cv", number=10, classProbs = TRUE)


###########################################################################
###########################################################################
## Model1: C5.0

## Fit the model
set.seed(333)
fit.cart <- train(usedTrainGroup~., data=Data_train, method="C5.0", trControl=control)
fit.cart

## Predict classes
test_pred_classes_cart <- predict(fit.cart, Testing, type="raw")
table(test_pred_classes_cart)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_cart <- confusionMatrix(test_pred_classes_cart, usedTestGroup, positive = "Yes")
Confusion_test_cart

## ROC/AUC in the Testing set
test_pred_prob_cart <- predict(fit.cart, Testing, type = "prob")

png("./Figs/All_genes/ROC_test_cart.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_cart[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main="ROC Test Cart")
dev.off()

## Top 100 predictors
Importance_cart <- varImp(fit.cart, scale = FALSE)
Importance_cart <- Importance_cart$importance
Importance_cart <- as.matrix(Importance_cart)
Importance_cart <- Importance_cart[order(Importance_cart[,1], decreasing = TRUE),]
genes_cart <- names(Importance_cart)[1:100]
save(genes_cart, file = "./Objs/All_genes/genes_cart.rdata")

#######################################################################
########################################################################
## Model 2: GLMnet

set.seed(333)
fit.GLMnet <- train(usedTrainGroup~., data=Data_train, method="glmnet", trControl = control)
fit.GLMnet

## Predict classes
test_pred_classes_GLMnet <- predict(fit.GLMnet, Testing, type="raw")
table(test_pred_classes_GLMnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_GLMnet <- confusionMatrix(test_pred_classes_GLMnet, usedTestGroup, positive = "Yes")
Confusion_test_GLMnet

## ROC/AUC in the Testing set
test_pred_prob_GLMnet <- predict(fit.GLMnet, Testing, type = "prob")

png("./Figs/All_genes/ROC_test_GLMnet.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_GLMnet[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main="ROC Test GLMnet")
dev.off()

## Top 100 predictors
Importance_GLMnet <- varImp(fit.GLMnet, scale = FALSE)
plot(Importance_GLMnet, top = 30)
Importance_GLMnet <- Importance_GLMnet$importance
Importance_GLMnet <- as.matrix(Importance_GLMnet)
Importance_GLMnet <- Importance_GLMnet[order(Importance_GLMnet[,1], decreasing = TRUE),]
genes_GLMnet <- names(Importance_GLMnet)[1:100]
save(genes_GLMnet, file = "./Objs/All_genes/genes_GLMnet.rdata")

################################################################################
################################################################################
## Model 3: LDA

set.seed(333)
fit.lda <- train(usedTrainGroup~., data=Data_train, method="lda2", trControl=control)
fit.lda

test_pred_classes_lda <- predict(fit.lda, Testing, type="raw")
table(test_pred_classes_lda)
#test_pred_classes_lda <- as.vector(test_pred_classes_lda)
#names(test_pred_classes_lda) <- names_prediction
#test_pred_classes_lda <- ordered(test_pred_classes_lda, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_lda <- confusionMatrix(test_pred_classes_lda, usedTestGroup, positive = "Yes")
Confusion_test_lda

## Ranking features by importance
Importance_LDA <- varImp(fit.lda, scale = FALSE)
Importance_LDA <- Importance_LDA$importance
Importance_LDA <- Importance_LDA[order(Importance_LDA$Mets, decreasing = TRUE),]
genes_LDA <- rownames(Importance_LDA)[1:100]
save(genes_LDA, file = "./Objs/All_genes/genes_LDA.rdata")
## ROC/AUC in the Testing set
test_pred_prob_lda <- predict(fit.lda, Testing, type = "prob")

png("./Figs/All_genes/ROC_test_LDA.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_lda[,2], plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="ROC Test LDA")
dev.off()

#######################################################################
########################################################################
## Model 4: High Dimensional Discriminant Analysis (HDDA)		

set.seed(333)
fit.hdda <- train(usedTrainGroup~., data=Data_train, method="hdda", trControl=control)
fit.hdda

test_pred_classes_hdda <- predict(fit.hdda, Testing, type="raw")
table(test_pred_classes_hdda)
#test_pred_classes_lda <- as.vector(test_pred_classes_lda)
#names(test_pred_classes_lda) <- names_prediction
#test_pred_classes_lda <- ordered(test_pred_classes_lda, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_hdda <- confusionMatrix(test_pred_classes_hdda, usedTestGroup, positive = "Yes")
Confusion_test_hdda

## Ranking features by importance
Importance_HDDA <- varImp(fit.hdda, scale = FALSE)
Importance_HDDA <- Importance_HDDA$importance
Importance_HDDA <- Importance_HDDA[order(Importance_HDDA$Mets, decreasing = TRUE),]
genes_HDDA <- rownames(Importance_HDDA)[1:100]
save(genes_HDDA, file = "./Objs/All_genes/genes_HDDA.rdata")
## ROC/AUC in the Testing set
test_pred_prob_hdda <- predict(fit.hdda, Testing, type = "prob")

png("./Figs/All_genes/ROC_test_HDDA.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_hdda[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main="ROC Test HDDA")
dev.off()

#######################################################################
########################################################################
## Model 5: SVM Poly

set.seed(333)
fit.svmPoly <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=control)
fit.svmPoly

test_pred_classes_svmPoly <- predict(fit.svmPoly, Testing, type="raw")
table(test_pred_classes_svmPoly)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly <- confusionMatrix(test_pred_classes_svmPoly, usedTestGroup, positive = "Yes")
Confusion_test_svmPoly

## ROC/AUC in the Testing set
test_pred_prob_svmPoly <- predict(fit.svmPoly, Testing, type = "prob")
png("./Figs/All_genes/ROC_test_SVMPoly.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_svmPoly[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly")
dev.off()

## Top 100 predictors
Importance_SVMPoly <- varImp(fit.svmPoly, scale = FALSE)
plot(Importance_SVMPoly, top = 50)
Importance_SVMPoly <- Importance_SVMPoly$importance
Importance_SVMPoly <- Importance_SVMPoly[order(Importance_SVMPoly$Mets, decreasing = TRUE),]
genes_SVMPoly <- rownames(Importance_SVMPoly)[1:100]
save(genes_SVMPoly, file = "./Objs/All_genes/genes_SVMPloy.rdata")

#######################################################################
########################################################################
## Model 6: SVM Radial

set.seed(333)
fit.svmRadial <- train(usedTrainGroup~., data=Data_train, method="svmRadial", trControl=control)
fit.svmRadial

test_pred_classes_svmRadial <- predict(fit.svmRadial, Testing, type="raw")
table(test_pred_classes_svmRadial)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmRadial <- confusionMatrix(test_pred_classes_svmRadial, usedTestGroup, positive = "Yes")
Confusion_test_svmRadial

## ROC/AUC in the Testing set
test_pred_prob_svmRadial <- predict(fit.svmRadial, Testing, type = "prob")
png("./Figs/All_genes/ROC_test_SVMRadial.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_svmRadial[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Radial")
dev.off()

## Top 100 predictors
Importance_SVMRadial <- varImp(fit.svmRadial, scale = FALSE)
Importance_SVMRadial <- Importance_SVMRadial$importance
Importance_SVMRadial <- Importance_SVMRadial[order(Importance_SVMRadial$Mets, decreasing = TRUE),]
genes_SVMRadial <- rownames(Importance_SVMRadial)[1:100]
save(genes_SVMRadial, file = "./Objs/All_genes/genes_SVMRadial.rdata")

#######################################################################
########################################################################
## Model 7: Linear Distance Weighted Discrimination		

set.seed(333)
fit.dwd <- train(usedTrainGroup~., data=Data_train, method="dwdLinear", trControl=control)
fit.dwd

test_pred_classes_dwd	<- predict(fit.dwd, Testing, type="raw")
table(test_pred_classes_dwd)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_dwd <- confusionMatrix(test_pred_classes_dwd, usedTestGroup, positive = "Yes")
Confusion_test_dwd

## ROC/AUC in the Testing set
test_pred_prob_dwd <- predict(fit.dwd, Testing, type = "prob")
png("./Figs/All_genes/ROC_test_dwd.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_dwd[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main = "ROC Test DWD")
dev.off()

## Top 100 predictors
Importance_DWD <- varImp(fit.dwd, scale = FALSE)
Importance_DWD <- Importance_DWD$importance
Importance_DWD <- Importance_DWD[order(Importance_DWD$Mets, decreasing = TRUE),]
genes_DWD <- rownames(Importance_DWD)[1:100]
save(genes_DWD, file = "./Objs/All_genes/genes_DWD.rdata")

#######################################################################
########################################################################
## Model 8: Naïve Bayes

set.seed(333)
fit.nb <- train(usedTrainGroup~., data=Over_Train, method="nb", trControl=control)
fit.nb

test_pred_classes_nb <- predict(fit.nb, Testing, type="raw")
table(test_pred_classes_nb)
#test_pred_classes_rf <- as.vector(test_pred_classes_rf)
#names(test_pred_classes_rf) <- names_prediction
#test_pred_classes_rf <- ordered(test_pred_classes_rf, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_nb <- confusionMatrix(test_pred_classes_nb, usedTestGroup, positive = "Mets")
Confusion_test_nb

## ROC/AUC in the Testing set
test_pred_prob_nb <- predict(fit.nb, Testing, type = "prob")

png("./Figs/All_genes/ROC_test_NB.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_nb[,2], plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "ROC Test Naive Bayes")
dev.off()

## Top 100 predictors
Importance_NB <- varImp(fit.nb, scale = FALSE)
Importance_NB <- Importance_NB$importance
Importance_NB <- Importance_NB[order(Importance_NB$Mets, decreasing = TRUE), ]
genes_NB <- rownames(Importance_NB)[1:100]
save(genes_NB, file = "./Objs/All_genes/genes_NB.rda")

################################################################################
## Model 9: Neural Network

set.seed(333)
fit.nnet <- train(usedTrainGroup~., data=Data_train, method="nnet", trControl=control)
fit.nnet

test_pred_classes_nnet <- predict(fit.nnet, Testing, type="raw")
table(test_pred_classes_nnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_nnet <- confusionMatrix(test_pred_classes_nnet, usedTestGroup, positive = "Yes")
Confusion_test_nnet

## ROC/AUC in the Testing set
test_pred_prob_nnet <- predict(fit.nnet, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_NNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_nnet[,1], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

## Top 100 predictors
Importance_nnet <- varImp(fit.nnet, scale = FALSE)
Importance_nnet <- Importance_nnet$importance
Importance_nnet <- as.matrix(Importance_nnet)
Importance_nnet <- Importance_nnet[order(Importance_nnet[,1], decreasing = TRUE),]
genes_nnet <- names(Importance_nnet)[1:100]
save(genes_nnet, file = "./Objs/genes_nnet.rda")

###########################################################################
#########################################################################
## Model 10: Model Averaged Neural Network	

set.seed(333)
fit.AVnnet <- train(usedTrainGroup~., data=Data_train, method="avNNet", trControl=control)
fit.AVnnet

test_pred_classes_AVnnet <- predict(fit.AVnnet, Testing, type="raw")
table(test_pred_classes_AVnnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_AVnnet <- confusionMatrix(test_pred_classes_AVnnet, usedTestGroup, positive = "Mets")
Confusion_test_AVnnet

## ROC/AUC in the Testing set
test_pred_prob_AVnnet <- predict(fit.AVnnet, Testing, type = "prob")

png("./Figs/ROC_test_AVNNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_AVnnet[,1], plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

## Top 100 predictors
Importance_AVnnet <- varImp(fit.AVnnet, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_AVnnet, top = 20)
Importance_AVnnet <- Importance_AVnnet$importance
Importance_AVnnet <- Importance_AVnnet[order(Importance_AVnnet$Mets, decreasing = TRUE),]
genes_AVnnet <- rownames(Importance_AVnnet)[1:100]
save(genes_AVnnet, file = "./Objs/genes_AVnnet.rda")

#############################################################################
#############################################################################
## Model 11: Neural Networks with Feature Extraction	

set.seed(333)
fit.PCAnnet <- train(usedTrainGroup~., data=Data_train, method="pcaNNet", trControl=control)
fit.PCAnnet

test_pred_classes_PCAnnet <- predict(fit.PCAnnet, Testing, type="raw")
table(test_pred_classes_PCAnnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_PCAnnet <- confusionMatrix(test_pred_classes_PCAnnet, usedTestGroup, positive = "Mets")
Confusion_test_PCAnnet

## ROC/AUC in the Testing set
test_pred_prob_PCAnnet <- predict(fit.PCAnnet, Testing, type = "prob")

png("./Figs/ROC_test_PCA_NNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_PCAnnet[,2], plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

## Top 100 predictors
Importance_PCAnnet <- varImp(fit.PCAnnet, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_PCAnnet, top = 20)
Importance_PCAnnet <- Importance_PCAnnet$importance
Importance_PCAnnet <- Importance_PCAnnet[order(Importance_PCAnnet$Mets, decreasing = TRUE),]
genes_PCAnnet <- rownames(Importance_PCAnnet)[1:100]
save(genes_PCAnnet, file = "./Objs/genes_PCAnnet.rda")

#############################################################################
#############################################################################
## Model 12: Random Forest

usedTrainGroup <- Data_train$usedTrainGroup
Data_train$usedTrainGroup <- NULL
PredictorData <- Data_train

set.seed(333)
fit.rf <- train(x = PredictorData, y = usedTrainGroup, method="rf", trControl=control, tuneLength = 5)
fit.rf

test_pred_classes_rf <- predict(fit.rf, Testing, type="raw")
table(test_pred_classes_rf)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_rf <- confusionMatrix(test_pred_classes_rf, usedTestGroup, positive = "Mets")
Confusion_test_rf

## ROC/AUC in the Testing set
test_pred_prob_rf <- predict(fit.rf, Testing, type = "prob")

png("./Figs/All_genes/ROC_test_RF.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_rf[,2], plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Random Forest")
dev.off()

Importance_RF <- varImp(fit.rf, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_RF, top = 50)
Importance_RF <- Importance_RF$importance
Importance_RF <- as.matrix(Importance_RF)
Importance_RF <- Importance_RF[order(Importance_RF[,1], decreasing = TRUE),]
genes_RF <- names(Importance_RF)[1:100]
save(genes_RF, file = "./Objs/All_genes/genes_RF.rdata")

###########################################################################
## Model 14: Stochastic Gradient Boosting
control$sampling <- "up"

set.seed(333)
fit.gbm <- train(x = t(usedTrainMat), y= usedTrainGroup, method="gbm", trControl=control)
fit.gbm
Model_Performance <- gbm.perf(fit.gbm$finalModel)

test_pred_classes_gbm <- predict(fit.gbm, t(usedTestMat), type="raw")
table(test_pred_classes_gbm)


## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_gbm <- confusionMatrix(test_pred_classes_gbm, usedTestGroup, positive = "Yes")
Confusion_test_gbm

## ROC/AUC in the Testing set
test_pred_prob_gbm <- predict(fit.gbm, t(usedTestMat), type = "prob")

png("./Figs/All_genes/ROC_test_gbm.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_gbm[,2], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main= "ROC Test GBM")
dev.off()

Importance_gbm <- varImp(fit.gbm, scale = TRUE)
## Plot the top 20 predictors
plot(Importance_gbm, top = 50)
Importance_gbm <- Importance_gbm$importance
Importance_gbm <- as.matrix(Importance_gbm)
Importance_gbm <- Importance_gbm[order(Importance_gbm[,1], decreasing = TRUE),]
genes_gbm <- names(Importance_gbm)[1:500]
save(genes_gbm, file = "./Objs/genes_gbm.rdata")

################################################################################
### Comparison between the models

# collect resamples
results <- resamples(list(CART=fit.cart, GLM_net = fit.GLMnet, LDA=fit.lda, HDDA = fit.hdda, SVM_Poly=fit.svmPoly, SVM_Radial = fit.svmRadial, LDWD = fit.dwd, NB=fit.nb, NNET=fit.nnet, NNET_AV = fit.AVnnet, NNET_PCA = fit.PCAnnet, RF=fit.rf, GBM=fit.gbm))
summary(results)

# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/Prognosis_genes/BWplot.png", width = 3000, height = 3000, res = 300)
bwplot(results, scales=scales)
dev.off()

# density plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/Prognosis_genes/DensPlot.png", width = 3000, height = 3000, res = 300)
densityplot(results, scales=scales, pch = "|")
dev.off()

# dot plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/Prognosis_genes/DotPlot.png", width = 3000, height = 3000, res = 300)
dotplot(results, scales=scales)
dev.off()

# pair-wise scatterplots of predictions to compare models
png("./Figs/Prognosis_genes/SPlom.png", width = 3000, height = 3000, res = 300)
splom(results)
dev.off()

# difference in model predictions
diffs <- diff(results)
# summarize p-values for pair-wise comparisons
summary(diffs)

