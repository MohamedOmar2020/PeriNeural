#################################################################################
### Mohamed Omar
### 2/08/2019
### GOAL: Creating a rondom forest classifier for Perineural invasion
### 
#################################################################################

# Clean the work space
rm(list = ls())

## settng the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")

## Load necessary libraries
library(randomForest)
library(pROC)
library(ROCR)
library(caret)
library(limma)
library(genefilter)


## Load the data
load("./Objs/PeriNeuralData.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

all(names(usedTrainGroup) == colnames(usedTrainMat))
all(names(usedTestGroup) == colnames(usedTestMat))
################################################################################
################################################################################
################################################################################
### Creating the classifier using the training set

## transpose the matrix
predictor_data <- t(usedTrainMat)
predictor_names <- c(as.vector(rownames(usedTrainMat))) #gene symbol
colnames(predictor_data) <- predictor_names

Data_Train <- cbind(t(usedTrainMat), usedTrainGroup)
## Setting the variable we are trying to predict as our target variable. In this case, it is Progression status.
## train group here is just the column containing the phenptype of interest (Progression vs NoProgression) from the phenotype table

#usedTrainGroup <- ordered(usedTrainGroup, levels=c("NoProgression", "Progression"))
target <- usedTrainGroup

## Finally we run the RF algorithm. 
## NOTE: use an ODD number for ntree. This is because when the forest is used on test data, ties are broken randomly. Having an odd number of trees avoids this issue.
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

######


# Tunning RF model (to find the best mtry)
set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 1000, proximity = TRUE, sampsize = sampsizes)  


## Train the Random Forest Model
set.seed(333)
rf_output <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 1000, proximity=TRUE, na.action = na.omit, sampsize = sampsizes, strata = usedTrainGroup, nodesize = 500, maxnodes = 100)  

####

## Print a summary of the Random Forest Model
print(rf_output)

##
plot(rf_output)

## Importance measures for each variable.
rf_importances <- importance(rf_output, scale=FALSE)

rf_importances <- rf_importances[order(rf_importances[,"MeanDecreaseGini"], decreasing = TRUE), ]
rf_importances <- rf_importances[rf_importances[,"MeanDecreaseGini"] != 0, ]

######
## Predict in the training data
Prediction_Train <- predict(rf_output, predictor_data)
confusion_train <- confusionMatrix(Prediction_Train, usedTrainGroup, positive = "Yes")
confusion_train

## Create a representation of the top 50 variables categorized by importance.
png(filename = "./Figs/RandomForest/VarImp.png", width = 2000, height = 2000, res = 300)
varImpPlot(rf_output, type=2, n.var=30, sort = TRUE, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

## An MDS plot provides a sense of the separation of classes.
png(filename="./Figs/RandomForest/MDS_plot.png", width = 2000, height = 2000, res = 300)
target_labels=as.vector(usedTrainGroup)
MDSplot(rf_output, usedTrainGroup, k=2, pch=target_labels, palette=c("red", "blue"), main="MDS plot")
dev.off()



train_pred <- predict(rf_output, newdata = predictor_data, type = "vote")
roc(usedTrainGroup, train_pred[,1], plot = TRUE, print.auc=TRUE, levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, auc = TRUE, ci = TRUE, percent = TRUE)


################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(usedTestMat)
predictor_names2 <- c(as.vector(rownames(usedTestMat))) #gene symbol
colnames(predictor_data2) <- predictor_names2


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(rf_output$importance)
predictor_data2 <- predictor_data2[,colnames(predictor_data)]

## Run the test data through forest!
RF_predictions_responses <- predict(rf_output, predictor_data2, type="response")
RF_predictions_votes <- predict(rf_output, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses, usedTestGroup, positive = "Yes")
confusion_test

## ROC curve and AUC
png(filename = "./Figs/ROC_Test.png", width = 2000, height = 2000, res = 300)
roc(usedTestGroup, RF_predictions_votes[,1], plot = TRUE, print.auc = TRUE, levels = c("No", "Yes"), col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)
dev.off()

## Save RF classifier
save(rf_output, file = "./Objs/RandomForestClassifier.rda")

## Save RF importances
save(rf_importances, file = "./Objs/RandomForest_Importance.rda")

