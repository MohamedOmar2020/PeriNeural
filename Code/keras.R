############################################################################
## Moahmed Omar
## 19/05/2019
## Goal: Creating a Deep Neural Network to predict the progression of Non-muscle invasive bladder cancer to Muscle-invasive 
## Using Keras/Tensor Flow
## Using all genes
#############################################################################

## Clean the Working Space
rm(list = ls())

## Setting the Working Directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Keras")

## Load Necessary Libraries
library(limma)
library(caret)
library(pROC)
library(DMwR)
library(genefilter)
library(keras)
#install_keras()
####################################

## Load data
load("./Objs/PeriNeuralData.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

## Filter out any variables (genes) that are not expressed or do not have enough variance to be informative in classification. 
## We will first take the values and un-log2 them, then filter out any genes according to following criteria: (1) At least 20% of samples should have raw intensity greater than 100; (2) The coefficient of variation (sd/mean) is between 0.7 and 10.
#X <- usedTrainMat
#ffun <- filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))

#filt <- genefilter(2^X,ffun)
#usedTrainMat <- usedTrainMat[filt, ]
#usedTestMat <- usedTestMat[filt, ]


### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
names_train <- c(as.vector(rownames(usedTrainMat)))
colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))

## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
names_Test <- c(as.vector(rownames(usedTestMat)))
colnames(Testing) <- names_Test
names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train$usedTrainGroup)  
Data_train$usedTrainGroup <- as.integer(Data_train$usedTrainGroup)

Train_label <- Data_train$usedTrainGroup-1
table(Train_label)


## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
Data_test$usedTestGroup <- as.integer(Data_test$usedTestGroup)
table(Data_test$usedTestGroup)
Test_label <- Data_test$usedTestGroup-1
table(Test_label)

################################################################
## Before making the NN we need to put all Data in a Matrix Format

Data_train <- as.matrix(Data_train)
Data_train[,1:4442] <- normalize(Data_train[,1:4442])
Data_train[,4443] <- as.numeric(Data_train[,4443])-1
table(Data_train[,4443])  ## 0: Progression, 1: NoProgression

## The same For the Test Data
Data_test <- as.matrix(Data_test)
Data_test[,1:4442] <- normalize(Data_test[,1:4442])
Data_test[,4443] <- as.numeric(Data_test[,4443])-1
table(Data_test[,4443])

Train_target <- Data_train[,4443]
Test_target <- Data_test[,4443]

## One Hot Encoding
TrainLabels <- to_categorical(Train_target) 
TestLabels <- to_categorical(Test_target)

## Performing Oversampling using SMOTE (to compensate for the unbalanced classes)
# set.seed(333)
# Data_train <- as.data.frame(Data_train)
# Data_train[,4997] <- as.factor(Data_train[,4997])
# Data_train <- SMOTE(usedTrainGroup~., Data_train, perc.over = 300, perc.under = 134)
# table(Data_train$usedTrainGroup)
# 
# Data_train <- as.matrix(Data_train)
# Data_train[,4997] <- as.numeric(Data_train[,4997])
# table(Data_train[,4997])
# 
# Train_target <- Data_train[,4997]
# TrainLabels <- to_categorical(Train_target) 

########################################################################
## Creating the Network

set.seed(333)
Model <- keras_model_sequential()
## Making 2Layers: the first layer is 200, the second is 2 (==The number of classes) 
Model %>% 
  layer_dense(units = 1000, activation = "relu", input_shape = c(4442)) %>%
  layer_dropout(0.2) %>%
  layer_dense(units = 600, activation = "relu") %>%
  layer_dense(units = 200, activation = "relu") %>%
  layer_dropout(0.2) %>%
  layer_batch_normalization() %>%
  layer_dense(units = 2, activation = "softmax")
summary(Model)

Model %>% 
  compile(loss="binary_crossentropy", optimizer= "adam", metrics ="accuracy")

History <- Model %>% fit(Data_train[,1:4442], TrainLabels, epoch=200, batch_size=32, validation_split = 0.2)

## Plot the model (epoch/iterations VS accuracy)
png("./Figs/Model_accuracy_allGnsFilt.png", width = 3000, height = 1500, res = 400)
plot(History)
dev.off()

## Evaluation of the NN Model on the Testing set
Model %>% evaluate(Data_test[,1:4442], TestLabels, batch_size=32) 

## Predict Probabilities in the Test Data
prob <- Model %>% predict_proba(Data_test[,1:4442], batch_size = 32)

## Predict classes in the Test Data
pred <- Model %>% predict_classes(Data_test[,1:4442], batch_size = 32)
table(pred)

## Create a confusion table 
table1 <- table(Predicted = pred, Actual=Test_target)
table1

## ROC and AUC in the Testing Data
png("./Figs/ROC_Test_allGnsFilt.png", width = 2000, height = 2000, res = 300)
roc(Test_target, pred, plot = TRUE, print.auc=TRUE, levels = c("1", "0"), col="blue", lwd=2, grid=TRUE, main="Keras ROC test using all genes (2913)")
dev.off()
