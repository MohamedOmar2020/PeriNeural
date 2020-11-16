#################################################################################
## Mohamed Omar
## 07/08/2019
## Goal: Creating a predictive model for perineural invasion in cancer using XGB
##############################################################################


rm(list = ls())

#setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")

############################################################################

## Load necessary packages
library(xgboost)
library(pdp)
library(lime)
library(preprocessCore)
library(limma)
library(pROC)
library(caret)
library(DiagrammeR)
library(ggplot2)
library(xgboostExplainer)
library(dplyr)
library(Ckmeans.1d.dp)
library(mltools)

## Load data
load("./Objs/PeriNeuralData.rda")

############################################################################
## Load the selected genes *(TF-MiR from Lotte)
AdhesionGns <- read.delim("/Volumes/Macintosh/MechPred/MechPred/USER/Mohamed/MechanisticModels/Genes/GO_Adhesion.txt")
AdhesionGns <- AdhesionGns$GO_BIOLOGICAL_ADHESION

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Common genes
keepGns <- intersect(as.vector(AdhesionGns), rownames(usedTrainMat))
usedTrainMat <- usedTrainMat[keepGns, ]
usedTestMat <- usedTestMat[keepGns, ]

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

#####################
## Here, WE divide the training data into "actual training" and another validation
# This validation data will be used in the "WatchList" parameter. It is independent of the testing data.
set.seed(333)
inds <- createDataPartition(usedTrainGroup, p = 0.7, list = F)
Training1 <- Training[inds, ]
Validation <- Training[-inds, ]

usedTrainGroup1 <- usedTrainGroup[inds]
usedValGroup <- usedTrainGroup[-inds]

table(usedTrainGroup1)
table(usedValGroup)


## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training1) == names(usedTrainGroup1))
all(rownames(Validation) == names(usedValGroup))

## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training1)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train <- cbind(Training, usedTrainGroup1)

## The same for validation
Validation <- as.data.frame(Validation)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_val <- cbind(Validation, usedValGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train$usedTrainGroup1)  
levels(Data_train$usedTrainGroup1) <- c(1,0) 
Data_train$usedTrainGroup1 <- factor(Data_train$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
Train_label <- Data_train$usedTrainGroup1
Train_label <- as.vector(Train_label)
table(Train_label)

## The same for validation
table(Data_val$usedValGroup)  
levels(Data_val$usedValGroup) <- c(1,0) 
Data_val$usedValGroup <- factor(Data_val$usedValGroup, levels = c(0,1)) # 0=No,1= Yes 
Val_label <- Data_val$usedValGroup
Val_label <- as.vector(Val_label)
table(Val_label)

## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_test$usedTestGroup)  
levels(Data_test$usedTestGroup) <- c(1,0)
Data_test$usedTestGroup <- factor(Data_test$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
Test_label <- Data_test$usedTestGroup
Test_label <- as.vector(Test_label)
table(Test_label)

## Convert to xgb.DMatrix
DataTrain <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataVal <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
DataTest <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain, test = DataVal)



## Make a 5-fold CV model to determine the best number of trees/iterations
# hyper_grid <- expand.grid(
#   eta = 0.3,
#   max_depth = 1,
#   min_child_weight = 1,
#   subsample = 1,
#   colsample_bytree = 1,
#   gamma = c(0,1,2),
#   lambda = c(0,0.3,0.5, 1),
#   alpha = c(0, 0.3, 0.5, 1),
#   optimal_trees = 0,               # a place to dump results
#   max_AUC = 0                     # a place to dump results
# )
# 
# # grid search
# for(i in 1:nrow(hyper_grid)) {
# 
#   # create parameter list
#   params <- list(
#     eta = hyper_grid$eta[i],
#     max_depth = hyper_grid$max_depth[i],
#     min_child_weight = hyper_grid$min_child_weight[i],
#     subsample = hyper_grid$subsample[i],
#     colsample_bytree = hyper_grid$colsample_bytree[i],
#     gamma = hyper_grid$gamma[i],
#     lambda = hyper_grid$lambda[i],
#     alpha = hyper_grid$alpha[i]
#   )
# 
#   # reproducibility
#   set.seed(1234)
# 
#   # train model
#   xgb.tune <- xgb.cv(
#     params = params,
#     data = DataTrain,
#     nrounds = 500,
#     nfold = 5,
#     objective = "binary:logistic",
#     eval_metric        = "auc",# for regression models
#     verbose = 1,               # silent,
#     early_stopping_rounds = 50, 
#     scale_pos_weight = No/Yes
#   )
# 
#   # add max training AUC and trees to grid
#   hyper_grid$optimal_trees[i] <- which.max(xgb.tune$evaluation_log$test_auc_mean)
#   hyper_grid$max_AUC[i] <- max(xgb.tune$evaluation_log$test_auc_mean)
# }
# 
# hyper_grid %>% arrange(max_AUC) %>% head(100)
# write.csv(hyper_grid, file = "./Objs/hyper_grid.csv")

##########################
Yes <- sum(Train_label == 1)
No <- sum(Train_label == 0)



## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #1   # default = 0,   range: [0,∞]
  max_depth          = 2,             # 1
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,           # 0.4      # default = 1,   range: (0,1]
  colsample_bytree   = 1,                 # default = 1,   range: (0,1]
  colsample_bylevel  = 1,                 # default = 1,   range: (0,1]
  lambda             = 1,             # 0   # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

## Make the final model  1411  507
xgb.model <- xgb.train(parameters, DataTrain, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = No/Yes)

################################################

## Plot the iterations VS the AUC
d <- xgb.model$evaluation_log
n <- nrow(d)
v <- xgb.model$best_iteration
df <- data.frame(x=rep(d$iter, 2), val=c(d$train_auc, d$test_auc), 
                 set=rep(c("train", "test"), each=n))
ggplot(data = df, aes(x=x, y=val)) + 
  geom_line(aes(colour=set)) + 
  geom_vline(xintercept=v) + 
  theme_bw() +
  labs(title="AUC values for XGBoost", x="Iteration", y="AUC values (train, test)")

##################################################3
################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train <- predict(xgb.model, DataTrain)
ROC_Train <- roc(Train_label, XGB_prob_Train, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")
ROC_Train

thr_Train <- coords(roc(Train_label, XGB_prob_Train, levels = c("0", "1"), direction = "<"), transpose = TRUE, "best")["threshold"]
thr_Train

## Convert predicted probabilities to binary outcome
prediction_Train <- as.numeric(XGB_prob_Train > thr_Train)
print(head(prediction_Train))

Train_label <- factor(Train_label, levels = c(0,1))
prediction_Train<- factor(prediction_Train, levels = c(0,1))

# Confusion matrix in the training data
CM_Train <- confusionMatrix(prediction_Train, Train_label, positive = "1")
CM_Train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = prediction_Train, actuals = Train_label)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROC_Train$ci, CM_Train$overall["Accuracy"], CM_Train$byClass["Balanced Accuracy"], CM_Train$byClass["Sensitivity"], CM_Train$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#################################################

## Predict in the Test data
xgb_prob_test <- predict(xgb.model, DataTest)

## Convert predicted probabilities to binary outcome
prediction_Test <- as.numeric(xgb_prob_test > thr_Train)
print(head(prediction_Test))

Test_label <- factor(Test_label, levels = c(0,1))
prediction_Test <- factor(prediction_Test, levels = c(0,1))

## Confusion matrix
CM_Test <- confusionMatrix(prediction_Test, Test_label, positive = "1")
CM_Test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = prediction_Test, actuals = Test_label)
MCC_Test

## ROC curve and AUC
ROCTest <- roc(Test_label, xgb_prob_test, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the testing cohort (Agnostic)")
ROCTest

## Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, CM_Test$overall["Accuracy"], CM_Test$byClass["Balanced Accuracy"], CM_Test$byClass["Sensitivity"], CM_Test$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
XGB_Performance <- cbind(TrainPerf, TestPerf)

# Save
#save(XGB_Performance, file = "./Objs/XGB_Performance")


### basic plots of various metrics
#library(precrec)
#precrec_obj <- evalmod(scores = xgb.predict_prob, labels = Test_label, mode="basic")
#autoplot(precrec_obj)   


##############################################
## Classification error rate
# err <- mean(as.numeric(xgb.predict_prob > thr) != Test_label)
# err
# 
# ##############
# ## Plot the prediction results
# 
# results <- data.frame(actual = Test_label, Prediction = prediction)
# 
# results$correct <- ifelse(results$actual == results$Prediction, TRUE, FALSE)
# 
# ggplot(results, aes(x = Prediction, fill = correct)) +
#   geom_bar(position = "dodge")
# 
# ## Confusion matrix
# ggplot(results, aes(x = Prediction, y = Test_label, colour = correct)) +
#   geom_jitter()
# 
# plot(xgb.model$evaluation_log)
# 
# png(filename = "./Figs/IterVsAUC.png", width = 2000, height = 1500, res = 300)
# GG <- ggplot(data = xgb.model$evaluation_log, aes(x = iter, y = test_auc)) + geom_jitter(size = 2, alpha = 0.6)
# GG + geom_vline(xintercept=xgb.model$best_iteration, color="black", linetype="dashed")
# dev.off()
# 
# ##############
# 
# # create importance matrix
# importance_matrix <- xgb.importance(model = xgb.model)
# write.csv(importance_matrix, file = "./Objs/XGB_importanceMatrix.csv")
# Keep <- importance_matrix$Feature
# save(Keep, file = "./Objs/Keep.rda")
# 
# # variable importance plot
# importance_matrix <- importance_matrix[order(importance_matrix$Gain, decreasing = TRUE), ]
# 
# png(filename = "./Figs/XGB_Importance.png", width = 3000, height = 2000, res = 300)
# ggplot(data=importance_matrix[1:20,], aes(x = reorder(Feature, -Gain), y = Gain)) +
#   geom_bar(stat="identity", colour = "black", fill = "lightgray") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "XG Boosted Feature Importance (Top 20)", x = "Features", y = "Information Gain")
# dev.off()
# 
# ####################
# xgb.plot.multi.trees(model = xgb.model, features_keep = 10)
# xgb.plot.deepness(model = xgb.model, which='max.depth', pch=16, col=rgb(0,0,1,0.3), cex=2)
# ##############
# 
# #### The XGBoost Explainer
# explainer = buildExplainer(xgb.model ,DataTrain, type="binary", base_score = 0.5, trees = NULL)
# pred.breakdown = explainPredictions(xgb.model, explainer, DataTest)
# 
# showWaterfall(xgb.model, explainer, DataTest, as.matrix(Testing),  2, type = "binary")
# showWaterfall(xgb.model, explainer, DataTest, as.matrix(Testing),  8, type = "binary")
# 


###############################################################################
## Gene set enrichment analysis
# library(enrichR)
# dbs <- listEnrichrDbs()
# dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018" , "ChEA_2016" ,"KEGG_2019_Human")
# Enriched_Genes <- enrichr(genes = importance_matrix$Feature, databases = dbs)
# printEnrich(Enriched_Genes, "XGB_gene_enrichment.txt" , sep = "\t", columns = c(1:9))
# 
# XGB_GO_BP <- Enriched_Genes["GO_Biological_Process_2018"]
# XGB_GO_BP <- XGB_GO_BP$GO_Biological_Process_2018
# XGB_GO_BP <- XGB_GO_BP[XGB_GO_BP$P.value <= 0.05, ]
# write.csv(XGB_GO_BP, file = "./Objs/GO.csv")
# 
# XGB_KEGG <- Enriched_Genes["KEGG_2019_Human"]
# XGB_KEGG <- XGB_KEGG$KEGG_2019_Human
# XGB_KEGG <- XGB_KEGG[XGB_KEGG$P.value <= 0.05, ]
# write.csv(XGB_KEGG, file = "./Objs/KEGG.csv")
# 
# ## Bar plot of KEGG Pathways
# png(filename = "./Figs/KeggPathways.png", width = 2500, height = 2000, res = 300)
# ggplot(data=XGB_KEGG, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
#   geom_bar(stat="identity", colour = "black", fill = "lightgray") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Enriched pathways", x = "Pathway", y = "Score")
# dev.off()