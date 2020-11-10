###############################################################################
### Mohamed Omar
### 14/07/2019
### Goal : Creating the Agnostic K-TSP classifier for prostate cancer metastasis
#################################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)

###################################################################
## Load Metastasis Data
load("./Objs/PeriNeuralData.rda")

## Load the RF important genes (ordered according to MeanDecreaseGini)
load("./Objs/Genes_XGB.rda")

## Keep only the top 500 genes
keepGns <- Genes_XGB #300

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(mixTrainMat)[keepGns, ]
boxplot(usedTrainMat, outline = FALSE)

usedTestMat <- normalizeBetweenArrays(mixTestMat)
boxplot(usedTestMat, outline = FALSE)

#
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Train a classifier using the default filter function
ktspPredictorUnRes <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, krange = c(1:25), FilterFunc = NULL)
ktspPredictorUnRes

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainUnRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTrainUnRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("Yes", "No"), ), transpose = TRUE, "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("Yes", "No"), ), transpose = TRUE ,"local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, plot = TRUE, print.thres=thr, print.thres.adj=c(0.01,1.25), print.auc=TRUE, print.auc.col="black", levels = c("Yes", "No"), col="blue", lwd=2, grid=TRUE, main="Agnostic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionUnRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorUnRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionUnRes, usedTrainGroup, positive = "Yes")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestUnRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTestUnRes$statistics)

## Threshold
thr_test <- coords(roc(usedTestGroup, ktspStatsTestUnRes$statistics, levels = c("Yes", "No"),),"best", transpose = FALSE)["threshold"]
thr_test

## Print ROC curve local maximas
coords(roc(usedTestGroup, ktspStatsTestUnRes$statistics, levels = c("Yes", "No"),), transpose = TRUE, "local maximas")

## Plot curve
png("./Figs/Combined_KTSP_RF/ROC_Test.png", width = 2000, height = 2000, res = 300)
roc(usedTestGroup, ktspStatsTestUnRes$statistics, plot = TRUE, print.thres=thr_test, print.auc=TRUE, print.auc.col="black", levels = c("Yes", "No"), col="blue", lwd=2, grid=TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main= "KTSP performance in the testing set")
dev.off()

### Get predictions based on best threshold from ROC curve
usedTestPredictionUnRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorUnRes, DecisionFunc = function(x) sum(x) > thr_test)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionUnRes, usedTestGroup, positive = "Yes")

## OR (calculate accuracy, sens, spec and AUC)
SWAP.GetKTSP.PredictionStats(usedTestPredictionUnRes, usedTestGroup, decision_values = ktspStatsTestUnRes$statistics)

############################################################################
###########################################################################
### Plot genes in the training set

## Which TSPs
i <- 1:nrow(ktspPredictorUnRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorUnRes$TSPs)

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=usedTrainMat, g=usedTrainGroup)

## Reduce
datTrain <- Reduce("rbind", dfTspTrain)

## Rename columns
colnames(datTrain)[colnames(datTrain) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
png("./Figs/Combined_KTSP_RF/trainKTSPexprs.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTrain), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Group), position = position_jitterdodge(1.2), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

########################################################################
#######################################################################
### Plot genes in the testing set

## Which TSPs
i <- 1:nrow(ktspPredictorUnRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorUnRes$TSPs)

## Assemble
dfTspTest <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=usedTestMat, g=usedTestGroup)

## Reduce
datTest <- Reduce("rbind", dfTspTest)

## Rename columns
colnames(datTest)[colnames(datTest) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
png("./Figs/Combined_KTSP_RF/testKTSPexprs.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTest), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Group), position = position_jitterdodge(1.2), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

#######################################################################
#### SAVE
save(list = ls(pattern = "^ktsp"), file = "./Objs/Combined_KTSP_RF/Combined_ktspPredictor.rda")


#################################################################################
##################################################################################
##################################################################################

