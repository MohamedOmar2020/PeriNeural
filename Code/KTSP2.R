###################################################################################
### Mohamed Omar
### 5/5/2019
### ### Goal: Creating the restricted ktsp classifier.
### TF_MIR genes
#################################################################################

###### 
# Clean Work space
rm(list = ls())
# Set work directory
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

## Load data
load("./Objs/PeriNeuralData.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup  


###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
featNo <- nrow(usedTrainMat)
ktsp <- c(1:25)

### Train a classifier using default filtering function based on Wilcoxon
ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo = 500)
ktspPredictorRes

### Check how many TSP
xx <- SWAP.Filter.Wilcoxon(usedTrainGroup,usedTrainMat,featureNo=featN)
dim(myTSPs[myTSPs[,1] %in% xx & myTSPs[,2] %in% xx ,])
choose(78,2)

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No", "Yes")), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No", "Yes"), ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
#bitmap("./Figs/TF_MiR/mechanistic_trainROC.png", width = 7.5, height = 7.5, res = 400)
roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = TRUE, print.thres=thr, print.auc=TRUE, print.auc.col="black", levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
#dev.off()

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Yes")


#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Threshold
#thr_test <- coords(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No", "Yes"),),"best")["threshold"]
#thr_test

## Print ROC curve local maximas
#coords(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No", "Yes"),),"local maximas")

## Plot curve
#bitmap("./Figs/TF_MiR/mechanistic_testROC_TF_MiRgenes.png", width = 7.5, height = 7.5, res = 400)
roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("No", "Yes"), col="blue", lwd=2, grid=TRUE)
#dev.off()

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Yes")

## OR (calculate accuracy, sens, spec and AUC)
SWAP.GetKTSP.PredictionStats(usedTestPredictionRes, usedTestGroup, decision_values = ktspStatsTestRes$statistics)

############################################################################
###########################################################################
### Plot genes in the training set

## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorRes$TSPs)

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
png("./Figs/TF_MiR/mechanistic.trainKTSPexprs.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTrain), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

########################################################################
#######################################################################
### Plot genes in the testing set

## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorRes$TSPs)

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
png("./Figs/TF_MiR/mechanistic.testKTSPexprs.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTest), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()
#########################################

## Make a heatmap of the TSPs expression in the test data
usedTestMat <- as.data.frame(usedTestMat)
TSP <- usedTestMat[which(rownames(usedTestMat)%in% ktspPredictorRes$TSPs),]
TSP <- as.matrix(TSP)
heatmap(TSP, scale = "row")

#OR
SWAP.PlotKTSP.Votes(ktspPredictorRes, as.matrix(usedTestMat), usedTestGroup)


##############################
TSPs <- ktspPredictorRes$TSPs
TSPs <- as.data.frame(TSPs)

GoodGns <- as.vector(TSPs$gene1)
BadGns <- as.vector(TSPs$gene2)

dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018" , "ChEA_2016" ,"KEGG_2019_Human")
Enriched_GoodGns <- enrichr(genes = GoodGns, databases = dbs)
Enriched_BadGns <- enrichr(genes = BadGns, databases = dbs)
printEnrich(EnrichedGoodGns, "PosGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
printEnrich(Enriched_BadGns, "NegGnsEnrichment.txt" , sep = "\t", columns = c(1:9))

Good_GO_BP <- Enriched_GoodGns["GO_Biological_Process_2018"]
Good_GO_BP <- Good_GO_BP$GO_Biological_Process_2018
Good_GO_BP <- Good_GO_BP[Good_GO_BP$P.value <= 0.05, ]

Bad_GO_BP <- Enriched_BadGns["GO_Biological_Process_2018"]
Bad_GO_BP <- Bad_GO_BP$GO_Biological_Process_2018
Bad_GO_BP <- Bad_GO_BP[Bad_GO_BP$P.value <= 0.05, ]


Good_KEGG <- Enriched_GoodGns["KEGG_2019_Human"]
Good_KEGG <- Good_KEGG$KEGG_2019_Human

Bad_KEGG <- Enriched_BadGns["KEGG_2019_Human"]
Bad_KEGG <- Bad_KEGG$KEGG_2019_Human

###########################


#######################################################################
#### SAVE
save(list = ls(pattern = "^ktsp"), file = "./Objs/TF_MiR/mechanistic.ktspPredictor.rda")

######
# html report
#rmarkdown::render(input = "./Code/TF_MiR/10KTSP_TF_MiR.R", output_dir = "./HTML", output_format = "html_document")