## Load the data
load("./Objs/PeriNeuralData.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


## transpose the matrix
predictor_data <- t(usedTrainMat)
predictor_names <- c(as.vector(rownames(usedTrainMat))) #gene symbol
colnames(predictor_data) <- predictor_names

#######################################################################

PCA_Train <- pca(predictor_data, scale = T, center = T)
plotIndiv(PCA_Train)
plotVar(PCA_Train)

PCA_Train <- spca(predictor_data, keepX = c(25,25))
plotIndiv(PCA_Train)
plotVar(PCA_Train, cex = 3)

plotIndiv(PCA_Train, group = usedTrainGroup, legend = TRUE)
plotIndiv(PCA_Train, group = usedTrainGroup, ind.names = FALSE, pch = as.factor(usedTrainGroup), legend = TRUE)

plotLoadings(PCA_Train, ndisplay = 50)

SPCA_Train <- spca(predictor_data, ncomp = 2, keepX = c(25,25))
plotIndiv(SPCA_Train, group = usedTrainGroup, legend = TRUE)
plotVar(SPCA_Train, cex = 3)
selectVar(SPCA_Train, comp = 1)$name


SPLSDA <- splsda(predictor_data, usedTrainGroup, keepX = c(25,25), scale = T)
plotVar(SPLSDA)

plotIndiv(SPLSDA , ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
plotVar(SPLSDA, cutoff = 0.5)

background <- background.predict(SPLSDA, comp.predicted = 2, dist = "max.dist")

auroc(SPLSDA)

##############################################################################

PLSDA_Train <- plsda(predictor_data, usedTrainGroup, ncomp = 4, scale = TRUE)
Perf_PLSDA <- perf(PLSDA_Train, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)

matplot(Perf_PLSDA$error.rate$BER, type = "l", lty = 1, col = color.mixo(1:3))
legend('topright', 
       c('max.dist', 'centroids.dist', 'mahalanobis.dist'), 
       lty = 1,
       col = color.mixo(1:3))

PLSDA_Train <- plsda(predictor_data, usedTrainGroup, ncomp = 4, scale = TRUE)

plotVar(PLSDA_Train, comp = c(3,4))

auroc(PLSDA_Train, roc.comp = 4)

Pred <- predict(PLSDA_Train, newdata = t(usedTestMat))
Confusion <- get.confusion_matrix(truth = usedTestGroup,predicted = Pred$MajorityVote$max.dist[,1])
Confusion
get.BER(Confusion)
#################################################################################
set.seed(333) 

list.keepX <- c(1:10,  seq(10, 100, 5))

Study <- list(Study1 = rep(2,44), Study2 = rep(3, 85), Study3 = rep(4, 545))
Study <- unlist(Study)
Study <- as.factor(Study)

tune.splsda.srbct <- tune.splsda(predictor_data, usedTrainGroup, ncomp = 4, dist = 'max.dist', progressBar = TRUE,
                                      measure = "BER", test.keepX = list.keepX,  auc = TRUE, scale = FALSE)  

error <- tune.splsda.srbct$error.rate
error
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
ncomp
select.keepX <- tune.splsda.srbct$choice.keepX[1:4]
select.keepX
plot(tune.splsda.srbct)

MyResult.splsda.final <- splsda(predictor_data, usedTrainGroup, ncomp = 4, keepX = select.keepX, scale = FALSE, max.iter = 200, near.zero.var = FALSE)


plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="SPLS-DA, Final result", study = "global")

plotVar(MyResult.splsda.final, comp = c(1,2), cex = 2, style = "ggplot2")


selectVar(MyResult.splsda.final, comp = 2)$name

auroc(MyResult.splsda.final, roc.comp = 1)

Pred <- predict(MyResult.splsda.final, newdata = t(usedTestMat))
Confusion <- get.confusion_matrix(truth = usedTestGroup,predicted = Pred$MajorityVote$max.dist[,3])
Confusion
get.BER(Confusion)


auroc(MyResult.splsda.final, newdata = t(usedTestMat), outcome.test = usedTestGroup,roc.comp = 1)

######################

