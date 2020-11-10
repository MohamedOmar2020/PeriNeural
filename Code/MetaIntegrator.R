Dataset2$expr <- expr2
Dataset2$pheno <- Pheno2
Dataset2$keys <- rownames(expr2)

Dataset4$expr <- expr4
Dataset4$pheno <- Pheno4
Dataset4$keys <- rownames(expr4)

Dataset5 <- list()
Dataset5$expr <- expr5
Dataset5$pheno <- Pheno5
Dataset5$keys <- rownames(expr5)

Dataset6 <- list()
Dataset6$expr <- expr6
Dataset6$pheno <- Pheno6
Dataset6$keys <- rownames(expr6)



########
PerineuralData <- list()
PerineuralData$originalData <- list(Dataset2, Dataset4, Dataset5, Dataset6)
names(PerineuralData$originalData) <- c("GSE103479", "GSE7055", "hnsc_tcga", "chol_tcga")

PerineuralData$originalData$GSE103479 <- classFunction(PerineuralData$originalData$GSE103479, column = "Perineural", diseaseTerms = c("Yes"))
PerineuralData$originalData$GSE7055 <- classFunction(PerineuralData$originalData$GSE7055, column = "Perineural", diseaseTerms = c("Yes"))
PerineuralData$originalData$hnsc_tcga <- classFunction(PerineuralData$originalData$hnsc_tcga, column = "Perineural", diseaseTerms = c("Yes"))
PerineuralData$originalData$chol_tcga <- classFunction(PerineuralData$originalData$chol_tcga, column = "Perineural", diseaseTerms = c("Yes"))

PerineuralMeta <- PerineuralData

PerineuralMeta$originalData$hnsc_tcga$formattedName <- "hnsc_tcga"
PerineuralMeta$originalData$chol_tcga$formattedName <- "chol_tcga"

PerineuralMeta$originalData$hnsc_tcga$expr <- as.matrix(PerineuralMeta$originalData$hnsc_tcga$expr)
PerineuralMeta$originalData$chol_tcga$expr <- as.matrix(PerineuralMeta$originalData$chol_tcga$expr)

checkDataObject(PerineuralMeta, "Meta", "Pre-Analysis")

## Run the meta analysis
Perineural_Metaanalysis <- runMetaAnalysis(PerineuralMeta, runLeaveOneOutAnalysis = F, maxCores = 3)

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
Perineural_Metaanalysis <- filterGenes(Perineural_Metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.1, numberStudiesThresh = 4, heterogeneityPvalThresh = 0.05)

## Assigning a name to the filter
filter <- Perineural_Metaanalysis$filterResults[[1]]

## Summarize filter results
filter_summary <- summarizeFilterResults(metaObject = Perineural_Metaanalysis, getMostRecentFilter(Perineural_Metaanalysis))

## Modify the gene signature for more accuracy and AUC
# Using forward search 
New_filter <- forwardSearch(metaObject = Perineural_Metaanalysis, filterObject = filter)

## Replace the old filter with the new smaller one
Perineural_Metaanalysis$filterResults$FDR0.1_es0_nStudies4_looaFALSE_hetero0.05$posGeneNames <- New_filter$posGeneNames
Perineural_Metaanalysis$filterResults$FDR0.1_es0_nStudies4_looaFALSE_hetero0.05$negGeneNames <- New_filter$negGeneNames

New_filter <- Perineural_Metaanalysis$filterResults[[1]]
New_filter_summary <- summarizeFilterResults(metaObject = Perineural_Metaanalysis, getMostRecentFilter(Perineural_Metaanalysis))

#####################################

Dataset7 <- list()
Dataset7$expr <- expr7
Dataset7$pheno <- Pheno7
Dataset7$keys <- rownames(expr7)
Dataset7$expr <- as.matrix(Dataset7$expr)
Dataset7$formattedName <- "coadread_tcga"

## Label the samples
Dataset7 <- classFunction(Dataset7, column = "Perineural", diseaseTerms = c("Yes"))

rocPlot(datasetObject = Dataset7, filterObject = New_filter)

