###############################################################################
### Mohamed Omar
### 02/08/2019
### Goal : Collecting Perineural invasion data
#################################################################################

rm(list = ls())

#setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")


############################################################################
### Load library
require(GEOquery)
require(Biobase)
require(limma)
require(caret)
require(reshape)
require(sampling)
require(genefilter)
library(MetaIntegrator)
library(data.table)
library(survJamda)




# Dataset1 <- getGEOData("GSE102238")
# Dataset1 <- Dataset1$originalData$GSE102238

# ## Colon (GEO)
# Dataset2 <- getGEOData("GSE103479")
# Dataset2 <- Dataset2$originalData$GSE103479
#   
# # Skin SCC (GEO)
# Dataset3 <- getGEOData("GSE86544")
# Dataset3 <- Dataset3$originalData$GSE86544
# 
# # Prostate (GEO)
# Dataset4 <- getGEOData("GSE7055")
# Dataset4 <- Dataset4$originalData$GSE7055_GPL571
# 
# save(Dataset2, Dataset3, Dataset4, file = "./Data/Datasets.rda")

load("./Data/Datasets.rda")

# ## Head and Neck TCGA pheno
# X <- read.delim("/Users/mohamedomar/Desktop/Research/hnsc_tcga_clinical_data.tsv")
# 
# # CholangioCarcinoma
# X <- read.delim("/Users/mohamedomar/Desktop/Research/chol_tcga_clinical_data.tsv")
# 
# ## Colon
# X <- read.delim("/Users/mohamedomar/Desktop/Research/coadread_tcga_clinical_data.tsv")

## Get expression
#expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr
# expr5 (HNSCC)
expr5 <- read.delim("/Volumes/Macintosh/Research/Projects/PeriNeural/Data/hnsc_tcga/data_RNA_Seq_v2_expression_median.txt")
# expr6 (CholangioCarcinoma)
expr6 <- read.delim("/Volumes/Macintosh/Research/Projects/PeriNeural/Data/chol_tcga/data_RNA_Seq_v2_expression_median.txt")
# expr7 (Colon)
expr7 <- read.delim("/Volumes/Macintosh/Research/Projects/PeriNeural/Data/coadread_tcga/data_RNA_Seq_v2_expression_median.txt")


## Check if expressions are normalized and log-scaled
# boxplot(expr1)
# boxplot(expr2)
# boxplot(expr3) ## Z-transformed
# boxplot(expr4)
# 
## Get phenotype
#Pheno1 <- Dataset1$pheno
Pheno2 <- Dataset2$pheno
Pheno3 <- Dataset3$pheno
Pheno4 <- Dataset4$pheno
Pheno5 <- read.delim("/Volumes/Macintosh/Research/Projects/PeriNeural/Data/hnsc_tcga/hnsc_tcga_clinical_data.tsv")
Pheno6 <- read.delim("/Volumes/Macintosh/Research/Projects/PeriNeural/Data/chol_tcga/chol_tcga_clinical_data.tsv")
Pheno7 <- read.delim("/Volumes/Macintosh/Research/Projects/PeriNeural/Data/coadread_tcga/coadread_tcga_clinical_data.tsv")

## Annotation

# expr2
rownames(expr2) <- Dataset2$keys
summary(is.na(rownames(expr2)))
expr2 <- expr2[!is.na(rownames(expr2)), ]
expr2 <- expr2[!(rownames(expr2) == ""), ]

rownames(expr2) <- gsub("-","", rownames(expr2))
rownames(expr2) <- gsub("_","",rownames(expr2))
sel <- which(apply(expr2, 1, function(x) all(is.finite(x)) ))
expr2 <- expr2[sel, ]
dim(expr2)


X2 <- expr2
ffun <- filterfun(pOverA(p = 0.8, A = 100))
filt2 <- genefilter(2^X2, ffun)
expr2 <- expr2[filt2, ]
dim(expr2)

expr2 <- t(scale(t(expr2), center = T, scale = T))
#####################

# expr3
rownames(expr3) <- Dataset3$keys
summary(is.na(rownames(expr3)))
expr3 <- expr3[!is.na(rownames(expr3)), ]
expr3 <- expr3[!(rownames(expr3) == ""), ]

rownames(expr3) <- gsub("-","", rownames(expr3))
rownames(expr3) <- gsub("_","",rownames(expr3))
sel <- which(apply(expr3, 1, function(x) all(is.finite(x)) ))
expr3 <- expr3[sel, ]

dim(expr3)
# ALready scaled and centered
#########################

# expr4
rownames(expr4) <- Dataset4$keys
summary(is.na(rownames(expr4)))
expr4 <- expr4[!is.na(rownames(expr4)), ]
expr4 <- expr4[!(rownames(expr4) == ""), ]

rownames(expr4) <- gsub("-","", rownames(expr4))
rownames(expr4) <- gsub("_","",rownames(expr4))
sel <- which(apply(expr4, 1, function(x) all(is.finite(x)) ))
expr4 <- expr4[sel, ]

dim(expr4)

X4 <- expr4
filt4 <- genefilter(2^X4, ffun)
expr4 <- expr4[filt4, ]

expr4 <- t(scale(t(expr4), center = T, scale = T))

# expr5
expr5 <- expr5[!duplicated(expr5$Hugo_Symbol), ]
expr5 <- expr5[!is.na(expr5$Hugo_Symbol), ]
rownames(expr5) <- expr5$Hugo_Symbol
expr5$Hugo_Symbol <- NULL
expr5$Entrez_Gene_Id <- NULL
expr5 <- expr5[!is.na(rownames(expr5)), ]
expr5 <- expr5[!(rownames(expr5) == ""), ]

expr5 <- expr5[unique(gsub("-","", rownames(expr5))), ]
rownames(expr5) <- gsub("_","",rownames(expr5))
sel <- which(apply(expr5, 1, function(x) all(is.finite(x)) ))
expr5 <- expr5[sel, ]

# C <- as.matrix(expr5)
# plot(density(log2(C)), xlim=c(-5,20))

# Filter expr9 to include only genes with > 1 RSEM in 50% of samples
X <- expr5
ffun <- filterfun(pOverA(p = 0.8, A = 1))

filt <- genefilter(X,ffun)
expr5 <- expr5[filt,]
dim(expr5)

expr5 <- log2(expr5 + 1)
#boxplot(expr5)
expr5 <- t(scale(t(expr5), center = T, scale = T))

# expr6
expr6 <- expr6[!duplicated(expr6$Hugo_Symbol), ]
expr6 <- expr6[!is.na(expr6$Hugo_Symbol), ]
rownames(expr6) <- expr6$Hugo_Symbol
expr6$Hugo_Symbol <- NULL
expr6$Entrez_Gene_Id <- NULL
expr6 <- expr6[!is.na(rownames(expr6)), ]
expr6 <- expr6[!(rownames(expr6) == ""), ]

expr6 <- expr6[unique(gsub("-","", rownames(expr6))), ]
rownames(expr6) <- gsub("_","",rownames(expr6))
sel <- which(apply(expr6, 1, function(x) all(is.finite(x)) ))
expr6 <- expr6[sel, ]

# C <- as.matrix(expr6)
# plot(density(log2(C)), xlim=c(-5,20))

# Filter expr9 to include only genes with > 1 FPKM in all samples
X <- expr6
filt <- genefilter(X,ffun)
expr6 <- expr6[filt,]
dim(expr6)

expr6 <- log2(expr6 + 1)
#boxplot(expr6)
expr6 <- t(scale(t(expr6), center = T, scale = T))

# expr7
expr7 <- expr7[!duplicated(expr7$Hugo_Symbol), ]
expr7 <- expr7[!is.na(expr7$Hugo_Symbol), ]
rownames(expr7) <- expr7$Hugo_Symbol
expr7$Hugo_Symbol <- NULL
expr7$Entrez_Gene_Id <- NULL
expr7 <- expr7[!is.na(rownames(expr7)), ]
expr7 <- expr7[!(rownames(expr7) == ""), ]

expr7 <- expr7[unique(gsub("-","", rownames(expr7))), ]
rownames(expr7) <- gsub("_","",rownames(expr7))
sel <- which(apply(expr7, 1, function(x) all(is.finite(x)) ))
expr7 <- expr7[sel, ]

#C <- as.matrix(expr7)
#plot(density(log2(C)), xlim=c(-5,20))

# Filter expr9 to include only genes with > 1 FPKM in all samples
X <- expr7
filt <- genefilter(X,ffun)
expr7 <- expr7[filt,]
dim(expr7)

expr7 <- log2(expr7 + 1)
#boxplot(expr7)
expr7 <- t(scale(t(expr7), center = T, scale = T))
#############################

## Modify Pheno 2
Pheno2$Perineural <- Pheno2$`perineural_invasion:ch1`
Pheno2 <- Pheno2[!(Pheno2$Perineural == "NA"), ]
Pheno2$Perineural <- factor(Pheno2$Perineural, levels = c("Yes", "No"))


expr2 <- expr2[, colnames(expr2) %in% rownames(Pheno2)]
dim(expr2)
all(rownames(Pheno2) == colnames(expr2))
#########

## Modify Pheno3
Pheno3$Perineural <- Pheno3$`pni:ch1`
Pheno3$Perineural[Pheno3$Perineural == "no neural involvement"] <- "No"
Pheno3$Perineural[Pheno3$Perineural == "clinical perineural invasion"] <- "Yes"
Pheno3$Perineural[Pheno3$Perineural == "incidental perineural invasion"] <- "Yes"

Pheno3$Perineural <- factor(Pheno3$Perineural, levels = c("Yes", "No"))
############

## Modify Pheno4
Pheno4$Perineural <- Pheno4$`perineural invasion:ch1`
Pheno4$Perineural <- factor(Pheno4$Perineural, levels = c("Yes", "No"))
##################

## Modify Pheno5
Pheno5$Sample.ID <- gsub("-", ".", Pheno5$Sample.ID)
Pheno5 <- Pheno5[match(colnames(expr5), Pheno5$Sample.ID), ]
rownames(Pheno5) <- Pheno5$Sample.ID
all(rownames(Pheno5) == colnames(expr5))

Pheno5$Perineural <- Pheno5$Perineural.Invasion
Pheno5 <- Pheno5[!is.na(Pheno5$Perineural), ]
Pheno5$Perineural <- as.factor(Pheno5$Perineural)
levels(Pheno5$Perineural) <- c("No", "Yes")
table(Pheno5$Perineural)
Pheno5$Perineural <- factor(Pheno5$Perineural, levels = c("Yes", "No"))

expr5 <- expr5[, colnames(expr5) %in% rownames(Pheno5)]
all(rownames(Pheno5) == colnames(expr5))

## Modify Pheno6
Pheno6$Sample.ID <- gsub("-", ".", Pheno6$Sample.ID)
Pheno6 <- Pheno6[Pheno6$Sample.ID %in% colnames(expr6), ]
rownames(Pheno6) <- Pheno6$Sample.ID
all(rownames(Pheno6) == colnames(expr6))

Pheno6$Perineural <- Pheno6$Perineural.Invasion
Pheno6 <- Pheno6[!is.na(Pheno6$Perineural), ]
Pheno6$Perineural <- as.factor(Pheno6$Perineural)
levels(Pheno6$Perineural) <- c("No", "Yes")
table(Pheno6$Perineural)
Pheno6$Perineural <- factor(Pheno6$Perineural, levels = c("Yes", "No"))

expr6 <- expr6[, colnames(expr6) %in% rownames(Pheno6)]
all(rownames(Pheno6) == colnames(expr6))


## Modify Pheno7 
Pheno7$Sample.ID <- gsub("-", ".", Pheno7$Sample.ID)
Pheno7 <- Pheno7[match(colnames(expr7), Pheno7$Sample.ID), ]
rownames(Pheno7) <- Pheno7$Sample.ID
all(rownames(Pheno7) == colnames(expr7))

Pheno7$Perineural <- Pheno7$Perineural.Invasion
Pheno7 <- Pheno7[!is.na(Pheno7$Perineural), ]
Pheno7$Perineural <- as.factor(Pheno7$Perineural)
levels(Pheno7$Perineural) <- c("No", "Yes")
table(Pheno7$Perineural)
Pheno7$Perineural <- factor(Pheno7$Perineural, levels = c("Yes", "No"))

expr7 <- expr7[, colnames(expr7) %in% rownames(Pheno7)]
all(rownames(Pheno7) == colnames(expr7))

#######################################################

#########################################################################

## Combine expression and phenotype

AllPheno <- list(Pheno2, Pheno4, Pheno5, Pheno6, Pheno7)
names(AllPheno) <- c("GSE103479", "GSE7055", "hnsc_tcga", "chol_tcga", "coadread_tcga")

AllExpr <- list(expr2, expr4, expr5, expr6, expr7)
names(AllExpr) <- c("GSE103479", "GSE7055", "hnsc_tcga", "chol_tcga", "coadread_tcga")

Group_Perineural <- list(Pheno2$Perineural, Pheno4$Perineural, Pheno5$Perineural, Pheno6$Perineural, Pheno7$Perineural)
#names(Group_Perineural) <- c("GSE103479", "GSE7055", "hnsc_tcga", "chol_tcga", "coadread_tcga")
Group_Perineural <- unlist(Group_Perineural)
table(Group_Perineural)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(AllExpr, rownames))

### Filter expression for the common genes
exprsMetastasis <- mapply(x=AllExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

##########



# Train <- c("GSE103479", "GSE7055", "chol_tcga","hnsc_tcga")
# Test <- c("coadread_tcga")
# # 
# # ### Pheno
# Pheno_Train <- AllPheno[Train]
# Pheno_Test <- AllPheno[Test]
# 
# ## Training Mat and Group
# mixTrainMat <- exprsMetastasis[Train]
# mixTrainMat <- do.call("cbind", mixTrainMat)
# 
# mixTrainGroup <- Group_Perineural[Train]
# mixTrainGroup <- unlist(mixTrainGroup)
# table(mixTrainGroup)
# 
# names(mixTrainGroup) <- colnames(mixTrainMat)
# 
# ## Testing Mat and Group
# mixTestMat <- exprsMetastasis[Test]
# mixTestMat <- mixTestMat$coadread_tcga
# 
# mixTestGroup <- Group_Perineural[Test]
# mixTestGroup <- unlist(mixTestGroup)
# table(mixTestGroup)
# 
# names(mixTestGroup) <- colnames(mixTestMat)


## Assemble in one data frame
AllMat <- do.call("cbind", exprsMetastasis)

## Check if sample names are identical
all(colnames(AllMat) == names(Group_Perineural))
#
# ##############################################
# ## T_Stage

AllPheno$GSE103479$t_Stage <- as.factor(AllPheno$GSE103479$`tstage:ch1`)

AllPheno$hnsc_tcga$t_Stage <- as.factor(AllPheno$hnsc_tcga$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
levels(AllPheno$hnsc_tcga$t_Stage) <- c("pT0", "pT1", "pT2", "pT3", "pT4", "pT4a", "pT4b", "TX")

AllPheno$chol_tcga$t_Stage <- as.factor(AllPheno$chol_tcga$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
levels(AllPheno$chol_tcga$t_Stage) <- c("pT1", "pT2", "pT2a", "pT2b", "pT3", "pT4")

AllPheno$coadread_tcga$t_Stage <- as.factor(AllPheno$coadread_tcga$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
levels(AllPheno$coadread_tcga$t_Stage) <- c("pT1", "pT2", "pT3", "pT4", "pT4a", "pT4b", "pTis")


# ### Covariates of relevance select complete cases: T_STAGE
All_t_Stage <- lapply(AllPheno, function(x) {
  i <- grep("t_Stage", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

All_t_Stage <- factor(unlist(All_t_Stage))
levels(All_t_Stage)[levels(All_t_Stage) == "TX"] <- NA
All_t_Stage[All_t_Stage == ""] <- NA

#####################################################
# N_Stage
AllPheno$GSE103479$n_Stage <- as.factor(AllPheno$GSE103479$`nstage:ch1`)

AllPheno$hnsc_tcga$n_Stage <- as.factor(AllPheno$hnsc_tcga$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Regional.Lymph.Node.N.Stage)
levels(AllPheno$hnsc_tcga$n_Stage) <- c("pN0", "pN1", "pN2", "pN2a", "pN2b", "pN2c", "pN3", "NX")

AllPheno$chol_tcga$n_Stage <- as.factor(AllPheno$chol_tcga$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code)
levels(AllPheno$chol_tcga$n_Stage) <- c("pN0", "pN1", "NX")

AllPheno$coadread_tcga$n_Stage <- as.factor(AllPheno$coadread_tcga$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code)
levels(AllPheno$coadread_tcga$n_Stage) <- c("pN0", "pN1", "pN1a", "pN1b", "pN1c", "pN2", "pN2a", "pN2b", "NX")



# ### Covariates of relevance select complete cases: N_STAGE
All_n_Stage <- lapply(AllPheno, function(x) {
  i <- grep("n_Stage", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

All_n_Stage <- factor(unlist(All_n_Stage))
levels(All_n_Stage)[levels(All_n_Stage) == "NX"] <- NA
All_n_Stage[All_n_Stage == ""] <- NA

#####################################################
# M_Stage
AllPheno$GSE103479$m_Stage <- as.factor(AllPheno$GSE103479$`mstage:ch1`)

AllPheno$hnsc_tcga$m_Stage <- as.factor(AllPheno$hnsc_tcga$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code)
levels(AllPheno$hnsc_tcga$m_Stage) <- c("pM0", "pM1", "pMX")

AllPheno$chol_tcga$m_Stage <- as.factor(AllPheno$chol_tcga$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code)
levels(AllPheno$chol_tcga$m_Stage) <- c("pM0", "pM1", "pMX")

AllPheno$coadread_tcga$m_Stage <- as.factor(AllPheno$coadread_tcga$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code)
levels(AllPheno$coadread_tcga$m_Stage) <- c("pM0", "pM1", "pM1", "pM1", "pMX")



# ### Covariates of relevance select complete cases: M_STAGE
All_m_Stage <- lapply(AllPheno, function(x) {
  i <- grep("m_Stage", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

All_m_Stage <- factor(unlist(All_m_Stage))
levels(All_m_Stage)[levels(All_m_Stage) == "pMX"] <- NA
All_m_Stage[All_m_Stage == ""] <- NA

####################################################################
# # Age
AllPheno$GSE103479$Age <- as.character(AllPheno$GSE103479$`age diagnosis:ch1`)
AllPheno$hnsc_tcga$Age <- as.character(AllPheno$hnsc_tcga$Diagnosis.Age)
AllPheno$chol_tcga$Age <- as.character(AllPheno$chol_tcga$Diagnosis.Age)
AllPheno$coadread_tcga$Age <- as.character(AllPheno$coadread_tcga$Diagnosis.Age)
#
# ### Covariates of relevance select complete cases: AGE
allAGE <- lapply(AllPheno, function(x) {
  i <- grep("^Age", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})
#
allAGE <- unlist(allAGE)
allAGE[allAGE == ""] <- NA
# ######################################################
#
# ## Gender
AllPheno$GSE103479$Sex <- as.factor(AllPheno$GSE103479$`gender:ch1`)

AllPheno$GSE7055$Sex <- rep("Male", 57)
AllPheno$GSE7055$Sex <- as.factor(AllPheno$GSE7055$Sex)

AllPheno$hnsc_tcga$Sex <- as.factor(AllPheno$hnsc_tcga$Sex)
AllPheno$chol_tcga$Sex <- as.factor(AllPheno$chol_tcga$Sex)
AllPheno$coadread_tcga$Sex <- as.factor(AllPheno$coadread_tcga$Sex)

# ### Covariates of relevance select complete cases: Sex
AllSex <- lapply(AllPheno, function(x) {
  i <- grep("Sex", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

AllSex <- factor(unlist(AllSex))
AllSex[AllSex == ""] <- NA

# ############################################################
# ### Assemble in one data.frame and turn numeric
covs <- data.frame(T_Stage = All_t_Stage, N_Stage = All_n_Stage, M_Stage =All_m_Stage, Age = allAGE, Sex = AllSex)
#
# ### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))))
#
# ###########################################################################
# ###SAMPLING
#
# ### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(Group_Perineural == "Yes"),
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs) * 0.3),
  comment=TRUE, method=1)

### Show
apply(covs[, -ncol(covs),drop=FALSE], 2, table, Group_Perineural, trainingOrTesting)

# ### Subset Training
mixTrainMat <- AllMat[ , trainingOrTesting == 0]
mixTrainGroup <- Group_Perineural[ trainingOrTesting == 0]
table(mixTrainGroup)
### Subset Testing
mixTestMat <- AllMat[ , trainingOrTesting == 1]
mixTestGroup <- Group_Perineural[ trainingOrTesting == 1]
table(mixTestGroup)
# #


##########################################################################
## Divide data into training and testing data
# set.seed(333)
# ind <- createDataPartition(y = Group_Perineural, p = 0.5, list = FALSE)
# 
# mixTrainGroup <- Group_Perineural[ind]
# mixTestGroup <- Group_Perineural[-ind]
# 
# mixTrainMat <- AllMat[, ind]
# mixTestMat <- AllMat[, -ind]

## Save
save(exprsMetastasis, Group_Perineural, AllPheno, mixTestMat, mixTrainMat, mixTrainGroup, mixTestGroup, file = "./Objs/PeriNeuralData.rda")
