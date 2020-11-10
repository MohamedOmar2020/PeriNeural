#################################################################################
## Mohamed Omar
## 07/08/2019
## Goal: Survival by genes
##############################################################################


rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")

############################################################################

## Load necessary packages
library(survival)
library(survminer)
library(maxstat)



## Load data
load("./Objs/PeriNeuralData.rda")

usedTestMat <- normalize.quantiles(mixTestMat)
rownames(usedTestMat) <- rownames(mixTestMat)
colnames(usedTestMat) <- colnames(mixTestMat)

load("./Objs/val_pheno_prob.rda")
load("./Objs/keep.rda")

val_Pheno$Time <- as.numeric(val_Pheno$Overall.Survival..Months.)

val_Pheno$Event <- as.numeric(val_Pheno$Overall.Survival.Status)
val_Pheno$Event[val_Pheno$Event == 2] <- 0 # 0 = Alive, 1 = Dead
table(val_Pheno$Event)

## Keep only relevant information from the phenotype table
Pheno_surv <- val_Pheno[,c("Time","Event")]

# Keep only the top 20 genes
Keep <- Keep[1:20]
expr_surv <- usedTestMat[Keep, ]

# create a merged pdata and Z-scores object
surv_data <- data.frame(Pheno_surv, t(expr_surv))

#################################################################################
## Define optimal cutpoints for each gene (converting the absolute expression into categorical low/high expression)

CutPoint_Genes <- surv_cutpoint(data = surv_data, time = "Time", event = "Event", variables = Keep)
CutPoint_Genes

SurvData_Genes <- surv_categorize(CutPoint_Genes)

########################################################################  
## Fit genes

Fit_ITGA5 <- survfit(Surv(Time, Event) ~ ITGA5 , data = SurvData_Genes)
Fit_FAM13A <- survfit(Surv(Time, Event) ~ FAM13A, data = SurvData_Genes)
Fit_SRRD <- survfit(Surv(Time, Event) ~ SRRD, data = SurvData_Genes)
Fit_TIMP3 <- survfit(Surv(Time, Event) ~ TIMP3, data = SurvData_Genes)
Fit_IL17RB <- survfit(Surv(Time, Event) ~ IL17RB, data = SurvData_Genes)
Fit_BCR <- survfit(Surv(Time, Event) ~ BCR, data = SurvData_Genes)
Fit_LGALS8 <- survfit(Surv(Time, Event) ~ LGALS8, data = SurvData_Genes)
Fit_CBFA2T2 <- survfit(Surv(Time, Event) ~ CBFA2T2, data = SurvData_Genes)
Fit_LAMB3 <- survfit(Surv(Time, Event) ~ LAMB3, data = SurvData_Genes)
Fit_POMT1 <- survfit(Surv(Time, Event) ~ POMT1, data = SurvData_Genes)
Fit_TM7SF3 <- survfit(Surv(Time, Event) ~ TM7SF3, data = SurvData_Genes)
Fit_GATM <- survfit(Surv(Time, Event) ~ GATM, data = SurvData_Genes)
Fit_ERAP1 <- survfit(Surv(Time, Event) ~ ERAP1, data = SurvData_Genes)
Fit_NDC1 <- survfit(Surv(Time, Event) ~ NDC1, data = SurvData_Genes)
Fit_SLC27A5 <- survfit(Surv(Time, Event) ~ SLC27A5, data = SurvData_Genes)
Fit_HES1 <- survfit(Surv(Time, Event) ~ HES1, data = SurvData_Genes)
Fit_EFS <- survfit(Surv(Time, Event) ~ EFS, data = SurvData_Genes)
Fit_CUL9 <- survfit(Surv(Time, Event) ~ CUL9, data = SurvData_Genes)
Fit_IGSF3 <- survfit(Surv(Time, Event) ~ IGSF3, data = SurvData_Genes)
Fit_DYNLT1 <- survfit(Surv(Time, Event) ~ DYNLT1, data = SurvData_Genes)


##################################################################

## Plot genes

Plot1 <- ggsurvplot(Fit_ITGA5,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "ITGA5")

Plot2 <- ggsurvplot(Fit_FAM13A,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "FAM13A")

Plot3  <-  ggsurvplot(Fit_SRRD,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "SRRD")

Plot4 <- ggsurvplot(Fit_TIMP3,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "TIMP3")

Plot5 <-  ggsurvplot(Fit_IL17RB,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "IL17RB")

Plot6 <-  ggsurvplot(Fit_BCR,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "BCR")

Plot7 <-  ggsurvplot(Fit_LGALS8,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "LGALS8")

Plot8 <-  ggsurvplot(Fit_CBFA2T2,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "CBFA2T2")

Plot9 <-  ggsurvplot(Fit_LAMB3,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "LAMB3")

Plot10 <-  ggsurvplot(Fit_POMT1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "POMT1")

Plot11 <-  ggsurvplot(Fit_TM7SF3,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "TM7SF3")

Plot12 <-  ggsurvplot(Fit_GATM,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "GATM")

Plot13 <-  ggsurvplot(Fit_ERAP1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "ERAP1")

Plot14 <-  ggsurvplot(Fit_NDC1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "NDC1")

Plot15 <-  ggsurvplot(Fit_SLC27A5,
                      risk.table = FALSE,
                      pval = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.y.text.col = FALSE,
                      risk.table.y.text = FALSE, title = "SLC27A5")

Plot16 <-  ggsurvplot(Fit_HES1,
                      risk.table = FALSE,
                      pval = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.y.text.col = FALSE,
                      risk.table.y.text = FALSE, title = "HES1")

Plot17 <-  ggsurvplot(Fit_EFS,
                      risk.table = FALSE,
                      pval = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.y.text.col = FALSE,
                      risk.table.y.text = FALSE, title = "EFS")

Plot18 <-  ggsurvplot(Fit_CUL9,
                      risk.table = FALSE,
                      pval = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.y.text.col = FALSE,
                      risk.table.y.text = FALSE, title = "CUL9")

Plot19 <-  ggsurvplot(Fit_IGSF3,
                      risk.table = FALSE,
                      pval = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.y.text.col = FALSE,
                      risk.table.y.text = FALSE, title = "IGSF3")

Plot20 <-  ggsurvplot(Fit_DYNLT1,
                      risk.table = FALSE,
                      pval = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.y.text.col = FALSE,
                      risk.table.y.text = FALSE, title = "DYNLT1")


PlotList <- list(Plot1, Plot2, Plot3, Plot4, Plot5, Plot6, Plot7, Plot8, Plot9, Plot10, Plot11, Plot12, Plot13, Plot14, Plot15, Plot16, Plot17, Plot18, Plot19, Plot20)
names(PlotList) <- Keep

Splot <- arrange_ggsurvplots(PlotList, title = "Survival plots of XGB top 20 genes",ncol = 5, nrow = 4)
ggsave("Surv_XGB_Genes.pdf", Splot, width = 20, height = 10)


################################################################
## Plot Cox Proportional Hazard Model
SurvData_Genes$ITGA5 <- factor(SurvData_Genes$ITGA5, levels = c("low", "high"))
SurvData_Genes$FAM13A <- factor(SurvData_Genes$FAM13A, levels = c("low", "high"))
SurvData_Genes$SRRD <- factor(SurvData_Genes$SRRD, levels = c("low", "high"))
SurvData_Genes$TIMP3 <- factor(SurvData_Genes$TIMP3, levels = c("low", "high"))
SurvData_Genes$IL17RB <- factor(SurvData_Genes$IL17RB, levels = c("low", "high"))
SurvData_Genes$BCR <- factor(SurvData_Genes$BCR, levels = c("low", "high"))
SurvData_Genes$LGALS8 <- factor(SurvData_Genes$LGALS8, levels = c("low", "high"))
SurvData_Genes$CBFA2T2 <- factor(SurvData_Genes$CBFA2T2, levels = c("low", "high"))
SurvData_Genes$LAMB3 <- factor(SurvData_Genes$LAMB3, levels = c("low", "high"))
SurvData_Genes$POMT1 <- factor(SurvData_Genes$POMT1, levels = c("low", "high"))
SurvData_Genes$TM7SF3 <- factor(SurvData_Genes$TM7SF3, levels = c("low", "high"))
SurvData_Genes$GATM <- factor(SurvData_Genes$GATM, levels = c("low", "high"))
SurvData_Genes$ERAP1 <- factor(SurvData_Genes$ERAP1, levels = c("low", "high"))
SurvData_Genes$NDC1 <- factor(SurvData_Genes$NDC1, levels = c("low", "high"))
SurvData_Genes$SLC27A5 <- factor(SurvData_Genes$SLC27A5, levels = c("low", "high"))
SurvData_Genes$HES1 <- factor(SurvData_Genes$HES1, levels = c("low", "high"))
SurvData_Genes$EFS <- factor(SurvData_Genes$EFS, levels = c("low", "high"))
SurvData_Genes$CUL9 <- factor(SurvData_Genes$CUL9, levels = c("low", "high"))
SurvData_Genes$IGSF3 <- factor(SurvData_Genes$IGSF3, levels = c("low", "high"))
SurvData_Genes$DYNLT1 <- factor(SurvData_Genes$DYNLT1, levels = c("low", "high"))


# fit multivariable Cox model 
CoxModel <- coxph(Surv(Time, Event) ~ ITGA5+FAM13A+SRRD+TIMP3+IL17RB+BCR+LGALS8+CBFA2T2+LAMB3+POMT1+TM7SF3+GATM+ERAP1+NDC1+SLC27A5+HES1+EFS+CUL9+IGSF3+DYNLT1, data = SurvData_Genes)

png("./Figs/CoxModel_Genes.png", width = 2000, height = 2000, res = 300)
ggforest(CoxModel, data = SurvData_Genes)
dev.off()

## Univariate Cox model
covariates <- c("Prob", "Lymphovascular_invasion",  "Age", "T_Stage", "N_Stage", "M_Stage", "Surgical_Margin")

univ_formulas <- sapply(Keep,
                        function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = SurvData_Genes)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res <- t(as.data.frame(univ_results, check.names = F))
as.data.frame(res)


