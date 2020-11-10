#################################################################################
## Mohamed Omar
## 07/08/2019
## Goal: Survival analysis and Cox proportional hazards model based on probablity of perineural invasion
##############################################################################
setwd("/Users/mohamedomar/Documents/Research/Projects/PeriNeural")

rm(list = ls())


library(survival)
library(survminer)
library(maxstat)


load("./Objs/val_pheno_prob.rda")

## Divide the dataset into quantiles based on the Meta score risk.
# quantiles <- quantile(val_Pheno$Prob, probs=c(0.33333,0.66667))
# val_Pheno$Probability <- val_Pheno$Prob
# val_Pheno[which(val_Pheno$Prob <= quantiles[1]), "Probability"] = "low"
# val_Pheno[which(val_Pheno$Prob > quantiles[1] &  val_Pheno$Prob <= quantiles[2]),"Probability"] = "int"
# val_Pheno[which(val_Pheno$Prob > quantiles[2]),"Probability"] = "high"
# table(val_Pheno$Probability)
# 


val_Pheno$Time <- as.numeric(val_Pheno$Overall.Survival..Months.)

val_Pheno$Event <- as.numeric(val_Pheno$Overall.Survival.Status)
val_Pheno$Event[val_Pheno$Event == 2] <- 0 # 0 = Alive, 1 = Dead
table(val_Pheno$Event)

## Keep only relevant information from the phenotype table
surv_data <- val_Pheno[,c("Time","Event","Prob","Surgical.Margin.Resection.Status", "Lymphovascular.invasion.indicator", "American.Joint.Committee.on.Cancer.Tumor.Stage.Code", "Diagnosis.Age", "Sex", "Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code", "American.Joint.Committee.on.Cancer.Metastasis.Stage.Code")]

CutPoint_Prob <- surv_cutpoint(data = surv_data, time = "Time", event = "Event", variables = "Prob")
CutPoint_Prob

SurvData_categ <- surv_categorize(CutPoint_Prob)
surv_data <- cbind(SurvData_categ, surv_data[,4:length(colnames(surv_data))])

surv_data$Prob <- factor(surv_data$Prob, levels = c("low","high"))
table(surv_data$Prob)

# Create a survival object
surv_data.surv <-  with(surv_data, Surv(Time, Event == 1))

#Calculate p-value
survdifftest <- survdiff(surv_data.surv ~ Prob, data = surv_data)
survpvalue <- 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)

png(filename = "./Figs/KM_Survival.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Prob, data = surv_data),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE,
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Overall survival by the PNI probability in the testing data", xlab = "Time in months", ylab = "Overall Survival")
dev.off()

#############################################################################
## Cox proportional hazards model

Cox_data <- surv_data
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data$Age <- ifelse(Cox_data$Diagnosis.Age >= 75, ">=75 years", "<75 years")
Cox_data$Age <- factor(Cox_data$Age, levels = c("<75 years", ">=75 years"))
Cox_data$Diagnosis.Age <- NULL

Cox_data$Sex <- factor(Cox_data$Sex, levels = c("Male", "Female"))
table(Cox_data$Sex)

Cox_data$Surgical_Margin <- as.factor(Cox_data$Surgical.Margin.Resection.Status)
levels(Cox_data$Surgical_Margin) <- c("Negative", "Positive", "Positive", "UnKnown")
Cox_data$Surgical_Margin <- factor(Cox_data$Surgical_Margin, levels = c("UnKnown", "Negative", "Positive"))
Cox_data$Surgical.Margin.Resection.Status <- NULL
levels(Cox_data$Surgical_Margin)[levels(Cox_data$Surgical_Margin) == "UnKnown"] <- NA
table(Cox_data$Surgical_Margin)

Cox_data$Lymphovascular_invasion <- as.factor(Cox_data$Lymphovascular.invasion.indicator)
Cox_data$Lymphovascular.invasion.indicator <- NULL

Cox_data$T_Stage <- as.factor(Cox_data$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
levels(Cox_data$T_Stage) <- c("Tis-T1-T2", "Tis-T1-T2", "T3-T4", "T3-T4", "T3-T4", "T3-T4","Tis-T1-T2")
table(Cox_data$T_Stage)

Cox_data$N_Stage <- as.factor(Cox_data$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code)
levels(Cox_data$N_Stage) <- c("N0", "N1-N2", "N1-N2", "N1-N2", "N1-N2", "N1-N2", "N1-N2", "N1-N2", "NX")
levels(Cox_data$N_Stage)[levels(Cox_data$N_Stage) == "NX"] <- NA
table(Cox_data$N_Stage)


Cox_data$M_Stage <- as.factor(Cox_data$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code)
levels(Cox_data$M_Stage) <- c("M0", "M1", "M1", "M1", "MX")
levels(Cox_data$M_Stage)[levels(Cox_data$M_Stage) == "MX"] <- NA
table(Cox_data$M_Stage)

## Univariate Cox model
covariates <- c("Prob", "Lymphovascular_invasion",  "Age", "T_Stage", "N_Stage", "M_Stage", "Surgical_Margin")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = Cox_data)})
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

## Multivariate Cox model
CoxModel <- coxph(Surv(Time, Event)~ Prob+Lymphovascular_invasion+T_Stage+N_Stage+M_Stage+Surgical_Margin+Age, data=Cox_data)

png(filename = "./Figs/CoxModel.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel, data = Cox_data, main = "Cox proportional hazards model")
dev.off()

#########
cox.zph.fit <- cox.zph(CoxModel)

ggcoxzph(cox.zph.fit, font.main = 8)
