######################### TEST MWASTools ###################################
library(RUnit)

## Load data
data(metabo_SE)
data(targetMetabo_SE)

#### Test MWAS_stats ####

## Test for association between diabetes and target_metabolites
T2D_model1 = MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         assoc_method = "logistic", mt_method = "BH",
                         output = "pvalues")

T2D_model2 = MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         assoc_method = "logistic", mt_method = "BY",
                         output = "pvalues")

T2D_model3 = MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         assoc_method = "logistic", mt_method = "bonferroni",
                         output = "pvalues")

T2D_model4 = MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         assoc_method = "logistic", mt_method = "none",
                         output = "pvalues")

## Test for association between diabetes and target_metabolites (age & gender adjusted)
T2D_model5 = MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         confounder_ids = c("Age", "Gender"),
                         assoc_method = "logistic", mt_method = "BH",
                         output = "pvalues")

## Test for association between BMI and target_metabolites
BMI_model1 = MWAS_stats (targetMetabo_SE, disease_id = "BMI",
                        assoc_method = "spearman", mt_method = "BH",
                        output = "pvalues")

## Test for association between BMI and target_metabolites (age & gender adjusted)
BMI_model2 = MWAS_stats (targetMetabo_SE, disease_id = "BMI",
                         confounder_ids = c("Age", "Gender"),
                         assoc_method = "spearman", mt_method = "BH",
                         output = "pvalues")

## Test for association between BMI and target_metabolites (output = models)
BMI_model3 = MWAS_stats (targetMetabo_SE, disease_id = "BMI",
                         confounder_ids = c("Age", "Gender"),
                         assoc_method = "spearman", mt_method = "BH",
                         output = "models")

## Attempt to do logistic regression with a non-binary predictive variable
BMI_model4 = try(MWAS_stats (targetMetabo_SE, disease_id = "BMI",
                             assoc_method = "logistic"), silent = TRUE)

## Attempto to apply MWAS_stats with invalid assoc_method
BMI_model5 = try(MWAS_stats (targetMetabo_SE, disease_id = "BMI", assoc_method = "X"),
                 silent = TRUE)

## Attempto to apply MWAS_stats with invalid assoc_method
BMI_model6 = try(MWAS_stats (targetMetabo_SE, disease_id = "BMI", assoc_method = "spearman",
                             mt_method = "X"), silent = TRUE)

## Tests
test.MWAS_stats <- function() {

    checkTrue(is.matrix(T2D_model1))
    checkTrue(is.matrix(T2D_model2))
    checkTrue(is.matrix(T2D_model3))
    checkTrue(is.matrix(T2D_model4))
    checkTrue(is.matrix(T2D_model5))
    checkTrue(is.matrix(BMI_model1))
    checkTrue(is.matrix(BMI_model2))
    checkTrue(is.list(BMI_model3))
    checkException(BMI_model4)
    checkException(BMI_model5)
    checkException(BMI_model6)
}
test.MWAS_stats()

################################################################################

#### Test QC_CV ####
metabo_CV =  QC_CV (metabo_SE)

## Tests
test.metabo_CV <- function() {
  checkTrue(is.vector(metabo_CV))
}
test.metabo_CV()

################################################################################

#### Test QC_PCA ####
PCA_model =  QC_PCA (targetMetabo_SE)

## Tests
test.PCA_model <- function() {
  checkTrue(is.list(PCA_model))
}
test.PCA_model()

################################################################################

#### Test MWAS_filter ####

## Test for association between diabetes and target_metabolites
T2D_model1 <- MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         assoc_method = "logistic")

## Filter T2D_model1 by pvalue
pvalue_filter1 <- MWAS_filter(T2D_model1, type = "pvalue", alpha_th = 0.001)

## Attempt to filter T2D_model1 with too stringent criteria
pvalue_filter2 <- try(MWAS_filter(T2D_model1, type = "pvalue", alpha_th = 0), silent = TRUE)

## Tests
test.MWAS_filter <- function() {
  checkTrue(is.matrix(pvalue_filter1))
  checkException(pvalue_filter2)
}
test.MWAS_filter()

################################################################################

#### Test CV_filter ####

## Calculate CVs
CV_metabo <-  QC_CV (metabo_SE)

## Filter metabolic_data by CV
metabo_CVfiltered <- CV_filter(metabo_SE, CV_metabo, CV_th = 0.30)
metabo_CVfiltered2 <- CV_filter(metabo_SE, CV_metabo, CV_th = 0.15)

## Tests
test.CV_filter <- function() {
  checkTrue(nrow(metabo_SE) > nrow(metabo_CVfiltered))
  checkTrue(nrow(metabo_CVfiltered) > nrow(metabo_CVfiltered2))
}
test.CV_filter()

