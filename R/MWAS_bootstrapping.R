## INTERNAL FUNCTION ##

### model_bootstrap ##
model_bootstrap <- function(all_var, indices, assoc_method) {

    all_var = all_var[indices, ]
    rownames(all_var) = 1:nrow(all_var)

    if (ncol(all_var) > 2) {
        CF_matrix = all_var[, 3:ncol(all_var)]
    } else {
        CF_matrix = NULL
    }

    metabolite = all_var[, 1]
    disease = all_var[, 2]

    MWAS_results = assoc_test(metabolite = metabolite, disease = disease,
        CF = CF_matrix, output_AT = "pr_values", method = assoc_method)

    MWAS_results = MWAS_results[1]

    return(MWAS_results)
}

## EXTERNAL FUNCTION ##

### MWAS_bootstrapping ##
MWAS_bootstrapping = function(metabo_SE, metabolite_id, disease_id,
    confounder_ids = NULL, assoc_method, iterations = 10000) {

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    if ((metabolite_id[1] %in% rownames(metabo_SE)) == FALSE) {
        # only 1 metabolite_id is valid
        stop("metabolite_id is not included in metabolic_data")
    }
    ind_metabo = which(rownames(metabo_SE) == metabolite_id[1])
    metabolite = (assays(metabo_SE)$metabolic_data)[ind_metabo, ]

    clinical_variables = as.matrix(colData(metabo_SE))
    if ((disease_id[1] %in% colnames(clinical_variables)) ==
        FALSE) {
        stop("disease_id is not included in clinical_data")
    }
    disease = clinical_variables[, disease_id[1]]  # only 1 disease is valid.

    if (!is.null(confounder_ids)) {
        if (length(setdiff(confounder_ids, colnames(clinical_variables))) >
            0) {
            stop("One or more confounders are not included in clinical_data")
        }
        CF_matrix = clinical_variables[, confounder_ids]
        CF_matrix = as.matrix(CF_matrix, nrow = nrow(clinical_variables))
    } else {
        CF_matrix = NULL
    }

    all_var = cbind(metabolite, disease, CF_matrix)

    # Check that method is correct
    possible_methods = c("pearson", "spearman", "logistic", "linear",
        "kendall")
    if (length(intersect(assoc_method, possible_methods)) == 0) {
        stop("Invalid method. Valid methods are: pearson, spearman,kendall, logistic, linear")
    }

    # Check the number of samples
    if (nrow(all_var) < 10) {
        stop("The number of samples is too small")
    }

    # Remove NA values
    if (nrow(na.omit(all_var)) < 10) {
        stop("The number of samples is too small after removing NA values")
    }

    ## Get NA index
    na_index = unique(as.numeric((which(is.na(all_var), arr.ind = TRUE)[, 1])))

    if (length(na_index) > 0) {
        metabolite = metabolite[-na_index]
        disease = disease[-na_index]
        CF_matrix = CF_matrix[-na_index, ]
        all_var = cbind(metabolite, disease, CF_matrix)
    }

    # Transform disease variable into a factor (logistic
    # regression)
    if (assoc_method == "logistic") {
        disease = as.factor(disease)

        if (length(levels(disease)) != 2) {
            to_print = paste("It is not possible to perform logistic regression",
                "with the current disease argument")
            stop(to_print)
        }
    }

    set.seed(12345)  # ensure results reproducibility

    MWAS_boot = boot(data = all_var, model_bootstrap, R = iterations,
        sim = "ordinary", assoc_method = assoc_method)

    MWAS_bootM = summary(MWAS_boot)
    rownames(MWAS_bootM) = c("estimate")
    colnames(MWAS_bootM) = c("iterations", "original", "bias",
        "std.error", "Median")

    ## Calculate CI
    CI_est = try(boot.ci(MWAS_boot, type = "bca", index = 1)$bca,
        silent = TRUE)
    if (grepl("Error", CI_est[1])) {
        stop("Couldn't calculate 95% CI: try with a higher number of iterations")
    }
    CI_est_interval = paste(CI_est[4], CI_est[5], sep = ", ")

    MWAS_list = vector(mode = "list", length = 3)
    names(MWAS_list) = c("Boot", "Boot summary", "95% CI estimate")
    MWAS_list[[1]] = MWAS_boot
    MWAS_list[[2]] = MWAS_bootM
    MWAS_list[[3]] = CI_est_interval

    return(MWAS_list)

}

