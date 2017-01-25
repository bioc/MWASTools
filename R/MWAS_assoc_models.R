## INTERNAL FUNCTIONS ##

### assoc_test ####
assoc_test = function(metabolite, disease, CF = NULL, method,
    output_AT = "pr_values") {

    cor_methods = c("pearson", "spearman", "kendall")

    if (method %in% cor_methods) {
        if (is.null(CF)) {
            model = suppressWarnings(cor.test(disease, metabolite,
                method = method))
            p_value = model$p.value
            r_value = model$estimate
        } else {
            model = pcor.test(disease, metabolite, CF, method = method)
            p_value = model["p.value"]
            r_value = model["estimate"]
        }
    } else if (method == "logistic") {
        if (is.null(CF)) {
            model = glm(disease ~ metabolite, family = binomial)
        } else {
            model = glm(disease ~ metabolite + CF, family = binomial)
        }
        coef = summary(model)$coefficients
        p_value = coef[2, 4]
        r_value = coef[2, 1]  # this is the beta coeff

    } else {
        # method = linear
        if (is.null(CF)) {
            model = glm(disease ~ metabolite, family = "gaussian")
        } else {
            model = glm(disease ~ metabolite + CF, family = "gaussian")
        }
        coef = summary(model)$coefficients
        p_value = coef[2, 4]
        r_value = coef[2, 1]  # this is the beta coeff
    }

    if (output_AT == "model") {
        return(model)
    } else if (output_AT == "pr_values") {
        pr_values = as.numeric(c(r_value, p_value))
        names(pr_values) = c("estimate", "p_value")
        return(pr_values)
    }
}

### MT_correction ####
MT_correction = function(pvalues, mt_method = mt_method) {
    if (mt_method == "qvalues") {
        bh_pvalues = p.adjust(pvalues, method = "BH")
        Q = try(qvalue(pvalues), silent = TRUE)
        if (grepl("Error", Q[1]) == TRUE) {
            stop("Could not calculate q-values. Chose another correction method")
        }
        hist(pvalues, prob = TRUE, col = "gray", xlab = "p-values",
            main = "p-values distribution")
        null_FP = try(head(Q, 2)$pi0, silent = TRUE)
        adjusted_pvalues = bh_pvalues * null_FP
        if (grepl("Error", adjusted_pvalues[1]) == TRUE) {
            stop("Could not calculate q-values. Chose another correction method")
        } else {
            adjusted_pvalues = bh_pvalues * null_FP
        }
    } else {
        adjusted_pvalues = p.adjust(pvalues, method = mt_method)
    }
    return(adjusted_pvalues)
}


## EXTERNAL FUNCTIONS ##

### MWAS_assoc ####
MWAS_stats = function(metabo_SE, disease_id, confounder_ids = NULL,
    assoc_method, mt_method = "BH", output = "pvalues", CV_metabo = NULL) {

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    metabo_matrix = t(assays(metabo_SE)$metabolic_data)
    metabo_ids = colnames(metabo_matrix)
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

    if (!is.null(CV_metabo)) {
        if (length(CV_metabo) != ncol(metabo_matrix)) {
            stop("CV_metabo length must be consistent with metabo_matrix dimension")
        }
    }

    all_var = cbind(metabo_matrix, disease, CF_matrix)

    ## Check that method is correct
    possible_methods = c("pearson", "spearman", "logistic", "linear",
        "kendall")
    if (length(intersect(assoc_method, possible_methods)) == 0) {
        stop("Invalid method. Valid methods are: pearson, spearman, kendall, logistic, linear")
    }

    ## Check that output is correct
    possible_outputs = c("models", "pvalues")
    if (length(intersect(output, possible_outputs)) == 0) {
        stop("Invalid method. Valid methods are: models,pvalues")
    }

    ## Check the number of samples
    if (nrow(all_var) < 10) {
        stop("The number of samples is too small")
    }

    ## Remove NA values
    if (nrow(na.omit(all_var)) < 10) {
        stop("The number of samples is too small after removing NA values")
    }

    ## Get NA index
    na_index = unique(as.numeric((which(is.na(all_var), arr.ind = TRUE)[,
        1])))

    if (length(na_index) > 0) {
        metabo_matrix = metabo_matrix[-c(na_index), ]
        disease = disease[-c(na_index)]
        CF_matrix = CF_matrix[-c(na_index), ]
    }

    ## Force metabo_matrix and CF_matrix to matrix (in case they
    ## are vectors)
    metabo_matrix = matrix(metabo_matrix, nrow = length(disease))

    if (!is.null(CF_matrix)) {
        CF_matrix = matrix(CF_matrix, nrow = length(disease))
    }

    ## Transform disease variable into a factor (logistic
    ## regression)
    if (assoc_method == "logistic") {
        disease = as.factor(disease)

        if (length(levels(disease)) != 2) {
            to_print = paste("Disease must be a binary numeric variable to perform",
                "logistic regression")
            stop(to_print)
        }
    }
    cols_metabo = split(t(metabo_matrix), row(t(metabo_matrix)))
    names(cols_metabo) = metabo_ids

    ## Run MWAS
    if (output == "models") {
        assoc_results = try(lapply(cols_metabo, assoc_test, disease = disease,
            CF = CF_matrix, method = assoc_method, output_AT = "model"),
            silent = TRUE)
        if (grepl("Error", assoc_results[1])) {
            stop("Association testing failed for at least one metabolic variable")
        }
        return(assoc_results)
    } else {
        assoc_results = try(lapply(cols_metabo, assoc_test, disease = disease,
            CF = CF_matrix, method = assoc_method, output_AT = "pr_values"),
            silent = TRUE)

        if (grepl("Error", assoc_results[1])) {
            stop("Association testing failed for at least one metabolic variable")
        }

        assoc_matrix = do.call(rbind, assoc_results)

        # Do mt_correction
        pvalues = as.numeric(assoc_matrix[, 2])
        mt_adjusted_pvalues = MT_correction(pvalues, mt_method = mt_method)
        assoc_matrix_adjusted = cbind(assoc_matrix, mt_adjusted_pvalues)

        if (is.null(CV_metabo)) {
            colnames(assoc_matrix_adjusted) = c("estimates",
                "pvalues", "adjusted_pvalues")
            return(assoc_matrix_adjusted)
        } else {
            assoc_matrix_adjusted = cbind(assoc_matrix_adjusted,
                CV_metabo)
            colnames(assoc_matrix_adjusted) = c("estimates",
                "pvalues", "adjusted_pvalues", "CV")
            return(assoc_matrix_adjusted)
        }
    }
}

### MWAS_filter ####

MWAS_filter = function(MWAS_matrix, type = "pvalue", alpha_th = 0.05,
                       CV_th = 0.3) {

    ## Check that input data are correct
    if (!is.matrix(MWAS_matrix) | !is.numeric(MWAS_matrix)) {
        stop("MWAS_matrix must be a numeric matrix. See function MWAS_stats()")
    }
    if (ncol(MWAS_matrix) < 3) {
        stop("MWAS_matrix must have at least 3 columns")
    }
    possible_types = c("pvalue", "CV", "all")
    if (type %in% possible_types == FALSE) {
        stop("Invalid type. Possible types are: pvalue, CV or all")
    }
    if (type != "pvalue") {
        if (ncol(MWAS_matrix) < 4) {
            stop("MWAS_matrix must have a CV column to apply CV filtering")
        }
    }

    ## Filter matrix
    MWAS_matrix = cbind(MWAS_matrix, original_index = 1:nrow(MWAS_matrix))
    index_CV = which(MWAS_matrix[, 4] >= CV_th)
    index_pval = which(MWAS_matrix[, 3] >= alpha_th)
    index_all = unique(c(index_CV, index_pval))

    if (type == "all") {
        if (length(index_all) == 0) {
            message ("All features satisfy the filtering parameters")
            return(MWAS_matrix)
        }
        if (length(index_all) == nrow(MWAS_matrix)) {
            stop("None of the metabolic features satisfies the filtering parameters")
        } else {
            MWAS_matrix = MWAS_matrix[-index_all, ]
        }
    } else if (type == "CV") {
        if (length(index_CV) == 0) {
            message ("All features satisfy the filtering parameters")
            return(MWAS_matrix)
        }
        if (length(index_CV) == nrow(MWAS_matrix)) {
            stop("None of the metabolic features satisfies the filtering parameters")
        } else {
            MWAS_matrix = MWAS_matrix[-index_CV, ]
        }
    } else {
        if (length(index_pval) == 0 ) {
            message ("All features satisfy the filtering parameters")
            return(MWAS_matrix)
        }
        if (length(index_pval) == nrow(MWAS_matrix)) {
            stop("None of the metabolic features satisfies the filtering parameters")
        } else {
            MWAS_matrix = MWAS_matrix[-index_pval, ]
        }
    }
    return(MWAS_matrix)
}

