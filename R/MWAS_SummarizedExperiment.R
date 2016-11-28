MWAS_SummarizedExperiment = function(metabo_matrix, clinical_matrix,
    sample_type) {

    ## Check that input data are correct
    if (!is.matrix(metabo_matrix) | !is.numeric(metabo_matrix) ) {
        stop("metabo_matrix must be a numeric matrix")
    }
    if (!is.matrix(clinical_matrix) | !is.numeric(clinical_matrix)) {
        stop("clinical_matrix must be a numeric matrix")
    }
    if (nrow(metabo_matrix) != nrow(clinical_matrix)) {
        stop("metabo_matrix nrow  must be consistent with clinical_matrix nrow")
    }
    if (is.null(rownames(metabo_matrix)) | is.null(rownames(clinical_matrix))) {
        stop("metabo_matrix or clinical_matrix rownames are missing ")
    }
    if (identical(rownames(metabo_matrix), rownames(clinical_matrix)) == FALSE) {
        stop("metabo_matrix and clinical_matrix rownames must be identical")
    }
    if (is.null(colnames(metabo_matrix)) | is.null(colnames(clinical_matrix))) {
        stop("metabo_matrix or clinical_matrix colnames are missing ")
    }
    if (!is.vector(sample_type)) {
        stop("sample_type must be a numeric vector")
    }
    if (length(setdiff(unique(sample_type), c(0, 1)) > 0)) {
        stop("sample_type must be a numeric vector with values 0 or 1")
    }
    if (length(sample_type) != nrow(clinical_matrix)) {
        stop("sample_type length must be consistent with clinical_matrix nrow")
    }
    clinical_matrixN = cbind(clinical_matrix, sample_type)
    colnames(clinical_matrixN) = c(colnames(clinical_matrix), "sample_type")

    experiment = SummarizedExperiment(assays = list(metabolic_data = t(metabo_matrix)),
                                      colData = clinical_matrixN)
    return(experiment)
}

