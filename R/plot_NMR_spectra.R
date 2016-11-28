## EXTERNAL FUNCTION ##

plot_spectraNMR = function(metabo_SE, type = "l", lty = 1,
                           xlab = "ppm", ylab = "intensity", xlim = NULL, ...) {

    ## Check that the input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    metabo_matrix = assays(metabo_SE)$metabolic_data
    ppm = rownames(metabo_SE)
    if (suppressWarnings(is.na(as.numeric(ppm[1])))) {
        stop("metabo_SE rownames seem not to correspond to a ppm scale")
    } else {
        ppm = as.numeric(ppm)
    }

    if (is.null(xlim)) {
        xlim = rev(range(ppm))
    } else {
        ind1 = which(ppm >= xlim[2])[1]
        ind2 = which(ppm >= xlim[1])[1]
        ppm = ppm[ind1:ind2]
        metabo_matrix = metabo_matrix[ind1:ind2, ]
    }

    matplot(ppm, metabo_matrix, type = type, xlim = xlim, xlab = xlab,
        ylab = ylab, lty = lty, ...)
}

