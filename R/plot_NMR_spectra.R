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
        if(ppm[1] < tail(ppm, n = 1)) {
            ind1 = tail(which(ppm <= min(xlim)), n = 1)
            ind2 = which(ppm >= max(xlim))[1]
        } else {
            ind1 = which(ppm <= min(xlim))[1]
            ind2 = tail(which(ppm <= max(xlim)), n = 1)
        }
        if (length(ind1) == 0  | length(ind2) == 0) {
            stop("xlim is not contained in metabo_SE rownames")
        }
        if (is.na(ind1) | is.na(ind2)) {
            stop("xlim is not contained in metabo_SE rownames")
        }
        ppm = ppm[ind1:ind2]
        metabo_matrix = metabo_matrix[ind1:ind2, ]

        if (is.null(ylim)) {
          if (min(metabo_matrix) > 0 ) {
              ylim = c(0.90*(min(metabo_matrix)),
                       1.10*(max(metabo_matrix)))
          } else {
              ylim = c(1.15*(min(metabo_matrix)),
                       1.10*(max(metabo_matrix)))
          }

        }
    }

    matplot(ppm, metabo_matrix, type = type, xlim = xlim, xlab = xlab,
        ylab = ylab, lty = lty, ...)
}

