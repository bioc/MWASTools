## EXTERNAL FUNCTIONS##

### QC_PCA ##

QC_PCA = function(metabo_SE, scale = FALSE, center = TRUE, ...) {

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }

    metabo_matrix = assays(metabo_SE)$metabolic_data

    if (any(is.na(metabo_matrix))) {
        warning("metabolic variables with NA values were removed")
        metabo_matrix = na.omit(metabo_matrix)
    }

    metabo_matrix = t(metabo_matrix)
    model = prcomp(metabo_matrix, scale = scale, center = center, ...)

    return(model)
}

### QC_PCA_scoreplot ##

QC_PCA_scoreplot = function(PCA_model, metabo_SE, plot_labels = FALSE,
    px = 1, py = 2, CI_level = 0.95, pch = 20, xlim = NULL, ylim = NULL,
    color_scale = c("cornflowerblue", "red"), grid = TRUE, ...) {

    ## Get scores
    PCA_scores = PCA_model$x

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }

    sample_ids = colnames(metabo_SE)
    class = colData(metabo_SE)[, "sample_type"]

    if (!is.null(ylim)) {
        if (length(ylim) != 2 | !is.numeric(ylim)) {
            stop("Invalid ylim: ylim must be a numeric vector with length 2")
        }
    }

    if (!is.null(xlim)) {
        if (length(xlim) != 2 | !is.numeric(xlim)) {
            stop("Invalid xlim: xlim must be a numeric vector with length 2")
        }
    }

    if (!is.numeric(CI_level)) {
        stop("CI_level must be a numeric values")
    }

    ## Get variance explained by each PC
    eig = (PCA_model$sdev)^2
    variance = eig * 100/sum(eig)
    cumvar = cumsum(variance)
    variance_PCA = data.frame(eig = eig, variance = variance, cumvariance = cumvar)
    variance_px = round(variance_PCA[px, 2], 1)
    variance_py = round(variance_PCA[py, 2], 1)

    ## Get PCs of interest
    PC_x = PCA_scores[, px]
    PC_y = PCA_scores[, py]

    ## Get plot limits
    if ( 0 %in% class) {  ## Remove QC to calculate ellipse
        index_experimental = which(class == 0)

        ellipse = dataEllipse(PC_x[index_experimental], PC_y[index_experimental],
                              levels = c(CI_level), draw = FALSE)
    } else { # e.g. do PCA of QC samples
        ellipse = dataEllipse(PC_x, PC_y, levels = c(CI_level), draw = FALSE)
    }

    if (is.null(xlim)) {
        xlim = c(min(min(PC_x), min(ellipse[, 1])), max(max(PC_x),
            max(ellipse[, 1])))
    }
    if (is.null(ylim)) {
        ylim = c(min(min(PC_y), min(ellipse[, 2])), max(max(PC_y),
            max(ellipse[, 2])))
    }

    color_vector = as.character(class)
    color_vector[color_vector == 0] = color_scale[1]
    color_vector[color_vector == 1] = color_scale[2]

    xlab = paste(paste("PC", px, "(", sep = ""), paste(variance_px, "%)", sep = ""),
                 sep = "")
    ylab = paste(paste("PC", py, "(", sep = ""), paste(variance_py, "%)", sep = ""),
                 sep = "")

    ## Plot scores
    if (pch == 21 | pch == 25) {
        bg = color_vector
        plot(PC_x, PC_y, pch = pch, bg = color_vector, xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...)
    } else {
        plot(PC_x, PC_y, pch = pch, col = color_vector, xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...)
    }

    ## remove QC to calculate ellipse
    if ( 0 %in% class) {
      index_experimental = which(class == 0)
      ellipse = dataEllipse(PC_x[index_experimental], PC_y[index_experimental],
                            levels = c(CI_level), add = TRUE, col = "black",
                            lwd = 0.6, plot.points = FALSE, center.cex = 0.2,
                            center.pch = NULL)
    } else { # e.g. do PCA of QC samples
      ellipse = dataEllipse(PC_x, PC_y, levels = c(CI_level), add = TRUE, col = "black",
                            lwd = 0.6, plot.points = FALSE, center.cex = 0.2,
                            center.pch = NULL)
    }

    if (plot_labels == TRUE) {
        text(PC_x, PC_y, labels = sample_ids, cex = 0.7)
    }

    if (grid == TRUE) {
        grid(lwd = 1.5, col = "gray", lty = 3, equilogs = TRUE)
    }
}
