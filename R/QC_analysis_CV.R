# INTERNAL FUNCTIONS#
CV_calculation = function(vector) {
    value = sd(vector)/mean(vector)
    return(value)
}

CV_features = function(vector, CV_th) {
    features_15 = round(sum(vector < 0.5 * CV_th)/length(vector), 2)
    features_30 = round(sum(vector < CV_th)/length(vector), 2)
    all_features = c(features_15, features_30)
    return(all_features)
}

# EXTERNAL FUNCTIONS#

## QC_CV##
QC_CV = function(metabo_SE, CV_th = 0.3, plot_hist = TRUE, hist_bw = 0.005,
    hist_col = "moccasin", size_lab = 12, size_axis = 12) {

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    if ((1 %in% metabo_SE$sample_type) == FALSE) {
        stop("metabo_SE does not contain QC samples (coded as sample_type = 1)")
    }

    QCmetabo_matrix = assays(metabo_SE[, metabo_SE$sample_type == 1])$metabolic_data
    metabo_ids = rownames(QCmetabo_matrix)

    ## Report number of variables
    to_print = paste(ncol(QCmetabo_matrix), "QC samples were detected", sep = " ")
    message(to_print)

    cols_metabo = split(QCmetabo_matrix, row(QCmetabo_matrix))
    CV_metabo = sapply(cols_metabo, CV_calculation)
    CV_metabo = abs(CV_metabo)  # absolute value of CV

    names(CV_metabo) = metabo_ids

    features_CV_metabo = CV_features(CV_metabo, CV_th = CV_th)

    if (is.na(features_CV_metabo[1])) {
        warning("CV calculation failed for at least some features")
        return(CV_metabo)
    } else {
        message("CV summary:")
        to_print = paste("   % metabolite features with CV <", 0.5 *
                        CV_th, ":", 100 * features_CV_metabo[1], sep = " ")
        message(to_print)
        to_print = paste("   % metabolite features with CV <", CV_th,
                         ":", 100 * features_CV_metabo[2], sep = " ")
        message(to_print, "\n")

        ## Set max CV = 1
        CV_metabo_backup = CV_metabo

        if (plot_hist == TRUE) {
            CV_metabo[CV_metabo > 1] = 1
            figure = ggplot(as.data.frame(CV_metabo), aes(CV_metabo)) +
              geom_histogram(fill = hist_col[1], binwidth = hist_bw,
                             colour = "black") +
              theme_bw() + labs(x = "CV", y = "count") +
              theme(panel.grid.major = element_blank(), panel.grid.minor =
                      element_blank(), axis.text = element_text(size = size_axis),
                    axis.title = element_text(size = size_lab), axis.title.y =
                      element_text(vjust = 0), axis.title.x = element_text(vjust = 0))
            plot(figure)
      }
        return(CV_metabo_backup)
    }

}

## QC_CV_specNMR ##

QC_CV_specNMR = function(metabo_SE, ref_sample, CV_th = 0.3,
    xlab = "ppm", ylab = "intensity", size_axis = 12, size_lab = 12,
    xlim = NULL, ylim = NULL, xbreaks = waiver(), xnames = waiver(),
    ybreaks = waiver(), ynames = waiver()) {

    ## Check that the input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }

    ppm = rownames(metabo_SE)
    if (suppressWarnings(is.na(as.numeric(ppm[1])))) {
        stop("metabo_SE rownames seem not to correspond to a ppm scale")
    } else {
        ppm = as.numeric(ppm)
    }

    ## Sort xlim in decreasing order (ppm is plotted in inverse scale)
    xlim = suppressWarnings(sort(xlim, decreasing = TRUE))

    if(!is.null(xlim)) {
        if(ppm[1] < tail(ppm, n = 1)) {
            ind1 = tail(which(ppm <= min(xlim)), n = 1)
            ind2 = which(ppm >= max(xlim))[1]
        } else {
            ind1 = which(ppm <= min(xlim))[1]
            ind2 = tail(which(ppm <= max(xlim)), n = 1)
        }
        if (length(ind1) == 0 | length(ind2) == 0) {
            stop("xlim is not contained in metabo_SE rownames")
        }
        if (is.na(ind1) | is.na(ind2)) {
            stop("xlim is not contained in metabo_SE rownames")
        }
        ind_range = sort(ind1:ind2)
        metabo_SE = metabo_SE[ind_range, ]
        ppm = ppm[ind_range]
    }

    metabo_matrix = t(assays(metabo_SE)$metabolic_data)
    if (ref_sample %in% colnames(metabo_SE) == FALSE) {
        stop ("ref_sample is not included in metabo_SE colnames")
    } else {
        metabo_vector = metabo_matrix[ref_sample, ]
    }

    if (is.null(ylim) & !is.null(xlim)) {
        if (min(metabo_vector) > 0) {
            ylim = c(0.90*(min(metabo_vector)),
                     1.10*(max(metabo_vector)))
        } else {
            ylim = c(1.15*(min(metabo_vector)),
                     1.10*(max(metabo_vector)))
      }
    }

    ## Calculate CV_metabo
    CV_metabo = QC_CV (metabo_SE, CV_th = CV_th, plot_hist = FALSE)

    ## Create a NMR spectrum colored by CV
    abs.CV = CV_metabo
    abs.CV[abs.CV >= CV_th] = CV_th

    data_CV = data.frame(ppm = ppm, metabo_vector = metabo_vector, abs.CV = abs.CV,
                         CV_raw = CV_metabo)
    # color_scale = c('dodgerblue2', 'green3', 'gold',
    # 'darkorange2','orangered2', 'red1')
    color_scale = c("green3", "dodgerblue2", "plum2", "purple", "purple4", "red1")
    col_values = c(0, 0.45, 0.55, 0.9, 0.9996666, 1)

    figure_spectrum = ggplot(data_CV, aes(ppm, metabo_vector,
        color = abs.CV)) + geom_line() + scale_colour_gradientn(colours = color_scale,
        values = col_values, space = "Lab", limits = c(0, CV_th),
        breaks = c(0, CV_th/2, CV_th)) + scale_x_reverse(limits = xlim,
        breaks = xbreaks, labels = xnames) + scale_y_continuous(limits = ylim,
        breaks = ybreaks, labels = ynames) + theme_bw() + labs(x = xlab,
        y = ylab) + theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text = element_text(size = size_axis),
        axis.title = element_text(size = size_lab, vjust = 0))

    plot(figure_spectrum)

    return(figure_spectrum)
}

## QC_CV_scatterplot ##

QC_CV_scatterplot = function(rt, mz, CV_metabo, CV_th = 0.30, xlab = "rt",
                             ylab = "mz", pch = 20, marker_size = 1, xlim = NULL,
                             ylim = NULL, size_axis = 10, size_lab = 10) {

    ## Check if input variables are correct
    if (is.vector(rt) == FALSE | is.numeric(rt) == FALSE) {
        stop ("rt must be a numeric vector")
    }

    if (is.vector(mz) == FALSE | is.numeric(mz) == FALSE) {
        stop ("mz must be a numeric vector")
    }

    if (is.vector(CV_metabo) == FALSE | is.numeric(CV_metabo) == FALSE) {
        stop ("CV_metabo must be a numeric vector")
    }

    if (length(CV_metabo) != length(rt) | length(rt) != length(mz)) {
        stop("lengths of rt, mz, and CV_metabo must be consistent")
    }

    CV_metabo = abs(CV_metabo)
    CV_metabo[CV_metabo >= CV_th] = CV_th

    data_CV = data.frame(rt = rt, mz = mz, abs.CV = CV_metabo)

    color_scale = c("green3", "dodgerblue2", "plum2", "purple", "purple4", "red1")
    col_values = c(0, 0.45, 0.55, 0.9, 0.9996666, 1)

    figure = ggplot (data_CV, aes(rt, mz, color = abs.CV)) +
        geom_point(shape = pch, size = marker_size) +
        scale_colour_gradientn(colours = color_scale, values = col_values,
                               space = "Lab", limits = c(0, CV_th),
                               breaks = c(0, CV_th/2, CV_th)) +
        theme_bw() +
        labs(x = xlab, y = ylab) +
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(limits = ylim) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text = element_text(size = size_axis),
              axis.title = element_text(size = size_lab, vjust = 0))

    return(figure)
}


## CV_filter ##
CV_filter = function(metabo_SE, CV_metabo, CV_th = 0.3) {

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }

    if (!is.vector(CV_metabo)) {
        stop("CV_metabo needs to be a numeric")
    }
    if (sum(is.na(CV_metabo)) > 0) {
        stop ("Cannot perform CV filtering with NA values in CV_metabo")
    }

    if (nrow(metabo_SE) != length(CV_metabo)) {
        stop("CV_metabo length must be consistent with metabolic_data dimension")
    }

    if (is.numeric(CV_th) == FALSE) {
        stop("CV_th must be a numeric value")
    }

    index_wanted = which(CV_metabo < CV_th)

    if (length(index_wanted) > 0) {
        metabo_SE = metabo_SE[index_wanted, ]
        return(metabo_SE)
    } else {
        stop("None of the metabolic features meets the filtering criteria")
    }
}

