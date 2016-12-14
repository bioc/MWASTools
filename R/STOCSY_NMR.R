## INTERNAL FUNCTIONS ##

### get_pval##
get_pval = function(x, y, method = method) {
    pval = cor.test(x, y, method = method)$p.val
    return(pval)
}

### BH_correct ##
BH_correct = function(x, y, alpha_th) {
    if (y > alpha_th) {
        x = 0
    }
    return(x)
}

## EXTERNAL FUNCTION ##

### STOCSY_NMR ##

STOCSY_NMR = function(metabo_SE, ppm_query, alpha_th = 0.05,
                      xlab = "ppm", ylab = "covariance", size_lab = 12, size_axis = 12,
                      xlim = NULL, ylim = NULL, xbreaks = waiver(), xnames = waiver(),
                      ynames = waiver(), ybreaks = waiver()) {

    ## Check that the input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    metabo_matrix = t(assays(metabo_SE)$metabolic_data)
    ppm = rownames(metabo_SE)
    if (suppressWarnings(is.na(as.numeric(ppm[1])))) {
      stop("metabo_SE rownames seem not to correspond to a ppm scale")
    } else {
      ppm = as.numeric(ppm)
    }

    ppm_query = ppm_query[1]  # in case the user enters more than 1 value
    if (identical(ppm_query, round(ppm_query, 1))) {
        stop("ppm_query must have at least two decimals")
    }

    ppm_index = grep(ppm_query, ppm, fixed = TRUE)[1]
    if (length(ppm_index) == 0 | is.na(ppm_index) == TRUE) {
        stop("Invalid ppm_query: make sure that ppm_query is contained in ppm")
    }

    ## Sort xlim in decreasing order (ppm is plotted in inverse scale)
    xlim = suppressWarnings(sort(xlim, decreasing = TRUE))

    ## Run STOCSY
    driver = metabo_matrix[, ppm_index] # NMR driver signal

    cols_metabo = split(t(metabo_matrix), row(t(metabo_matrix)))

    cov_metabo = as.vector(cov(driver, metabo_matrix, method = "pearson"))
    cor_metabo = as.vector(cor(driver, metabo_matrix, method = "pearson"))

    pval_metabo = sapply(cols_metabo, get_pval, y = driver, method = "pearson")
    BH_pvalue = p.adjust(pval_metabo, method = "BH")
    RT = list()
    RT[["alpha_th"]] = alpha_th
    cov_metabo_adjusted = mapply(BH_correct, cov_metabo, BH_pvalue, MoreArgs = RT,
                                 SIMPLIFY = FALSE)
    covar = unlist(cov_metabo_adjusted)
    cor_metabo_adjusted = mapply(BH_correct, cor_metabo, BH_pvalue, MoreArgs = RT,
                                 SIMPLIFY = FALSE)
    abs.r = abs(unlist(cor_metabo_adjusted))

    if(!is.null(xlim) & is.null(ylim)) {
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
        ylim = c(0.90*(min(covar[ind1:ind2])),
                 1.10*(max(covar[ind1:ind2])))
    }

    data_cov = data.frame(ppm = ppm, covar = covar, abs.r = abs.r)

    col_scale = c("red4", "red", "orangered", "darkorange4", "darkolivegreen",
                  "darkgreen", "seagreen", "cyan1", "dodgerblue", "darkblue")
    col_values = sort(seq(0, 1, 0.1), decreasing = TRUE)

    ## Plot STOCSY
    figure_STOCSY = ggplot(data_cov, aes(ppm, covar, color = abs.r)) +
                    geom_line() +
                    scale_colour_gradientn(colours = col_scale,
                                           values = col_values, space = "Lab",
                                           limits = c(0, 1)) +
                    theme_bw() +
                    scale_x_reverse(limits = xlim, breaks = xbreaks, labels = xnames) +
                    scale_y_continuous(limits = ylim, breaks = ybreaks, labels = ynames) +
                    labs(x = xlab, y = ylab) +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.text = element_text(size = size_axis),
                          axis.title = element_text(size = size_lab, vjust = 0))

    plot(figure_STOCSY)
    return(figure_STOCSY)
}

