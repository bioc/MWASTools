#### Clustering method ######

## User Parameters ##
# st: minimum cluster size
# ct: correlation threshold
# int: integration method (sum, median or max)
# ct: correlation method ("pearson" or "spearman")


## INTERNAL FUNCTIONS ##

## roundUp ##
roundUp = function(x,to = 1) {
    to*(x%/%to + as.logical(x%%to))
}

## get_mean ##
get_mean = function(cluster, cm) {
    cor_matrix = cor(cluster, cluster,
                     method = cm)
    mean_cor = mean(cor_matrix[upper.tri(cor_matrix)])
    return(mean_cor)
}

################################################################################

## find_max ##
find_max = function(peak_range, NMR_data, JBA_clusters, mean_clusters, th_value, cor_values) {
    ## Find other clusters (i.e. local maxima) in the peak range
    ## A given variable (i.e. mean cluster) is consider as a maximum when there are at least
    ## n_th/2 descending variables on the right and n_th/2 ascending variables on the left.

    for (i in 1:length(peak_range)){

        ind_in_clusters = unlist(strsplit(JBA_clusters[, "all_ind"], "_"))
        # list of index included in the clusters so far

        if(peak_range[i] - th_value <= 0) { # we are beggining of the ppm_scale
            cor_values_lower = cor_values[peak_range[i]]
        } else {
            cor_values_lower = cor_values[(peak_range[i] - th_value):peak_range[i]]
        }
        if(peak_range[i] + th_value >= ncol(NMR_data)) { # we are at the end of the ppm_scale
            cor_values_upper = cor_values[peak_range[i]]
        } else {
            cor_values_upper = cor_values[peak_range[i]:(peak_range[i] + th_value)]
        }

        ans_upper = identical(cor_values_upper,
                              sort(cor_values_upper, decreasing = T))
        ans_lower = identical(cor_values_lower,
                              sort(cor_values_lower, decreasing = F))

        ind_clust = unlist(strsplit(mean_clusters[peak_range[i],"all_ind" ], "_"))
        ans_clust = length(intersect(ind_in_clusters, ind_clust)) == 0

        if (ans_upper & ans_lower & ans_clust) {
            JBA_clusters = rbind(JBA_clusters, mean_clusters[peak_range[i], ])
        }
    }
    return(JBA_clusters)
}

################################################################################

## find_min ##
find_min = function(peak_range, JBA_clusters, mean_clusters, th_value, cor_values) {
    ## Find other clusters (i.e. local maxima) in the peak range
    ## To find the local maxima, we search for local min. A given variable
    ## (i.e. mean cluster) is consider as a minimum when there are at least
    ## n_th/2 descending variables on the right and n_th/2 ascending variables on the left.

    start = th_value + 1
    end = length(peak_range) - (th_value + 1)

    min_ans = rep(0, length(peak_range))

    for (i in start:end) {
        cor_values_lower = cor_values[(peak_range[i] - th_value):peak_range[i]]
        cor_values_upper = cor_values[peak_range[i]:(peak_range[i] + th_value)]
        ans_lower = identical(cor_values_lower,
                              sort(cor_values_lower, decreasing = T))
        ans_upper = identical(cor_values_upper,
                              sort(cor_values_upper, decreasing = F))

        if( ans_lower & ans_upper) {
            min_ans[i] = 1
        }
    }
    if (sum(min_ans == 1) > 0) { # found some min
        min_ans[1] = 1 # the first ind of the peak range
        min_ans[length(min_ans)] = 1 # the last ind of the peak range

        all_min_ind = which(min_ans == 1)

        for (x in 1:length(all_min_ind[-1])) {
            ind_start = all_min_ind[x]
            ind_end = all_min_ind[x+1]
            search_range = peak_range[ind_start:ind_end]
            ind_max = search_range[which(cor_values[search_range] == max(cor_values[search_range]))[1]]
            JBA_clusters = rbind(JBA_clusters, mean_clusters[ind_max, ])
        }
    }
    return(JBA_clusters)
}

################################################################################

## set_new_limits ##
set_new_limits = function(core_cluster, upper_limits, lower_limits, NMR_data,
                          all_ind, cm = "pearson", cov = FALSE, ct = ct,
                          dup = FALSE) {
    new_limits = c()
    value = min(length(upper_limits), length(lower_limits))
    # in case we are at the beggining or at the end of the ppm scale

    for(z in 1:value) {
        idxu = upper_limits[z] ## upper index
        idxl = lower_limits[z] ## lower index

        ## Check covariance or correlation goes in the right direction
        if (idxu < ncol(NMR_data)) {
            if (cov) { # Check the covariance
                cov_ansu = cov(NMR_data[, idxu], NMR_data[, idxu - 1]) >=
                    cov(NMR_data[, idxu], NMR_data[, idxu + 1])
            } else { # Check the correlation - round to 3, avoids cutting peaks when the cor spec is flat
                #at the top
                cov_ansu = round(cor(NMR_data[, idxu], NMR_data[, idxu - 1], method = cm), 3) >=
                    round(cor(NMR_data[, idxu], NMR_data[, idxu + 1], method = cm), 3)
            }
        } else {
            cov_ansu = TRUE
        }

        if (idxl > 0) {
            if (cov) { # Check the covariance
                cov_ansl = cov(NMR_data[, idxl], NMR_data[, idxl - 1]) <=
                    cov(NMR_data[, idxl], NMR_data[, idxl + 1])
            } else { # Check the correlation
                cov_ansl = round(cor(NMR_data[, idxl], NMR_data[, idxl - 1], method = cm), 3) <=
                    round(cor(NMR_data[, idxl], NMR_data[, idxl + 1], method = cm), 3)
            }
        } else {
            cov_ansl = TRUE
        }

        ## Check that the lower or upper index is not included in another cluster
        if(dup == FALSE) {
            if (idxl %in% all_ind) {
                cov_ansl = FALSE
            }
            if (idxu %in% all_ind) {
                cov_ansu = FALSE
            }
        }

        ## Redefine cluster
        candidate_cluster_up = cbind(core_cluster, NMR_data[, idxu])
        candidate_cluster_low = cbind(core_cluster, NMR_data[, idxl])

        if(get_mean(candidate_cluster_up, cm) >= ct & cov_ansu) {
            ## Add one downfield variable
            core_cluster = candidate_cluster_up # update the cluster
            new_limits = c(new_limits, idxu)
            ## Now try to add one upfield variable
            candidate_cluster_low = cbind(core_cluster, NMR_data[, idxl])
            if (get_mean(candidate_cluster_low, cm) >= ct & cov_ansl) {
                core_cluster = candidate_cluster_low # update the cluster
                new_limits = c(new_limits, idxl)
            } else {
                lower_limits = rep(lower_limits[z], length(lower_limits))
            }
        } else if (get_mean(candidate_cluster_low, cm) >= ct & cov_ansl) {
            ## Add one downfield variable
            core_cluster = candidate_cluster_low # update the cluster
            new_limits = c(new_limits, idxl)
            upper_limits = rep(upper_limits[z], length(upper_limits))
        } else {
            break
        }
    }
    return(new_limits)
}

################################################################################

## JBA_mergeClusters ## Merges highly correlated neighboring clusters
JBA_mergeClusters = function(NMR_JBA, NMR_data, cm = "pearson", mt = 0.90,
                             int = "sum") {

    # STEP 0a: Predifine parameters
    cov = FALSE # Merge clusters based on correlation (instead of covariance)

    # STEP 0b: Get data
    JBA_data = NMR_JBA$clustered_data
    JBA_clusters_exp = NMR_JBA$cluster_info
    rownames(JBA_clusters_exp) = JBA_clusters_exp[, 1]
    ppm = as.numeric(colnames(NMR_data))

    if(nrow(JBA_data) != nrow(NMR_data)) {
        stop("NMR_JBA and NMR_data have inconsistent dimensions")
    }
    mean_clusters = NMR_JBA$all_clusters

    # STEP 1: Decide which clusters can be merged

    ## Get correlation matrix of the JBA clusters
    cor_matrix = cor(JBA_data, JBA_data, method = cm)
    cor_pairs = matrix(nrow = ncol(cor_matrix)-1, ncol = 3)
    rownames(cor_pairs) = paste(colnames(cor_matrix)[-ncol(cor_matrix)],
                                colnames(cor_matrix)[seq(2,ncol(cor_matrix)-1, 1)],
                                sep = "_")
    colnames(cor_pairs) = c("cor", "limits", "decision")

    ## A given cluster in JBA_data will be merge with the next one (a + 1) if:
    ## cor(a, a+1) > ct; and they are actually adjacent based on JBA_clusters_exp.

    for (a in 1:(ncol(cor_matrix)-1)) {
        cor_pairs[a, 1] = as.numeric(cor_matrix[a, a+1])
        low_lim = as.numeric(unlist(strsplit(JBA_clusters_exp[a, "ind_limits"], "_")))
        up_lim = as.numeric(unlist(strsplit(JBA_clusters_exp[a +1, "ind_limits"], "_")))
        cor_pairs[a,2] = paste(low_lim[1], up_lim[2], sep = "_")

        if(as.numeric(cor_pairs[a, 1]) > mt & (up_lim[1] - low_lim[2] == 1)) {
            cor_pairs[a,3] = "merge"
        }
    }

    ind_merge = which(cor_pairs[, 3] == "merge")

    if(length(ind_merge) == 0) {
        message = "None of the clusters meets the criteria to be merged"
        print(message)
        return(NULL)
    }

    cor_pairs_cleaned = na.omit(cor_pairs)
    cor_pairs_cleaned = matrix(cor_pairs_cleaned, ncol = 3)
    rownames(cor_pairs_cleaned) = rownames(cor_pairs)[ind_merge]
    colnames(cor_pairs_cleaned) = colnames(cor_pairs)
    cluster_ids = unlist(strsplit(rownames(cor_pairs_cleaned), "_"))

    # Remove clusters that will be merged from original datasets
    JBA_data_merged = JBA_data[, setdiff(colnames(JBA_data), cluster_ids)]
    JBA_clusters_merged = JBA_clusters_exp[setdiff(colnames(JBA_data), cluster_ids), ]

    dup_idx = which(duplicated(cluster_ids))
    cor_pairs_final = cor_pairs_cleaned

    if (length(dup_idx) > 0) { ## More than 2 consecutive clusters will be merged
      for (idx in dup_idx) {
        dup_clusters = grep(cluster_ids[idx], rownames(cor_pairs_final), fixed = TRUE)
        if (length(dup_clusters) > 0) { ## This has been updated # Jun 2018
          clust_range = unlist(strsplit(cor_pairs_final[dup_clusters, 2], "_"))
          clust_lim = paste(min(clust_range), max(clust_range), sep = "_")
          rownames(cor_pairs_final)[dup_clusters[1]] =
            paste(rownames(cor_pairs_final)[dup_clusters], collapse = "_")
          cor_pairs_final[dup_clusters[1], 2] = clust_lim ## updates the merging limits
          cor_pairs_final = cor_pairs_final[-dup_clusters[-1], ]
        }
      }
    }

    # STEP 2: Merge selected clusters
    merged_clusters = matrix(ncol = nrow(cor_pairs_final), nrow = nrow(JBA_data))
    merged_clusters_info = matrix(nrow = nrow(cor_pairs_final), ncol = ncol(JBA_clusters_merged))
    colnames(merged_clusters_info) = colnames(JBA_clusters_merged)
    rownames(merged_clusters_info) = rownames(cor_pairs_final)
    colnames(merged_clusters) = rownames(cor_pairs_final)

    for (y in 1:nrow(cor_pairs_final)) {
        lim_idx = as.numeric(unlist(strsplit(cor_pairs_final[y, 2], "_")))
        merged_clusters_info[y, "mean_ppm"] = mean(ppm[lim_idx[1]:lim_idx[2]])
        merged_clusters_info[y, "mean_cor"] = get_mean(NMR_data[, lim_idx[1]:lim_idx[2]],
                                                       cm = cm)
        merged_clusters_info[y, "all_ppm"] = paste(ppm[lim_idx[1]:lim_idx[2]], collapse = "_")
        merged_clusters_info[y, "all_ind"] = paste(lim_idx[1]:lim_idx[2], collapse = "_")
        merged_clusters_info[y, "ppm_limits"] = paste(ppm[lim_idx[1]], ppm[lim_idx[2]], sep = "_")
        merged_clusters_info[y, "ind_limits"] = paste(lim_idx, collapse = "_")
        rownames(merged_clusters_info)[y] = as.character(merged_clusters_info[y, "mean_ppm"])

        colnames(merged_clusters)[y] = rownames(merged_clusters_info)[y]
        if(int == "max") {
            merged_clusters[, y] = apply(NMR_data[, lim_idx[1]:lim_idx[2]], 1, max)
        } else if (int == "median") {
            merged_clusters[, y] = apply(NMR_data[, lim_idx[1]:lim_idx[2]], 1, median)
        } else if (int == "mean") {
            merged_clusters[, y] = apply(NMR_data[, lim_idx[1]:lim_idx[2]], 1, mean)
        } else {
            merged_clusters[, y]  = apply(NMR_data[, lim_idx[1]:lim_idx[2]], 1, sum)
        }

    }

    # STEP 3: Add the merged clusters to the JBA_data_merged
    JBA_data_merged = cbind(JBA_data_merged, merged_clusters)
    JBA_data_merged = JBA_data_merged[, sort(colnames(JBA_data_merged), decreasing = F)]

    JBA_clusters_merged = rbind(JBA_clusters_merged, merged_clusters_info)
    JBA_clusters_merged = JBA_clusters_merged[sort(rownames(JBA_clusters_merged), decreasing = F), ]
    rownames(JBA_clusters_merged) = NULL

    res = vector(mode = "list", length = 2)
    names(res) = c("merged_clusters", "JBA_data_merged")
    res[[1]] = cor_pairs_cleaned ## all clusters to be merged
    res[[2]] = NMR_JBA$all_clusters ## all clusters in the input dataset
    res[[3]] = NMR_JBA$metabo_clusters ## all metabo in the input dataset
    res[[4]] = JBA_clusters_merged  ## merged cluster info
    res[[5]] = JBA_data_merged ## merged JBA data
    colnames(JBA_data_merged) = round(as.numeric(colnames(JBA_data_merged)), 3)
    res[[6]]  = JBA_data_merged # merged JBA data with rounded colnames
    names(res) = c("cor_pairs", "all_clusters","metabo_clusters", "cluster_info",
                   "clustered_data", "clustered_data_names_rounded")
    return(res)
}

################################################################################
################################################################################

## EXTERNAL FUNCTIONS ##

## JBA_corDistribution ##
JBA_corDistribution = function(NMR_data, st = 4, cm = "pearson",
                               metabo_range = c(3.50, 3.96),
                               noise_range = c(9.72, 9.99),
                               color_scale = c("lightcoral", "honeydew3")) {
    ## Check st value
    st = round(st, 0)
    if (st < 2) {
        stop ("st value should be at least 2")
    }

    ## Get st clusters
    mean_clusters = matrix(ncol = 4, nrow = (ncol(NMR_data) - st))
    ppm = as.numeric(colnames(NMR_data))

    for( i in 1:(ncol(NMR_data) - st)) {
        range = c(i:(i+ (st - 1)))
        mean_cor = get_mean(NMR_data[, range], cm = cm)
        mean_ppm = mean(ppm[range])
        all_ppm = paste(ppm[range], collapse = "_")
        mean_clusters[i, ] = c(mean_ppm, mean_cor, all_ppm, paste(range, collapse = "_"))
    }

    ## Get metabo and noise ranges
    metabo_range = sort(metabo_range, decreasing = FALSE)
    noise_range = sort(noise_range, decreasing = FALSE)

    ppm_clusters = as.numeric(mean_clusters[, 1])

    if (noise_range[2] > tail(ppm_clusters, 1)| noise_range[1] < ppm_clusters[1]) {
        stop("noise_range not included in ppm scale")
    }

    metabo_ppm = which(ppm_clusters >= metabo_range[1] & ppm_clusters <= metabo_range[2])
    noise_ppm = which(ppm_clusters >= noise_range[1] & ppm_clusters <= noise_range[2])

    ## Compare distributions of r coefficients ##
    metabo = as.numeric(mean_clusters[metabo_ppm, 2])
    metabo_fun = ecdf(metabo)
    noise = as.numeric(mean_clusters[noise_ppm, 2])
    noise_fun = ecdf(noise)

    metabo_cum = data.frame(r.coeff = metabo, cum = as.numeric(metabo_fun(metabo)))
    noise_cum = data.frame(r.coeff = noise, cum = as.numeric(noise_fun(noise)))
    noise_sorted = noise_cum[order(noise_cum[, "r.coeff"]), ]
    ct_val = noise_sorted[which(noise_sorted[, "cum"] == 1)[1], ]

    ylab = paste("Cumulative proportion of clusters (st = ", st, ")", sep = "")
    ylab = "Cumulative proportion of clusters"

    plot = ggplot(metabo_cum, aes(r.coeff, cum)) +
        geom_area(fill  = color_scale[1], alpha = 0.5, color = color_scale[1], size = 0.7) +
        geom_area(data = noise_cum, aes(r.coeff, cum), fill = color_scale[2],
                  color = color_scale[2], alpha = 0.6, size = 0.7) +
        labs(x = "r coefficient", y = ylab) + theme_bw() +
        theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10, vjust = 0),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(angle = 0, size = 10))

    ## Suggest st value
    print(paste("Suggested ct value:", round(ct_val[, 1], 5)))

    return(plot)
}

## JBA_binning ##
JBA_binning = function(NMR_data, st = 4, ct = 0.85, int = "sum",
                       cm = "pearson", ef = 2, merge = TRUE, mt = 0.9) {

    # STEP 0: Predefined parameters & check data
    search_max = TRUE # Find local maxima (instead of local minima)
    dup = FALSE # Not allow the same variable to be included in several clusters
    cov = FALSE # Expand the clusters based on correlation (not covariance)

    if (!is.matrix(NMR_data)| !is.numeric(NMR_data)) {
        stop ("NMR_data must be a numeric matrix")
    }

    # STEP 1: Get all possible mean clusters of size st

    mean_clusters = matrix(ncol = 4, nrow = (ncol(NMR_data) - st))
    ppm = as.numeric(colnames(NMR_data))

    for(i in 1:(ncol(NMR_data) - st)) {
        range = c(i:(i+ (st - 1)))
        mean_cor = get_mean(NMR_data[, range], cm = cm)
        mean_ppm = mean(ppm[range])
        all_ppm = paste(ppm[range], collapse = "_")
        mean_clusters[i, ] = c(mean_ppm, mean_cor, all_ppm, paste(range, collapse = "_"))
    }
    colnames(mean_clusters) = c("mean_ppm", "mean_cor", "all_ppm", "all_ind")

    ###########################################################################

    # STEP2: Identify the candidate clusters: i.e. mean clusters passing the ct threshold

    cor_values = as.numeric(mean_clusters[, "mean_cor"])
    cor_values_clean = as.numeric(mean_clusters[, "mean_cor"])
    cor_values_clean[cor_values_clean < ct] = 0 # cor_values_mean below st are set to 0
    index_zero = which(cor_values_clean == 0)
    index_zero = unique(c(1, index_zero)) # in case there is not a 0 at the beggining

    JBA_clusters = c()
    for(index in index_zero) {
        if(index < length(cor_values_clean)) {
            if(cor_values_clean[index + 1] != 0) {# it found a peak
                #print(index)
                remaining_ind = (index + 1):length(cor_values_clean)
                ## Might give a problem if there is not a minimum in the last peak
                end_ind = remaining_ind[(which(cor_values_clean[remaining_ind] == 0)[1]) - 1]
                # found another zero: cluster ends
                if (is.na(end_ind)) { # the last cluster never drops to 0
                    end_ind = length(cor_values_clean)
                }
                peak_range = (index + 1):end_ind
                best_SRV_ind = peak_range[which(cor_values[peak_range] == max(cor_values[peak_range]))[1]]

                JBA_clusters = rbind(JBA_clusters, mean_clusters[best_SRV_ind, ])

                ## STEP 2.1. Find other possible clusters in the peak range
                th_value = roundUp(st/2)

                if (search_max) { ## Find local maxima
                    if (length(peak_range) > 1) {
                        JBA_clusters = find_max(peak_range = peak_range, NMR_data = NMR_data,
                                                JBA_clusters = JBA_clusters,
                                                th_value = th_value, cor_values = cor_values,
                                                mean_clusters = mean_clusters)
                    }
                } else { ## Find local minima (gives almost same results as local maxima)
                    if (length(peak_range) >= (2*th_value + 2)) {
                        JBA_clusters = find_min(peak_range = peak_range, JBA_clusters = JBA_clusters,
                                                th_value = th_value, cor_values = cor_values,
                                                mean_clusters = mean_clusters)
                    }
                }
            }

        }
    }

    JBA_clusters = unique(JBA_clusters)
    JBA_clusters = matrix(JBA_clusters, ncol = ncol(mean_clusters))
    JBA_clusters = JBA_clusters[order(JBA_clusters[, 1]), ]
    colnames(JBA_clusters) = colnames(mean_clusters)

    ###########################################################################

    # STEP3: Expand clusters and get JBA_data

    JBA_data = matrix(nrow = nrow(NMR_data), ncol = nrow(JBA_clusters))
    colnames(JBA_data) = JBA_clusters[, 1]
    JBA_clusters_exp = cbind(JBA_clusters, "ppm_limits" = NA, "ind_limits" = NA)
    st = st*ef # establishes the number of upfield and downfield variables
    # that can be added
    all_ind = unlist(strsplit(JBA_clusters[, "all_ind"], "_"))

    for(i in 1:ncol(JBA_data)) {
        ppm_limits = as.numeric(unlist(strsplit(JBA_clusters[i, "all_ind"], "_")))

        ## STEP 3.1: Expand clusters
        ## Expand the cluster if by adding other adjacent variables, the mean cor
        # of the cluster is still above the ct

        core_cluster = NMR_data[, ppm_limits]

        # Will allow to add max st variables from each side
        lower_limits = sort((ppm_limits[1] - st) : (ppm_limits[1] - 1), decreasing = T)
        lower_limits = lower_limits[lower_limits > 0]

        if(length(lower_limits) == 0) { # passes the lower limit of ppm
            #lower_limits = ppm_limits[1]
            lower_limits = 1
        }
        upper_limits = (tail(ppm_limits, 1) + 1): (tail(ppm_limits, 1) + st)
        upper_limits = upper_limits[upper_limits <= ncol(NMR_data)]

        if(length(upper_limits) == 0) { # passes the upper limit of ppm
            upper_limits = ncol(NMR_data)
        }

        new_limits = set_new_limits(core_cluster = core_cluster, upper_limits = upper_limits,
                                    lower_limits = lower_limits, NMR_data = NMR_data,
                                    all_ind = all_ind, cm = cm, cov = cov, dup = dup,
                                    ct = ct)

        all_ind = unique(c(all_ind, new_limits)) # This has been updated # Jan-18

        ## Update ppm_limits and JBA_clusters
        ppm_limits = sort(c(new_limits, ppm_limits))
        JBA_clusters_exp[i, "mean_ppm"] = mean(ppm[ppm_limits])
        JBA_clusters_exp[i, "mean_cor"] = get_mean(core_cluster, cm)
        JBA_clusters_exp[i, "all_ppm"] = paste(ppm[ppm_limits], collapse = "_")
        JBA_clusters_exp[i, "all_ind"] = paste(ppm_limits, collapse = "_")
        JBA_clusters_exp[i, "ppm_limits"] = paste(ppm[ppm_limits][1], tail(ppm[ppm_limits], 1),
                                                  sep = "_")
        JBA_clusters_exp[i, "ind_limits"] = paste(ppm_limits[1], tail(ppm_limits, 1),
                                                  sep = "_")
        colnames(JBA_data)[i] = mean(ppm[ppm_limits])

        ## STEP 3.2: Calculate integrals
        if(int == "max") {
            JBA_data[, i] = apply(NMR_data[, ppm_limits], 1, max)
        } else if (int == "median") {
            JBA_data[, i] = apply(NMR_data[, ppm_limits], 1, median)
        } else if (int == "mean") {
            JBA_data[, i] = apply(NMR_data[, ppm_limits], 1, mean)
        } else {
            JBA_data[, i] = rowSums(NMR_data[, ppm_limits])
        }
    }

    ###########################################################################

    # STEP 4: Generate output
    res = vector(mode = "list", length = 5)
    mean_clusters[, "mean_ppm"] = format(round(as.numeric(mean_clusters[, "mean_ppm"]), 7),
                           nsmall = 7)
    JBA_clusters[, "mean_ppm"] = format(round(as.numeric(JBA_clusters[, "mean_ppm"]), 7),
                                        nsmall = 7)
    JBA_clusters_exp[, "mean_ppm"] = format(round(as.numeric(JBA_clusters_exp[, "mean_ppm"]), 7),
                                            nsmall = 7)
    colnames(JBA_data) = format(round(as.numeric(colnames(JBA_data)), 7),
                                nsmall = 7)

    res[[1]] = mean_clusters ## all clusters formed
    res[[2]] = JBA_clusters ## best clusters of each peak
    res[[3]] = JBA_clusters_exp ## clusters expanded
    rownames(JBA_data) = rownames(NMR_data)
    res[[4]] = JBA_data
    colnames(JBA_data) = round(as.numeric(colnames(JBA_data)), 4)
    res[[5]] = JBA_data # SRV data with rounded colnames
    names(res) = c("all_clusters", "metabo_clusters", "cluster_info",
                   "clustered_data", "clustered_data_names_rounded")

    if (merge == FALSE) {
        final_res = res[-4]
        names(final_res) = c("all_clusters", "JBA_seeds", "JBA_bins_expanded",
                             "JBA_data") #Updated in 08/18
        return(final_res)
    }
    # Merge adjacent clusters with intercluster correlation > mt
    merged_data = JBA_mergeClusters(NMR_JBA = res, NMR_data = NMR_data, cm = cm,
                                    mt = mt, int = int)

    if (is.null(merged_data)) { ## Did not find clusters to be merged
        final_res = res[-4]
        names(final_res) = c("all_clusters", "JBA_seeds", "JBA_bins_expanded",
                             "JBA_data") #Updated in 08/18
        return(final_res)
    }

    final_res = res[-4]
    final_res$cluster_info = merged_data$cluster_info # Update cluster info based on merged ones.
    final_res[[4]]= merged_data$clustered_data
    colnames(final_res[[4]]) = round(as.numeric(colnames(final_res[[4]])), 4)
    names(final_res) = c("all_clusters", "JBA_seeds", "JBA_bins_expanded",
                         "JBA_data")
    return(final_res)
}

###########################################################################
## JBA_plotBins ##
JBA_plotBins = function(NMR_JBA, NMR_data, ct = 0.85, ref_sample = 1, xlim = NULL,
                        ylim = NULL) {

    ## Define parameters
    color = "lightcoral" # color for the correlation plot
    color_spec = "gray40"
    color_start_bin = "midnightblue"
    color_end_bin = "cornflowerblue"
    size_line = 0.7
    size_point = 1
    size_text = 11

    ## Get information regarding bin edges
    info = NMR_JBA$JBA_bins_expanded

    ## Get bin limits
    cluster_lim = do.call(rbind, strsplit(info[, "ppm_limits"], "_"))

    type = cbind(rep(1, nrow(cluster_lim)), rep(1, nrow(cluster_lim)),
                 rep(3, nrow(cluster_lim)))

    clusters_DF = data.frame(limits = c(as.numeric(cluster_lim[,1]),
                                        as.numeric(info[, "mean_ppm"]),
                                        as.numeric(cluster_lim[ ,2])),
                             value = c(type[,1], type[, 3], type[, 2]))

    ## Full-resolution spectrum
    ppm = as.numeric(colnames(NMR_data))
    fr_DF = as.data.frame(cbind(intensity = NMR_data[ref_sample, ],
                                ppm = ppm))

    ## Correlation_based spectrum
    clean_all_mean = NMR_JBA$all_clusters[, "mean_cor"]
    clean_all_mean[clean_all_mean < ct] = 0 # will be used for the color
    clean_all_mean [clean_all_mean != 0] = 1 # will be used for the color
    names(clean_all_mean) = NMR_JBA$all_clusters[, "mean_ppm"]
    max_clust = NMR_JBA$JBA_seeds[, "mean_ppm"]
    clean_all_mean[max_clust] = 2

    cor_dataframe = data.frame(ppm = as.numeric(NMR_JBA$all_clusters[, "mean_ppm"]),
                               mean_cor = as.numeric(NMR_JBA$all_clusters[, "mean_cor"]),
                               clean_mean = as.numeric(clean_all_mean))

    if (is.null(xlim)) {
        xlim = as.numeric(c(min(cluster_lim[, 2]), max(cluster_lim[, 1])))
    }
    xlim = sort(xlim, decreasing = TRUE)

    if (is.null(ylim)) {
        idx = which(ppm >= xlim[2] & ppm <= xlim[1])
        ylim = c(0, max(NMR_data[ref_sample, idx]))
    }

    plot_clust = ggplot(clusters_DF, aes(limits, value, color = type)) +
        geom_line(aes(limits, value), color = color, size = size_line) +
        geom_point(aes(limits, value), color = color) +
        scale_x_reverse() +
        xlim (xlim) +
        xlab("ppm") +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

    plot_full_res = ggplot(fr_DF) +
        geom_line(aes(ppm, intensity), color  = color_spec, size = size_line) +
        scale_x_reverse(limits = xlim) +
        xlab ("ppm") +
        ylab("Intensity") +
        scale_y_continuous(limits = ylim) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = size_text),
            axis.title = element_text(size = size_text, vjust = 0)) +
        geom_vline(xintercept = as.numeric(cluster_lim[, 1]), linetype="dashed",
                   color = color_start_bin, size = 0.6) +
        geom_vline(xintercept = as.numeric(cluster_lim[, 2]), linetype="dashed",
                   color = color_end_bin, size = 0.6)

    plot_cor = ggplot(cor_dataframe, aes(ppm, mean_cor, color = clean_mean)) +
        geom_line(aes(ppm, mean_cor), color  =  "gray", size = size_line) +
        geom_point(aes(ppm, mean_cor), size = size_point) +
        scale_x_reverse(limits=xlim) +
        scale_y_continuous(limits = c(0.6, 1)) +
        ylab("Average correlation") +
        xlab("ppm") +
        geom_hline(yintercept = ct, colour = "red", linetype = "dashed") +
        scale_color_gradient(low = "gray", high = color,guide = FALSE) +
        theme_bw() +
        theme(#axis.title.x=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = size_text),
            axis.title = element_text(size = size_text, vjust = 0))

    plot_fr = ggplot_gtable(ggplot_build(plot_full_res))
    plot_clust = ggplot_gtable(ggplot_build(plot_clust))
    plot_cor = ggplot_gtable(ggplot_build(plot_cor))

    maxWidth = unit.pmax(plot_fr$widths[2:3], plot_clust$widths[2:3],
                         plot_clust$cor[2:3])

    plot_fr$widths[2:3] <- maxWidth
    plot_clust$widths[2:3] <- maxWidth
    plot_cor$widths[2:3] <- maxWidth

    #grid.arrange(plot_fr, plot_cor, plot_clust, heights = c(1, 1, 1))
    #g <- arrangeGrob(plot_fr, plot_cor, plot_clust, heights = c(1, 1, 1))

    grid.arrange(plot_fr, plot_cor, heights = c(1, 1))
    g <- arrangeGrob(plot_fr, plot_cor, heights = c(1, 1))

    return(g)
}
