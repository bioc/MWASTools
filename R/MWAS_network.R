### build_cornet_edge ####
build_cornet_edges = function(metabolite, index, CorMat) {

    index1 = index + 1
    rangeL = length(index1:ncol(CorMat))

    N1 = rep(metabolite, rangeL)
    N2 = colnames(CorMat)[index1:ncol(CorMat)]
    values = CorMat[index, index1:ncol(CorMat)]

    types = values
    types[types > 0] = "pos_assoc"
    types[types < 0] = "neg_assoc"

    edges = cbind(N1, N2, values, types)
    return(edges)
}

### find_score ####
find_score = function(metabolite, MWASscores) {
    ind_pvalue = which(names(MWASscores) == metabolite)[1]  #to be sure
    pvalue_met = MWASscores[ind_pvalue]
    res = c(metabolite, pvalue_met)
    return(res)
}

### match_igAttribute ###
match_igAttribute = function(metabolite, attributes_met, alpha_th) {
    score_th = -log10(alpha_th)
    ind_metabolite = which(attributes_met[, 1] == metabolite)
    att = as.numeric(attributes_met[ind_metabolite, 2])
    if (att < (-score_th)) {
        col = "cornflowerblue"
    } else if (att > score_th) {
        col = "firebrick1"
    } else {
        col = "gray"
    }
    res = c(att, col)
    return(res)
}

### remove_isolatedV ###
remove_isolatedV = function(metabolite, all_metabolites) {
    if (metabolite %in% all_metabolites == FALSE) {
        res = "isolated"
    } else {
        res = "connected"
    }
}

##### MWAS_network ###
MWAS_network = function(metabo_SE, MWAS_matrix, alpha_th = 0.05,
    cor_th = 0.25, file_name = "Correlation", res_cor = 2) {

    ## Check that input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    metabo_matrix = t(assays(metabo_SE)$metabolic_data)

    if (!is.matrix(MWAS_matrix)) {
        stop("MWAS_matrix must be a numeric matrix")
    }
    if (ncol(MWAS_matrix) < 3) {
        stop("MWAS_matrix seems not to have the correct format")
    }
    if (nrow(MWAS_matrix) != ncol(metabo_matrix)) {
        stop("metabo_SE and MWAS_matrix are not consistent")
    }
    metabo_ids = colnames(metabo_matrix)
    num_answer = suppressWarnings(!is.na(as.numeric(metabo_ids[1])))

    if (num_answer == TRUE) {
        metabo_ids = paste("m", metabo_ids, sep = "")
    }
    if (length(metabo_ids) != length(unique(metabo_ids))) {
        stop("metabo_ids must be unique")
    }

    ## Select submatrix of metabolites based on p-value th
    MWASpvalues = MWAS_matrix[, 3]
    MWASestimates = MWAS_matrix[, 1]
    MWASestimates[MWASestimates < 0] = -1
    MWASestimates[MWASestimates >= 0] = 1
    MWASscores = -log10(MWASpvalues) * MWASestimates

    colnames(metabo_matrix) = metabo_ids
    names(MWASscores) = metabo_ids

    ## Build correlation matrix##
    CorMat = cor(metabo_matrix, metabo_matrix, use = "complete.obs")
    CorMat = round(CorMat, res_cor)

    ## Build igraph network
    CorMat_W = CorMat
    CorMat_W[abs(CorMat_W) <= cor_th] = 0
    CorMat_W[row(CorMat_W) == col(CorMat_W)] = 0
    igNet = graph.adjacency(CorMat_W, mode = "undirected", weighted = TRUE)

    ## Build cytoscape network
    metabolites = colnames(CorMat)[-ncol(CorMat)]
    index_vector = seq_along(metabolites)
    NL = list()
    NL[["CorMat"]] = CorMat
    all_edges = mapply(build_cornet_edges, metabolites, index_vector,
        MoreArgs = NL, SIMPLIFY = FALSE)
    cor_network = do.call(rbind, all_edges)

    ## Filter cor_network
    index_net = which(abs(as.numeric(cor_network[, 3])) > cor_th)

    if (length(index_net) == 0) {
        stop("Impossible to build a network with the current cor_th")
    } else {
        cor_network = cor_network[index_net, ]
        cor_network = matrix(cor_network, ncol = 4)  # force it to be a matrix
        colnames(cor_network) = c("node1", "node2", "r.coeff",
            "type")
        rownames(cor_network) = NULL
    }

    ## Create attributes cytoscape
    all_metabolites = unique(as.vector(cor_network[, 1:2]))
    attributes_met = lapply(all_metabolites, find_score, MWASscores)
    attributes_met = do.call(rbind, attributes_met)
    colnames(attributes_met) = c("node", "MWASpvalue")

    all_metabolitesi = rownames(as.matrix(unlist(V(igNet))))
    V_type = sapply(all_metabolitesi, remove_isolatedV, all_metabolites)
    isolated_V = names(which(V_type == "isolated"))

    ## Remove isolated nodes
    if (length(isolated_V) > 0) {
        igNet = delete_vertices(igNet, isolated_V)
    }

    all_metabolitesi = rownames(as.matrix(unlist(V(igNet))))
    attributes_meti = lapply(all_metabolitesi, match_igAttribute,
        attributes_met, alpha_th = alpha_th)
    attributes_meti = do.call(rbind, attributes_meti)
    attributes_col = as.character(attributes_meti[, 2])
    attributes_score = as.numeric(attributes_meti[, 1])
    V(igNet)$color = attributes_col
    V(igNet)$score = attributes_score

    file_nameN = paste(file_name, "AttributeFile.txt", sep = "_")
    write.table(attributes_met, file_nameN, row.names = FALSE,
        sep = "\t", quote = FALSE, col.names = TRUE)

    ## Export network to cytoscape
    cytoscape_net = as.data.frame(cor_network, rownames = NULL)
    file_nameN = paste(file_name, "NetworkFile.txt", sep = "_")
    write.table(cytoscape_net, file_nameN, row.names = FALSE,
        sep = "\t", quote = FALSE, col.names = TRUE)

    return(igNet)
}

