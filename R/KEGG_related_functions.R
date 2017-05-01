globalVariables("KEGG_metabolic_paths")

## INTERNAL FUNCTIONS ##

## convertTable ## From MetaboSignal
convertTable = function(res) {
    if (nchar(res) == 0) {
        result = NULL
    } else {
        rows = strsplit(res, "\n")
        rows.len = length(rows[[1]])
        result = matrix(unlist(lapply(rows, strsplit, "\t")), nrow = rows.len,
                        byrow = TRUE)
  }
    return(result)
}

## link_cpd_name ## Transform one metabolite ID into common name
link_cpd_name = function(compound, compoundM) {
    name = compound # by default no change
    index = which(compoundM[, 1] == compound)
    if (length(index) > 0) {
        names = compoundM[index, 2]
        names = trimws(unlist(strsplit(names, "[;]")))
        #names = unlist(strsplit(names, "[,]"))
        if (nchar(names[1]) > 15) { ## name too long
            shortest_ind = which(nchar(names) == min(nchar(names)))[1]
            name = names[shortest_ind][1]
        } else {
            name = names[1]
        }
        name = gsub(" ", "-", name)
    }
    return(name)
}

## MWAS_ChangeNames ## Transform all metabolites ids into common names
MWAS_ChangeNames = function(metabolites) {
    ## Get compounds table
    file = "http://rest.kegg.jp/list/compound"
    response = getURL(file)
    compoundM = convertTable(response)
    all_names = sapply(metabolites, link_cpd_name, compoundM)
    return(all_names)
}

### find_index ####
find_index = function(name, all_names) {
    return(which(all_names == name))
}

### get_metabonetR #### Parse KEGG maps
get_metabonetR = function(path) {
    message(path)
    if (substr(path, 4, nchar(path)) == "01100") {
        # remove metabolic pathways map
        parsed_path = NULL
        message
    } else {
        # Check that the input path exists
        file = paste("http://rest.kegg.jp/get/", path, "/kgml",
            sep = "")
        pathway = try(getURL(file), silent = TRUE)
        reactions = try(getReactions(parseKGML(pathway)), silent = TRUE)

        if (grepl("Error", reactions[1]) == TRUE) {
            to_print = paste(path, "-path ID without XML:path removed",
                sep = "")
            message(to_print)
            parsed_path = NULL
        } else {
            parsed_path = capture.output(reactions, file = NULL)
        }
    }
    return(parsed_path)
}

### new_sr_edge #### create individual substrate_reaction edges
new_sr_edge = function(reaction, substrates, Type) {
    substrate_lines = cbind(substrates, rep(reaction, length(substrates)),
        rep(Type, length(substrates)))
    return(substrate_lines)
}

### add_sr_edge #### create all substrate_reaction edges for a given reaction
add_sr_edge = function(Reaction, Substrate, Type) {

    reactions = unlist(strsplit(Reaction, " "))  # Some reactions have several IDs
    substrates = unlist(strsplit(Substrate, ";"))
    substrates = unlist(strsplit(substrates, " "))

    # Correct for direction based on Duarte et al., 2007## this was updated
    for (reaction in reactions) {
        directionality_reactions = KEGG_metabolic_paths[[2]]
        a = which(directionality_reactions[, 1] == reaction)
        if (length(a) >= 1) {
            Type = directionality_reactions[, 2][a]
        }
    }

    # Build all substrate_reaction edges of a given reaction
    all_edges = lapply(reactions, new_sr_edge, substrates, Type)
    all_edges = unique(do.call(rbind, all_edges))
    return(all_edges)
}

### sr_network_from_path #### create substrate_reaction edges
### for all reactions
sr_network_from_path = function(path, list_parsed_paths) {
    lines = list_parsed_paths[[path]]
    # print(path) Select reactions
    global_lines_reactions = intersect(grep("rn", lines), grep("Name",
        lines))
    Reactions = lines[c(global_lines_reactions)]
    Reactions = substr(Reactions, 11, nchar(Reactions))

    ## Select substrates
    global_lines_substrates = grep("Substrate Name", lines)
    Substrates = lines[c(global_lines_substrates)]
    Substrates = substr(Substrates, 21, nchar(Substrates))

    # Select type
    global_lines_type = grep("Type", lines)
    Type = lines[c(global_lines_type)]
    Type = substr(Type, 11, nchar(Type))

    ## Build network edges
    all_edges = mapply(add_sr_edge, Reactions, Substrates, Type, SIMPLIFY = FALSE)
    all_edges = unique(do.call(rbind, all_edges))
    return(all_edges)
}

### new_sp_edge #### create individual substrate_product_reaction edges
new_sp_edge = function(row_substrate, mid_value) {
    if (identical(as.character(mid_value),
                  as.character(unlist(row_substrate)[1:2]))) {
        message(" -Process: 50% completed")
    }
    reaction = unlist(row_substrate)[2]
    substrate = unlist(row_substrate)[1]
    Type = unlist(row_substrate)[3]

    rp = keggGet(reaction)[[1]]$RCLASS
    rp = gsub("RC", "rc", rp)  # to avoid confusion with C of compound
    rp = unlist(strsplit(rp, "  "))
    #print(reaction)
    #print(rp)
    substrate_rf = gsub("cpd:", "", substrate)
    substrate_rf = gsub("gl:", "", substrate_rf)
    ind_rp = grep(substrate_rf, rp)
    edges_reaction = c()

    if (length(ind_rp) == 0) {
        return(edges_reaction)
    }

    for (z in 1:length(ind_rp)) {
        rp_mapped = unlist(strsplit(rp[ind_rp[z]], "_"))
        rp_mapped = unlist(strsplit(rp_mapped, " "))
        # print(rp_mapped) print(substrate_rf)
        ind_s = grep(substrate_rf, rp_mapped)
        # print(ind_s)
        if (ind_s == 1) {
            new_edge = c(substrate_rf, rp_mapped[2], reaction, Type)
        } else {
            new_edge = c(substrate_rf, rp_mapped[1], reaction, Type)
        }
        edges_reaction = rbind(edges_reaction, new_edge)
    }

    edges_reaction[, 1:2] = paste("cpd:", edges_reaction[, 1:2],
        sep = "")
    if (Type == "reversible") {
        reverse = cbind(edges_reaction[, 2], edges_reaction[,
            1], edges_reaction[, 3], edges_reaction[, 4])
        edges_reaction = rbind(edges_reaction, reverse)
    }
    rownames(edges_reaction) = NULL
    colnames(edges_reaction) = c("node1", "node2", "reaction", "type")
    return(edges_reaction)
}

### MWAS_build_reaction_network ####
MWAS_build_reaction_network = function(metabo_paths) {
    ### Get KGML files and transform them into reaction files####
    message("Reading paths from KGML files")
    list_parsed_paths = lapply(metabo_paths, get_metabonetR)
    names(list_parsed_paths) = metabo_paths
    path_names = metabo_paths

    if (length(list_parsed_paths) == 0) {
        to_print = ("Impossible to build a metabolic network")
        stop(to_print, "\n")

    } else {
        ## Remove empty paths: for example oxidative phosphorylation#
        length_paths = sapply(list_parsed_paths, length)
        empty_paths = which(length_paths <= 1)

        if (length(empty_paths) >= 1) {
            list_parsed_paths = list_parsed_paths[-c(empty_paths)]
            paths_included[path_names[empty_paths]] = 0
            path_names = path_names[-c(empty_paths)]
        }
        if (length(list_parsed_paths) == 0) {
            to_print = ("Impossible to build a metabolic network")
            stop(to_print, "\n")
        } else {
            ### Create metabolic_table #### 1) Create substrate_reaction
            ### table
            substrate_table = lapply(path_names, sr_network_from_path,
                list_parsed_paths)

            substrate_table = unique(do.call(rbind, substrate_table))

            colnames(substrate_table) = c("substrate", "reaction", "type")

            # 2) Create substrate_product_reaction table based on
            # reactant_pairs
            rows_substrate = split(substrate_table, row(substrate_table))
            message()
            message("Building reaction network: might take few minutes")
            ptm <- proc.time()
            mid_value = rows_substrate[[round(length(rows_substrate)/2,
                0)]][1:2]
            metabolic_table = lapply(rows_substrate, new_sp_edge, mid_value)
            metabolic_table = do.call(rbind, metabolic_table)
            proc.time() - ptm
            if (is.null(metabolic_table)) {
                stop("Impossible to build reaction network")
            }
            metabolic_table = matrix(metabolic_table, ncol = 4)
            colnames(metabolic_table) = c("node1", "node2", "reaction", "type")
            return(metabolic_table)
        }
    }
}

### ASP_paths ####
ASP_paths = function(ASP) {
    shortpath = rownames(as.matrix(unlist(ASP)))
    return(shortpath)
}

### path_as_network #### reformat shortest paths
path_as_network = function(path) {
    all_edges = c()
    for (i in 1:(length(path) - 1)) {
        edge = c(path[i], path[i + 1])
        all_edges = rbind(all_edges, edge)
    }
    return(all_edges)
}

### target_metabolite_SP #### get shortest path(s) from a
### source to a target
target_metabolite_SP = function(target_metabolite, source_metabolite,
    reaction_network_i, network_table, distance_table, distance_th,
    type) {

    index_source = which(rownames(distance_table) == source_metabolite)
    index_target = which(colnames(distance_table) == target_metabolite)
    distanceGM = distance_table[index_source, index_target]

    if (distanceGM < distance_th) {
        # If distance is not Inf,

        if (type == "first") {
            ASP = get.shortest.paths(reaction_network_i, source_metabolite,
                target_metabolite, mode = "out")
            ASP = ASP[[1]]
        } else {
            ASP = get.all.shortest.paths(reaction_network_i,
                source_metabolite, target_metabolite, mode = "out")
            ASP = ASP$res
        }

        all_paths = unique(do.call(rbind, lapply(ASP, ASP_paths)))
        rownames(all_paths) = NULL

        ## Transform the shortest path matrix into a network table
        all_pathsList = split(all_paths, row(all_paths))
        all_pathsnetworkList = lapply(all_pathsList, path_as_network)
        subnetwork = unique(do.call(rbind, all_pathsnetworkList))
        rownames(subnetwork) = NULL  # subnetwork is a network-table

        ## Add edges for reversible interactions
        reverseSubnetwork = cbind(subnetwork[, 2], subnetwork[, 1])
        supernetwork = rbind(reverseSubnetwork, network_table)
        intersection_reverseSubnetwork = supernetwork[duplicated(supernetwork), ]
        intersection_reverseSubnetwork = matrix(intersection_reverseSubnetwork,
                                                ncol = 2)

        if (nrow(intersection_reverseSubnetwork) > 0) {
            subnetwork = rbind(subnetwork, intersection_reverseSubnetwork)
        }
        pathsGM = unique(subnetwork)
        rownames(pathsGM) = NULL

    } else {
        pathsGM = NULL
    }
    return(pathsGM)
}

### source_metabolite_SP #### shortest paths from a source to
### several targets
source_metabolite_SP = function(source_metabolite, target_metabolites,
    reaction_network_i, network_table, distance_table, distance_th,
    type) {

    ind_source = which(target_metabolites == source_metabolite)
    target_metabolites_new = target_metabolites[-ind_source]

    metabo_pathsGM = lapply(target_metabolites_new, target_metabolite_SP,
        source_metabolite, reaction_network_i, network_table,
        distance_table, distance_th, type)

    metabo_pathsGM = unique(do.call(rbind, metabo_pathsGM))
    return(metabo_pathsGM)
}

### MWAS_subnetwork_SP #### create shortest path subnetwork between metabolites
MWAS_subnetwork_SP = function(network_table, metabolites, type = "all",
                              distance_th = Inf) {

    ## Get BU ##
    net_reaction = unique(network_table[, 1:3])

    ## Get distance table#
    network_table = unique(network_table[, 1:2])
    reaction_network_i = graph.data.frame(network_table, directed = TRUE)
    distance_table = distances(reaction_network_i, mode = "out")

    ## Check if the metabolites are included in the network_table
    common_metabolites = intersect(metabolites, unique(as.vector(network_table)))
    source_metabolites = common_metabolites
    target_metabolites = common_metabolites

    ## Calculate shortest path network
    message("Building shortest path network")
    all_pathsGM = lapply(source_metabolites, source_metabolite_SP,
        target_metabolites, reaction_network_i, network_table,
        distance_table, distance_th, type)
    all_pathsGM = unique(do.call(rbind, all_pathsGM))

    if (length(all_pathsGM) >= 1) {
        all_pathsGM = matrix(all_pathsGM, ncol = 2)
        colnames(all_pathsGM) = c("node1", "node2")
        return(all_pathsGM)

    } else {
        to_print = paste("Impossible to calculate shortest path: the metabolites
                         are not directly connected", sep = " ")
        stop(to_print)
    }
}

### reformat_network #### collapse reactions that are involved in the same reaction
### and have the same directionality
reformat_network = function(collapsed_edge, network_table_to_rf) {
    ind_edge = which(network_table_to_rf[, 1] == collapsed_edge)
    reactions_edge = paste(sort(network_table_to_rf[ind_edge, 2]), collapse = ";")
    uncollapsed_edge = unlist(strsplit(collapsed_edge, "_"))
    collapsed_cpd = paste(uncollapsed_edge[1], uncollapsed_edge[2], sep = "_")
    new_edge = c(collapsed_cpd, reactions_edge, uncollapsed_edge[3])
    return(new_edge)
}

### link_sp_reaction #### link edges of the shortest path subnetwork to reactions
link_sp_reaction = function(sp_edge, network_table_rf) {
    ind_edge = which(network_table_rf[, 1] == sp_edge)
    new_edges = cbind(rep(sp_edge, length(ind_edge)), network_table_rf[ind_edge, 2],
                      network_table_rf[ind_edge, 3])
    return(new_edges)

}

### collapse_directionality #### remove one of the reversible edges
collapse_directionality = function(sp_edge_rf) {

  uncollapsed_edge = unlist(strsplit(sp_edge_rf, "_"))

  if(uncollapsed_edge[4] == "irreversible") {
      return(uncollapsed_edge)
  } else {
      sorted_cpd = sort(uncollapsed_edge[1:2])
      new_edge = c(sorted_cpd, uncollapsed_edge[3], uncollapsed_edge[4])
      return(new_edge)
  }
}

### pathway_cpd_table #### link compounds to their pathways
pathway_cpd_table = function(compound, linked_table) {
    ind = which(linked_table[, 1] == compound)
    path_names = paste(linked_table[ind, 2], collapse = ";")
    res = c(compound, path_names)
    return(res)
}

### path_enrichment_test #### test for enrichment
path_enrichment_test = function(path, mapped_pathsM, cpd_pathM) {
    my_success = length(grep(path, mapped_pathsM[, "pathway_KEGG_ID"]))
    my_failure = length(unique(mapped_pathsM[, 1])) - my_success
    external_success = length(grep(path, cpd_pathM[, 2]))
    external_failure = nrow(cpd_pathM) - external_success
    #external_freq = length(grep(path, cpd_pathM[, 2]))/nrow(cpd_pathM)
    contingence_table = matrix(c(my_success, my_failure,
                                 external_success, external_failure),
                               ncol = 2)
    fisher_pval = fisher.test(contingence_table, alternative = "greater")$p.val
    path_name = mapped_pathsM[(mapped_pathsM[, 3] == path), 4][1]
    ans = c(path, path_name, fisher_pval)
    return(ans)
}


## EXTERNAL FUNCTIONS ##

### MWAS_KEGG_network ####
MWAS_KEGG_network = function(kegg_paths = NULL) {

    #### Check that impute data are correct
    if (is.null(kegg_paths)) {
        kegg_paths = KEGG_metabolic_paths[[1]][, 1] ## this was updated
    }

    if (!is.vector(kegg_paths)) {
        stop("paths must be a character vector")
    }

    ## Build reaction network
    network_table = MWAS_build_reaction_network(kegg_paths)
}

### MWAS_KEGG_shortestpaths ####
MWAS_KEGG_shortestpaths = function(network_table, metabolites,
    MWAS_matrix = NULL, type = "all", distance_th = "Inf", names = TRUE,
    file_name = "KeggSP") {

    ## Check if input data are correct
    if (!is.vector(metabolites)) {
        stop("metabolites must be a character vector")
    }

    if (!is.matrix(network_table)) {
        stop("network_table must be a matrix")
        if (ncol(network_table) != 4) {
            stop("network_table must have 4 columns")
        }
    }

    metabolites = unique(metabolites)
    common_metabolites = intersect(metabolites, unique(as.vector(network_table[,
        1:2])))
    if (length(common_metabolites) < 2) {
        to_print = paste("Impossible to calculate shortest paths: less than two",
            "metabolites were mapped onto the KEGG reaction network",
            sep = " ")
        stop(to_print)
    }

    if (!is.null(MWAS_matrix)) {
        if (!is.matrix(MWAS_matrix)) {
            stop("MWAS_matrix must be a matrix. See MWAS_stats")
        }
        if (ncol(MWAS_matrix) < 3) {
            stop("MWAS_matrix must have at least 3 columns. See MWAS_stats")
        }
        if (nrow(MWAS_matrix) != length(metabolites)) {
            stop("nrow of MWAS_matrix must be equal to metabolites length")
        }
        ## Get scores
        estimates = MWAS_matrix[, 1]
        estimates[estimates < 0] = -1
        estimates[estimates > 0] = 1
        logpvalue = -log10(MWAS_matrix[, 3])
        scores_metabolites = estimates * logpvalue
        names(scores_metabolites) = metabolites
    }

    ## Build shortest path subnetwork
    subnetwork = try(MWAS_subnetwork_SP(network_table, metabolites,
                                        type = type, distance_th = distance_th),
                     silent = TRUE)

    if (grepl("Error", subnetwork[1])) {
        subnetwork = NULL
        stop("Impossible to build shortest path subnetwork")
    }

    ## Reformat network_table ##
    network_table_to_rf = cbind(paste(network_table[, 1], network_table[, 2],
                                      network_table[, 4], sep = "_"),
                                network_table[, 3])
    network_table_rf = lapply(network_table_to_rf[, 1], reformat_network,
                              network_table_to_rf)
    network_table_rf = unique(do.call(rbind, network_table_rf))

    ## Reformat subnetwork ##
    subnetwork_to_rf = paste(subnetwork[, 1], subnetwork[, 2], sep = "_")
    subnetwork_rf = lapply(subnetwork_to_rf, link_sp_reaction, network_table_rf)
    subnetwork_rf = unique(do.call(rbind, subnetwork_rf))

    edges_subnetwork_rf = split(subnetwork_rf, row(subnetwork_rf))

    subnetwork_final = lapply(edges_subnetwork_rf, collapse_directionality)
    subnetwork_final = unique(do.call(rbind, subnetwork_final))

    colnames(subnetwork_final) = c("node1", "node2", "reaction", "direction")

    ## Change names of subnetwork ##
    if(names == TRUE) {
      subnetwork_final[, 1:2] = MWAS_ChangeNames(subnetwork_final[, 1:2])
    }

    ## Generate cytoscape file ##
    cytoscape_sp_net = as.data.frame(subnetwork_final, rownames = NULL)
    file_nameSP = paste(file_name, "NetworkFile.txt",
                       sep = "_")
    write.table(cytoscape_sp_net, file_nameSP, row.names = FALSE,
                sep = "\t", quote = FALSE, col.names = TRUE)

    ## Get scores
    if (!is.null(MWAS_matrix)) {
        scores_M = cbind(names(scores_metabolites), scores_metabolites)
        if (names == TRUE) {
          scores_M[, 1] = MWAS_ChangeNames(scores_M[, 1])
        }
        colnames(scores_M) = c("metabolite", "MWAS_score")
        cytoscape_score_net = as.data.frame(scores_M, rownames = NULL)
        file_nameS = paste(file_name, "AttributeFile.txt",
            sep = "_")
        write.table(cytoscape_score_net, file_nameS, row.names = FALSE,
            sep = "\t", quote = FALSE, col.names = TRUE)
    }

    if (length(setdiff(metabolites, as.vector(network_table))) > 0) {
      message ("Note: one or more metabolites were not mapped onto the network")
    }
    return(subnetwork_final)
}

### MWAS_KEGG_pathways ####
MWAS_KEGG_pathways = function(metabolites, MWAS_matrix = NULL,
                              enrichment_test = FALSE, mt_method = "bonferroni",
                              file_name = "KeggPaths") {

    ## Check if input data are correct
    if (!is.vector(metabolites)) {
        stop("metabolites must be a character vector")
    }
    metabolites = unique(metabolites)

    if (!is.null(MWAS_matrix)) {
        if (!is.matrix(MWAS_matrix)) {
            stop("MWAS_matrix must be a matrix. See MWAS_stats")
        }
        if (ncol(MWAS_matrix) < 3) {
            stop("MWAS_matrix must have at least 3 columns. See MWAS_stats")
        }
        if (nrow(MWAS_matrix) != length(metabolites)) {
            stop("nrow of MWAS_matrix must be equal to metabolites length")
        }
        ## Get scores
        estimates = MWAS_matrix[, 1]
        estimates[estimates < 0] = -1
        estimates[estimates > 0] = 1
        logpvalue = -log10(MWAS_matrix[, 3])
        scores_metabolites = estimates * logpvalue
        names(scores_metabolites) = metabolites
    }

    ## Build pathway subnetwork ##
    message("Building pathway subnetwork")
    mapped_paths = try(keggLink("pathway", metabolites), silent = TRUE)

    if (grepl("Error", mapped_paths[1])) {
        stop("None of the metabolites could be mapped onto the KEGG pathways")
    }

    mapped_paths = gsub("path:", "", mapped_paths)
    mapped_paths = cbind(names(mapped_paths), as.character(mapped_paths))
    unique_paths = unique(mapped_paths[, 2])
    unique_paths_names = vector(length = length(unique_paths))
    unique_paths_class = vector(length = length(unique_paths))
    unique_human_class = rep("Not_human", length(unique_paths))

    all_human_paths = unique(names(keggLink("hsa", "pathway")))
    all_human_maps = gsub("path:hsa", "map", all_human_paths)

    for (i in 1:length(unique_paths)) {
        path_info = keggGet(unique_paths[i])[[1]]
        unique_paths_names[i] = path_info$NAME
        if (is.null(path_info$CLASS)) {
            unique_paths_class[i] = "Not_determined"
        } else {
            unique_paths_class[i] = path_info$CLASS
        }
        if (unique_paths[i] %in% all_human_maps) {
            unique_human_class[i] = "Human"
        }
    }

    path_ind = sapply(mapped_paths[, 2], find_index, unique_paths)
    compound_names = MWAS_ChangeNames(mapped_paths[, 1])
    mapped_pathsM = cbind(mapped_paths[, 1], compound_names,
        mapped_paths[, 2], unique_paths_names[path_ind],
        unique_paths_class[path_ind], unique_human_class[path_ind])
    colnames(mapped_pathsM) = c("compound_KEGG_ID", "compound_name",
        "pathway_KEGG_ID", "pathway_name", "pathway_class", "pathway_organism")

    cytoscape_path_net = as.data.frame(mapped_pathsM, rownames = NULL)
    file_nameN = paste(file_name, "NetworkFile.txt", sep = "_")
    write.table(cytoscape_path_net, file_nameN, row.names = FALSE,
        sep = "\t", quote = FALSE, col.names = TRUE)

    ## Get scores
    if (!is.null(MWAS_matrix)) {
        scores_M = cbind(names(scores_metabolites), scores_metabolites)
        scores_M[, 1] = MWAS_ChangeNames(scores_M[, 1])
        colnames(scores_M) = c("metabolite", "MWAS_score")
        cytoscape_score_net = as.data.frame(scores_M, rownames = NULL)
        file_nameS = paste(file_name, "AttributeFile.txt",
            sep = "_")
        write.table(cytoscape_score_net, file_nameS, row.names = FALSE,
            sep = "\t", quote = FALSE, col.names = TRUE)
    }

    ## Enrichment test
    if(enrichment_test == FALSE) {
        return(mapped_pathsM)
    }

    ## Get compounds
    message()
    message("Running enrichment analysis")
    file = "http://rest.kegg.jp/link/pathway/cpd"
    response = getURL(file)
    linked_table = unique(convertTable(response))
    all_compounds = unique(linked_table[, 1])
    all_compounds = gsub("path:", "", all_compounds)

    ## Link compounds to paths
    cpd_pathM = unique(do.call(rbind, lapply(all_compounds,
                                             pathway_cpd_table, linked_table)))

    ## Do Fisher test
    paths_to_test = unique(mapped_pathsM[, "pathway_KEGG_ID"])
    paths_tested = do.call(rbind, lapply(paths_to_test, path_enrichment_test,
                                         mapped_pathsM, cpd_pathM))
    adjusted_pval = p.adjust(paths_tested[, 3], method = mt_method)
    paths_tested = cbind(paths_tested, adjusted_pval)
    colnames(paths_tested) = c("pathway_KEGG_ID", "pathway_name", "raw_pval",
                               "adjusted_pval")
    cytoscape_enrich_net = as.data.frame(paths_tested, rownames = NULL)
    file_nameE = paste(file_name, "EnrichmentAnalysis.txt", sep = "_")
    write.table(cytoscape_enrich_net, file_nameE, row.names = FALSE,
                sep = "\t", quote = FALSE, col.names = TRUE)

    return(mapped_pathsM)
}

