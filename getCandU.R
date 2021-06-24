library(phangorn)
library(foreach)
library(doParallel)

getCandU <- function(tree, t = NULL,
                    score_data, node_column,
                    p_column, sign_column,
                    threshold = 0.05,
                    cores = 1,
                    message = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }
    
    # t values
    if (is.null(t)) {
        t <- c(0,
               seq(0.01, 0.04, by = 0.01), 
               seq(0.05, 1, by = 0.05))
    }
    
    appCores=1
    if(cores != 1)
        appCores <- 4
    
    # a list to store levels under different t 
    # columns: p value, sign, node
    p_col <- score_data[[p_column]]
    sign_col <- score_data[[sign_column]]
    node_col <- score_data[[node_column]]
    
    # paths
    path <- matTree(tree = tree)
    keep <- !is.na(p_col)
    node_keep <- score_data[[node_column]][keep]
    node_all <- printNode(tree = tree, type = "all")
    node_in <- node_all[["nodeNum"]][!node_all[["isLeaf"]]]
    node_in <- intersect(node_keep, node_in)
    
    # nodes with p value available
    leaf <- .pseudoLeaf(tree = tree, score_data = score_data, 
                        node_column = node_column,
                        p_column = p_column, cores = appCores)
    
    for(i in seq_along(t))
    {
        name_q <- paste0("q_", t[i])
        score_data[[name_q]] <- ifelse(p_col > t[i], 0,
                                       1) * sign(sign_col)
    }
    
    chl_I <- Descendants(tree, node_in, type = "children")
    
    # for each t
    registerDoParallel(cores)
    level_list <- foreach(i=1:length(seq_along(t))) %dopar% {
        if (message) {
            message("Searching candidates on t =  ", t[i], " ...")
        }
        
        # add a column to store the result of search: keep
        if (any(colnames(score_data) == "keep")) {
            stop("The result will be output in the 'keep' column;
             Please use other name for the current 'keep' column.")
        }
        
        # For an internal node, if more than half of its direct child nodes has
        # NA score, it would not be picked.
        name_q <- paste0("q_", t[i])
        #chl_I <- share(chl_I)
        sel_0 <- mclapply(chl_I, FUN = function(x){
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            sum(!is.na(qx))/length(qx) > 0.5
        }, mc.cores = appCores)
        gc()
        sel_0 <- unlist(sel_0)
        node_0 <- node_in[sel_0]
        
        # For an internal nodes, if itself and all its descendant have q score
        # equals 1 or -1, pick the node
        br_I <- findDescendant(tree = tree, node = node_0, only.leaf = FALSE,
                               self.include = TRUE, use.alias = TRUE)
        #br_I <- share(br_I)
        sel_1 <- mclapply(br_I, FUN = function(x){
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            abs(mean(qx, na.rm = TRUE)) == 1
        }, mc.cores = appCores)
        gc()
        sel_1 <- unlist(sel_1)
        
        # filter by threshold
        sel_2 <- p_col[match(node_0, node_col)] <= threshold
        
        node_1 <- node_0[sel_1 & sel_2]
        
        # remove nodes whose ancestor is selected
        ind0 <- apply(path, 2, FUN = function(x) {
            x %in% node_1
        })
        ind0 <- which(ind0, arr.ind = TRUE)
        rs <- split(seq_len(nrow(ind0)), ind0[, "row"])
        rl <- unlist(lapply(rs, length)) > 1
        rs <- rs[rl]
        node_rm <- lapply(rs, FUN = function(x) {
            xx <- ind0[x, , drop = FALSE]
            mx <- xx[, "col"] < max(xx[, "col"])
            fx <- xx[mx, , drop = FALSE]
            path[fx]
        })
        
        node_rm <- unlist(node_rm)
        node_2 <- setdiff(node_1, node_rm)
        desd_2 <- findDescendant(tree = tree, node = node_2, 
                                 only.leaf = FALSE, self.include = TRUE)
        desd_2 <- unlist(desd_2)
        
        c(setdiff(leaf, desd_2), node_2)
        
    }
    gc()
    names(level_list) <- c(t)
    # q
    
    out <- list(candidate_list = level_list,
                score_data = score_data)    
    return(out)
    
}

.pseudoLeaf <- function(tree, score_data, node_column, p_column, cores = 1) {
    mat <- matTree(tree = tree)
    nd <- score_data[[node_column]][!is.na(score_data[[p_column]])]
    exist_mat <- apply(mat, 2, FUN = function(x) {x %in% nd})
    
    ww <- which(exist_mat, arr.ind = TRUE)
    ww <- ww[order(ww[, 1]), , drop = FALSE]
    loc_leaf <- ww[!duplicated(ww[, 1]), ]
    leaf_0 <- unique(mat[loc_leaf])
    #leaf_0 <- share(leaf_0)
    ind_0 <- mclapply(leaf_0, FUN = function(x) {
        xx <- which(mat == x, arr.ind = TRUE)
        ux <- xx[!duplicated(xx), , drop = FALSE]
        y0 <- nrow(ux) == 1
        if (nrow(ux) > 1) {
            ux[, "col"] <- ux[, "col"] - 1
            y1 <- all(!exist_mat[ux])
            y0 <- y0 | y1
        }
        return(y0)
    }, mc.cores = cores)
    gc()
    leaf_1 <- leaf_0[unlist(ind_0)]
    
    
    return(leaf_1)
}
