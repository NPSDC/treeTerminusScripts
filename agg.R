evalCand <- function(tree,
                     type = c("single", "multiple"),
                     levels = cand_list,
                     score_data = NULL,
                     node_column, p_column,
                     sign_column = sign_column,
                     feature_column = NULL,
                     method = "BH",
                     limit_rej = 0.05,
                     use_pseudo_leaf = FALSE,
                     message = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }
    
    type <- match.arg(type)
    if (type == "single") {
        score_data <- list(score_data)
        levels <- list(levels)
    }
    
    
    if (type == "multiple" & is.null(feature_column)) {
        warning("To distinct results from different features,
                feature_column is required")
    }
    
    node_list <- lapply(score_data, FUN = function(x) {
        x[[node_column]]
    })
    
    # ------------------------- the pseudo leaf level -------------------------
    # some nodes might not be included in the analysis step because they have no
    # enough data. In such case, an internal node would become a pseudo leaf
    # node if its descendant nodes are filtered due to lack of sufficient data.
    
    if (use_pseudo_leaf) {
        if (message) {
            message("collecting the pseudo leaf level for all features ...")
        }
        
        pseudo_leaf <- lapply(seq_along(score_data), FUN = function(x) {
            if (message) {
                message(x, " out of ", length(score_data),
                        " features finished", "\r", appendLF = FALSE)
                flush.console()}
            
            .pseudoLeaf(tree = tree, score_data = score_data[[x]],
                        node_column = node_column, p_column = p_column)
        })
        names(pseudo_leaf) <- names(score_data)
        
        
        if (message) {
            message("Calculating the number of pseudo-leaves of each node
                for all features ...")
        }
        info_nleaf <- lapply(seq_along(node_list), FUN = function(x) {
            if (message) {
                message(x, " out of ", length(node_list),
                        " features finished", "\r", appendLF = FALSE)
                flush.console()}
            
            xx <- node_list[[x]]
            ps.x <- pseudo_leaf[[x]]
            
            desd.x <- findDescendant(tree = tree, node = xx,
                                     only.leaf = FALSE, self.include = TRUE)
            leaf.x <- findDescendant(tree = tree, node = xx,
                                     only.leaf = TRUE, self.include = TRUE)
            psLeaf.x <- lapply(desd.x, FUN = function(x) {
                intersect(x, ps.x)})
            info <- cbind(n_leaf = unlist(lapply(leaf.x, length)),
                          n_pseudo_leaf = unlist(lapply(psLeaf.x, length)))
            
            return(info)
        })
        
        names(info_nleaf) <- names(score_data)
    } else {
        node_all <- showNode(tree = tree, only.leaf = FALSE)
        desc_all <- findDescendant(tree = tree, node = node_all,
                                   only.leaf = TRUE, self.include = TRUE)
        info_nleaf <- data.frame(
            node = node_all,
            n_leaf = unlist(lapply(desc_all, length)))
    }
    
    # add two columns in score_data
    # ---------------- info about the candidate level --------------------------
    # candidates in the candidate level
    # a data frame: t, br_size, candidate, method, limit_rej, level_name,
    # rej_leaf, rej_node, rej_pseudo_leaf
    
    if (message) {
        message("Evaluating candidates ... ")
    }
    tlist <- lapply(levels, names)
    t <- tlist[!duplicated(tlist)]
    if (length(t) > 1) {
        stop("the names of elements in 'levels' are different")
    }
    t <- unlist(t)
    t <- as.numeric(t)
    
    level_info <- data.frame(t = t, upper_t = NA,
                             is_valid = FALSE,
                             method = method,
                             limit_rej = limit_rej,
                             level_name = tlist[[1]],
                             best = FALSE,
                             rej_leaf = NA,
                             rej_node = NA,
                             rej_pseudo_leaf = NA,
                             rej_pseudo_node = NA)
    print("adss")
    sel <- vector("list", length(t))
    names(sel) <- tlist[[1]]
    for (i in seq_along(t)) {
        # message
        if (message) {
            message("working on ", i , " out of ",
                    length(t), " candidates \r", appendLF = FALSE)
            flush.console()
        }
        
        # get the candidate level at t[i]
        name_i <- as.character(level_info$level_name[i])
        level_i <- lapply(levels, FUN = function(x) {
            x[[i]]
        })
        sel_i <- mapply(function(x, y) {
            ii <- match(x, y[[node_column]])
        }, level_i, score_data, SIMPLIFY = FALSE)
        len_i <- lapply(sel_i, length)
        
        
        # adjust p-values
        p_i <- mapply(FUN = function(x, y) {
            x[y, p_column]
        }, score_data, sel_i, SIMPLIFY = FALSE)
        adp_i <- p.adjust(p = unlist(p_i), method = method)
        rej_i <- adp_i <= limit_rej
        
        # the largest p value that is rejected
        maxp_i <- max(c(-1, unlist(p_i)[rej_i]))
        
        # the number of branches
        path <- matTree(tree = tree)
        n_C <- mapply(FUN = function(x, y) {
            # nodes rejected in each feature
            xx <- x[y, c(node_column, sign_column, p_column)]
            xs <- xx[xx[[p_column]] <= maxp_i, ]
                
            # split nodes by sign
            sn <- split(xs[[node_column]], sign(xs[[sign_column]]))
            is_L <- lapply(sn, FUN = function(x) {
                isLeaf(tree = tree, node = x)})
            
            rej_L <- mapply(FUN = function(x, y) {
                unique(x[y])}, sn, is_L)
            rej_I <- mapply(FUN = function(x, y) {
                unique(x[!y]) }, sn, is_L)
            rej_L2 <- lapply(rej_L, FUN = function(x) {
                unique(path[path[, "L1"] %in% x, "L2"])})
            
            length(unlist(rej_I)) + length(unlist(rej_L2))
        }, score_data, sel_i, SIMPLIFY = FALSE)
        n_C <- sum(unlist(n_C))
        
        # The number of leaves
        if(use_pseudo_leaf) {
            rej_m1 <- mapply(FUN = function(x, y) {
                x[y, "n_pseudo_leaf"]
            }, info_nleaf, sel_i, SIMPLIFY = FALSE)
            n_m <- sum(unlist(rej_m1)[rej_i %in% TRUE])
            av_size <- n_m/max(n_C, 1)
            
        } else {
            node_i <- mapply(FUN = function(x, y) {
                x[y, node_column]
            }, score_data, sel_i, SIMPLIFY = FALSE)
            node_r <- unlist(node_i)[rej_i %in% TRUE]
            ind_r <- match(node_r, info_nleaf[["node"]])
            n_m <- sum(info_nleaf[ind_r, "n_leaf"])
            av_size <- n_m/max(n_C, 1)
        }
        
        # This is to avoid get TRUE from (2*0.05*(2.5-1)) > 0.15
        up_i <- min(2 * limit_rej * (max(av_size, 1) - 1), 1)
        up_i <- round(up_i, 10)
        
        
        
        #level_info$lower_t[i] <- low_i
        level_info$upper_t[i] <- up_i
        level_info$rej_leaf[i] <- n_m
        level_info$rej_node[i] <- sum(rej_i)
        
        if (use_pseudo_leaf) {
            level_info$rej_pseudo_leaf[i] <- n_m
            level_info$rej_pseudo_node[i] <- n_C
        }
        sel[[i]] <- sel_i
        
        level_info$is_valid[i] <- up_i > t[i] | t[i] == 0
    }
    print("adss")
    # candidates: levels that fullfil the requirement to control FDR on the
    # (pseudo) leaf level when multiple hypothesis correction is performed on it
    isB <- level_info %>%
        filter(is_valid) %>%
        filter(rej_leaf == max(rej_leaf)) %>%
        filter(rej_node == min(rej_node)) %>%
        select(level_name) %>%
        unlist() %>%
        as.character()
    level_info <- level_info %>%
        mutate(best = level_name %in% isB)
    level_b <- lapply(levels, FUN = function(x) {x[[isB[1]]]})
    
    # output the result on the best level
    if (message) {
        message("mulitple-hypothesis correction on the best candidate ...")
    }
    sel_b <- sel[[isB[1]]]
    #return(sel_b)
    outB <- lapply(seq_along(score_data), FUN = function(i) {
        si <- sel_b[[i]]
        score_data[[i]][si, , drop = FALSE]
    })
    
    outB <- rbindlist(outB)
    pv <- outB[[p_column]]
    apv <- p.adjust(pv, method = method)
    outB$adj.p <- apv
    outB$signal.node <- apv <= limit_rej
    
    if (message) {
        message("output the results ...")
    }
    out <- list(candidate_best = level_b, 
                output = outB,
                candidate_list = levels,
                level_info = level_info,
                FDR = limit_rej, 
                method = method,
                column_info = list("node_column" = node_column,
                                   "p_column"= p_column,
                                   "sign_column"= sign_column,
                                   "feature_column" = feature_column))
    return(out)
    
}

getCand <- function(tree, t = NULL,
                    score_data, node_column,
                    p_column, sign_column,
                    threshold = 0.05,
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
    
    # a list to store levels under different t 
    level_list <- vector("list", length(t))
    names(level_list) <- c(t)
    
    # columns: p value, sign, node
    p_col <- score_data[[p_column]]
    sign_col <- score_data[[sign_column]]
    node_col <- score_data[[node_column]]
    
    # paths
    path <- matTree(tree = tree)
    
    # for each t
    for (i in seq_along(t)) {
        if (message) {
            message("Searching candidates on t =  ", t[i], " ...")
        }
        # q
        name_q <- paste0("q_", t[i])
        score_data[[name_q]] <- ifelse(p_col > t[i], 0,
                                       1) * sign(sign_col)
        
        # add a column to store the result of search: keep
        if (any(colnames(score_data) == "keep")) {
            stop("The result will be output in the 'keep' column;
             Please use other name for the current 'keep' column.")
        }
        
        keep <- !is.na(p_col)
        node_keep <- score_data[[node_column]][keep]
        node_all <- printNode(tree = tree, type = "all")
        node_in <- node_all[["nodeNum"]][!node_all[["isLeaf"]]]
        node_in <- intersect(node_keep, node_in)
        
        # nodes with p value available
        leaf <- .pseudoLeaf(tree = tree, score_data = score_data, 
                            node_column = node_column,
                            p_column = p_column)
        
        
        # For an internal node, if more than half of its direct child nodes has
        # NA score, it would not be picked.
        chl_I <- findChild(tree = tree, node = node_in)
        sel_0 <- lapply(chl_I, FUN = function(x){
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            sum(!is.na(qx))/length(qx) > 0.5
        })
        sel_0 <- unlist(sel_0)
        node_0 <- node_in[sel_0]
        
        
        
        # For an internal nodes, if itself and all its descendant have q score
        # equals 1 or -1, pick the node
        br_I <- findDescendant(tree = tree, node = node_0, only.leaf = FALSE,
                               self.include = TRUE, use.alias = TRUE)
        
        sel_1 <- lapply(br_I, FUN = function(x){
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            abs(mean(qx, na.rm = TRUE)) == 1
        })
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
        
        level_list[[i]] <- c(setdiff(leaf, desd_2), node_2)
        
    }
    
    out <- list(candidate_list = level_list,
                score_data = score_data)
    return(out)
    
}

#' @importFrom TreeSummarizedExperiment matTree
#' @keywords internal
.pseudoLeaf <- function(tree, score_data, node_column, p_column) {
    mat <- matTree(tree = tree)
    nd <- score_data[[node_column]][!is.na(score_data[[p_column]])]
    exist_mat <- apply(mat, 2, FUN = function(x) {x %in% nd})
    
    ww <- which(exist_mat, arr.ind = TRUE)
    ww <- ww[order(ww[, 1]), , drop = FALSE]
    loc_leaf <- ww[!duplicated(ww[, 1]), ]
    leaf_0 <- unique(mat[loc_leaf])
    ind_0 <- lapply(leaf_0, FUN = function(x) {
        xx <- which(mat == x, arr.ind = TRUE)
        ux <- xx[!duplicated(xx), , drop = FALSE]
        y0 <- nrow(ux) == 1
        if (nrow(ux) > 1) {
            ux[, "col"] <- ux[, "col"] - 1
            y1 <- all(!exist_mat[ux])
            y0 <- y0 | y1
        }
        return(y0)
    })
    leaf_1 <- leaf_0[unlist(ind_0)]
    
    
    return(leaf_1)
}