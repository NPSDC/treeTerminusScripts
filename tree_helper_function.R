library(tximeta)
library(tximport)
library(ape)
library(phangorn)
library(fishpond)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)

## Given a tree of phylo, convert it into a data frame needed by BOUTH or treeRowData
getTreeDf <- function(tree) {
    nleaves <- length(tree$tip.label)
    maxDepth <- max(node.depth(tree, 2))
    
    df <- data.frame(matrix("unknown",nrow = nleaves, ncol = maxDepth, dimnames=list(c(tree$tip.label), paste("Group", c(1:maxDepth), sep="_"))))
    df[,maxDepth] <- tree$tip.label
    colnames(df)[maxDepth] <- "LEAF"
    des <- lapply(seq(1:nleaves), function(i) Ancestors(tree, i))
    df[,1] <- rep("Root", nrow(df))
    for(i in seq_along(des)) {
        l <- length(des[[i]])
        if(l > 1)
            df[i,2:l] <- rev(des[[i]])[2:l]
    }
    df
}

## Given a list of trees create one unified tree such all individual trees are children of one root node
mergeTree <- function(trees, updateInd = T, rowInd = T, tnames = NULL, se = NULL) {
    if(updateInd) {### Because of Salmon is 0 index and R is 1 index
        trees <- lapply(trees, function(tree) {
            tree$tip.label = as.character(as.numeric(tree$tip.label)+1)
            tree
        })    
    }
    
    nwk <- paste("(",paste(sapply(trees, function(tree) gsub(";", "", write.tree(tree))), collapse = ","), ");", sep = "")
    tree <- read.tree(text = nwk)
    if(!is.null(se))
        tree$tip.label <- rownames(se)[as.numeric(tree$tip.label)]
    return(tree)
}

## Given a tree run swish only on the filtered txps along with the group txps obtained by terminus
## type in c("union", "tree")
runSwishtree <- function(tree, ySwish, type) {
    if(!type %in% c("union", "tree"))
        stop(paste("invalid type", type))
    if(type == "union")
        tnames <- union(tree$tip.label, rownames(ySwish)[mcols(ySwish)$keep])
    if(type == "tree")
        tnames <- tree$tip.label
    ySwish <- ySwish[tnames,]
    mcols(ySwish)[,'keep'] <- TRUE
    set.seed(1)
    ySwish <- swish(ySwish, x="condition")
    return(ySwish)
}

## Given the updated swim seObject, update the indexes of tree that correspond to those indexes and also add the remaining leaf nodes 
mergeLeaves <- function(tree, ySwish, se) {
    #seInds <- match(rownames(ySwish), rownames(se))
    missing_txps <- setdiff(rownames(ySwish), tree$tip.label)
    if(length(missing_txps) > 0)
    {
        print(paste("Missing txps", length(missing_txps)))
        remLeaves <- paste(as.character(missing_txps), collapse = ",")
        nwk <- write.tree(tree)
        nwk <- substr(nwk, 2, nchar(nwk)-2)
        nwk <- paste("(", remLeaves, ",", nwk, ");", sep = "")
        tree <- read.tree(text = nwk)
    }
    #l <- length(setdiff(seInds, as.numeric(tree$tip.label)))
    # if(l != 0)
    #     stop(paste("Indexes do not match", l))
    
    # if(!is.null(se))
    #     tree <- convTree(tree, ySwish, se)
    return(tree)
}

### Arrange indices in the order of the tree that has been specified in the updated SE
convTree <- function(tree, ySwish, se) {
    tInds <- as.numeric(tree$tip.label)
    updatedInds <- match(rownames(se)[tInds], rownames(ySwish))
    tree$tip.label <- as.character(updatedInds)
    return(tree)
}

readTrueCounts <- function(files, sNames) {
    if(length(sNames) != length(files))
        stop("length of files not same as samples")
    ncol = length(length)
    txps = c()
    
    for(i in seq_along(files)){
        df <- read.table(files[i], header=T, stringsAsFactors = F)
        txps <- union(txps, df[,1])
    }
    cMat <- matrix(0, nrow=length(txps), ncol = length(sNames), dimnames = list(txps, sNames))
    for(i in seq_along(files)) {
        df <- read.table(files[i], header=T, stringsAsFactors = F)
        cMat[df[["Transcript"]],i] <- df[,"Counts"]
    }
    return(cMat)
}

createTreeSummarizedObject <- function(ySe, rTreeDf, rTree) {
    #rData <- cbind(rowData(ySe), rTreeDf)
    tse <- TreeSummarizedExperiment(assays = list(Count = assays(ySe)[["counts"]]),
                                        rowData = rTreeDf,
                                        colData = colData(ySe),
                                        rowTree = rTree,
                                        rowNodeLab = rTree$tip.label
                                    )
    return(tse)
}

## For a given node/s compute the aggreated sum across the leaves and then take the average across the replicates for each condition
computeAggNodes <- function(tree, nodeID, se_counts, group_inds = list(c(1),c(2))) {
    performColAgg <- function(counts, group_inds) {
        if(!is.null(dim(counts)))
            counts <- colSums(counts)
        vals <- sapply(group_inds, function(x) mean(counts[x]))
        vals
    }
    
    if(!is.numeric(nodeID))
    {
        nodeID <- as.numeric(nodeID)
        if(sum(is.na(nodeID)) > 0)
            stop("Node ids contain a non numeric")
    }
    df <- matrix(0, nrow = length(nodeID), ncol = length(group_inds))
    leaves <- which(nodeID <= nrow(se_counts))
    innNodes <- which(nodeID > nrow(se_counts))
    
    for(i in seq_along(leaves))
        df[i,] <- performColAgg(se_counts[i,], group_inds)
    
    lInds <- Descendants(tree, nodeID[innNodes], type = "tips")
    if(is.null(i))
        i <- 0
    for(j in seq_along(lInds))
        df[i+j,] <- performColAgg(se_counts[lInds[[j]],], group_inds)
    return(df)
}

