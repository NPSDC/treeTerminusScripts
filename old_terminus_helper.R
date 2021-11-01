library(fishpond)
library(fishpond)
library(SummarizedExperiment)

parseClustFile <- function(file, y)
{
    if(!file.exists(file))
        stop("Invalid file")
    df <- read.delim(file, header=F)
    groups <- lapply(strsplit(df$V1, split = ",", fixed = T), function(x) {
            inds <- match(unlist(x[2:length(x)]), rownames(y))
            #print(inds)
            if(sum(is.na(inds) > 0))
                stop("error ")
            inds
        })
    names(groups) <- as.character(nrow(y) + 1:length(groups))
    return(groups)
    
}

computeAggNodesU <- function(groups, nodeID, se_counts, group_inds = NULL) {
    performRowAgg <- function(counts, col_inds = NULL) {
        if(is.null(dim(counts)))
            stop("counts has to be matrix/dataframe")
        if(is.null(rownames(counts)))
            stop("rows must be named")
        
        if(is.null(col_inds))
            return(counts)
        
        else
        {
            
            df <- matrix(0, nrow = nrow(counts), ncol = length(col_inds), dimnames = list(rownames(counts)))
            for(i in seq_along(col_inds))
                df[,i] = rowMeans(counts[,col_inds[[i]]])
            return(df)
        }
    }
    performColAgg <- function(counts, row_inds = NULL)
    {
        if(is.null(dim(counts)))
            stop("counts has to be matrix/dataframe")
        if(is.null(row_inds))
            return(counts)
        if(is.null(names(row_inds)))
            stop("row indexes must be named")
        df <- matrix(0, nrow = length(row_inds), ncol = ncol(counts))
        for(i in seq_along(row_inds))
            df[i,] <- colSums(counts[row_inds[[i]],])
        rownames(df) <- names(row_inds)
        return(df)
    }
    
    if(!is.numeric(nodeID))
    {
        nodeID <- as.numeric(nodeID)
        if(sum(is.na(nodeID)) > 0)
            stop("Node ids contain a non numeric")
    }
    
    mat <- matrix(0, nrow=0, ncol=ncol(se_counts))
    
    leaves <- which(nodeID <= nrow(se_counts))
    innNodes <- which(nodeID > nrow(se_counts))
    
    lInds <- groups
    names(lInds) <- as.character(nodeID[innNodes])
    ls <- sapply(lInds, length)
    if(length(leaves) > 0)
        mat <- performColAgg(se_counts[nodeID[leaves],])
    if(length(innNodes) > 0)
        mat <- rbind(mat, performColAgg(se_counts, lInds))
    
    mat <- performRowAgg(mat, group_inds)       
    
    return(mat)
}

prepSwish <- function(y, inds, groups) 
{
    asList <- vector(mode = "list", length(assays(y)))
    #asList <- vector(mode = "list", 5)
    names(asList) <- assayNames(y)
    
    for(n in names(asList))
        asList[[n]] <- computeAggNodesU(groups, inds, assays(y)[[n]], NULL)
    
    y <- SummarizedExperiment(assays = asList, colData = colData(y), metadata = metadata(y))
    metadata(y)$infRepsScaled=F
    y
}

### All in y space
### mInds - number of missing txps in yAll compared to y, from mInd +1 groups will start
extractTxpGroup <- function(yAll, y, groups, diffTxpInds, mInds, txpSpace = c("y")) {
    dtes <- diffTxpInds[diffTxpInds <= mInds]
    dgroups <- rownames(yAll)[setdiff(diffTxpInds, dtes)]
    dtAll <- c(match(rownames(yAll)[dtes], rownames(y)), unlist(groups[dgroups]))
    return(dtAll)
}
    