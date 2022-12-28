library(fishpond)
library(fishpond)
library(SummarizedExperiment)

replaceMissingLength <- function(lengthMat, aveLengthSampGroup) {
    nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
    if (length(nanRows) > 0) {
        for (i in nanRows) {
            if (all(is.nan(lengthMat[i, ]))) {
                # if all samples have 0 abundances for all tx, use the simple average
                lengthMat[i, ] <- aveLengthSampGroup[i]
            } else {
                # otherwise use the geometric mean of the lengths from the other samples
                idx <- is.nan(lengthMat[i, ])
                lengthMat[i, idx] <- exp(mean(log(lengthMat[i, !idx]), na.rm = TRUE))
            }
        }
    }
    lengthMat
}


parseClustFile <- function(file, y, gsub=F)
{
    if(!file.exists(file))
        stop("Invalid file")
    df <- read.delim(file, header=F)
    groups <- lapply(strsplit(df$V1, split = ",", fixed = T), function(x) {
            if(gsub)
                inds <- match(unlist(gsub("\\.[0-9]+", "", x[2:length(x)])), rownames(y))
            else
                inds <- match(unlist(x[2:length(x)]), rownames(y))
            #print(inds)
            if(sum(is.na(inds) > 0))
                stop("error ")
            inds
        })
    names(groups) <- as.character(nrow(y) + 1:length(groups))
    return(groups)

}

computeOAggNodesU <- function(groups, nodeID, se_counts, group_inds = NULL) {
    performORowAgg <- function(counts, col_inds = NULL) {
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
    performOColAgg <- function(counts, row_inds = NULL)
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
        mat <- performOColAgg(se_counts[nodeID[leaves],])
    if(length(innNodes) > 0)
        mat <- rbind(mat, performOColAgg(se_counts, lInds))

    mat <- performORowAgg(mat, group_inds)

    return(mat)
}

prepOSwish <- function(y, inds, groups)
{
    asList <- vector(mode = "list", length(assays(y)))
    #asList <- vector(mode = "list", 5)
    names(asList) <- assayNames(y)

    for(n in names(asList)) {
        if(n!="length")
            asList[[n]] <- computeOAggNodesU(groups, inds, assays(y)[[n]], NULL)
        else {
            abundLength <- assays(y)[["abundance"]] * assays(y)[["length"]]
            weightedLength <- computeOAggNodesU(groups, inds, abundLength, NULL)
            lengthMat <- weightedLength / asList[["abundance"]]
            aveLengthSamp <- rowMeans(assays(y)[["length"]])
            gs <- c(lapply(seq(nrow(y)), function(i) i), groups)
            aveLengthSampGroup <- sapply(gs, function(ind) mean(aveLengthSamp[ind]))
            lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGroup)
            asList[[n]] <- lengthMat
        }

    }

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
