library(phangorn)
library(SummarizedExperiment)
library(fishpond)
computeInfRV <- function (y, pc = 5, shift = 0.01, rMean=F, addOne = F, subOne = F) 
{
    infReps <- assays(y)[grep("infRep", assayNames(y))]
    infReps <- abind::abind(as.list(infReps), along = 3)
    if(addOne)
        infReps <- infReps + 1
    infMean <- apply(infReps, 1:2, mean)
    if(subOne)
        infMean <- infMean - 1
    infVar <- apply(infReps, 1:2, var)
    
    assays(y)[["mean"]] <- infMean
    assays(y)[["variance"]] <- infVar
    
    InfRV <- pmax(infVar - infMean, 0)/(infMean + pc) + shift
    if(rMean)
        return(rowMeans(InfRV))
    return(list(infRV=InfRV, mean=assays(y)[["mean"]], variance=assays(y)[["variance"]]))
}

findOptSum <- function(tree, spl, node, lengths = NULL) {
    if(globArr[node]!=-100)
        return(globArr[node])
    if(node <= length(tree$tip)) {
        globArr[node] <<- spl[node]
        return(spl[node])
    }
    children <- Descendants(tree, node, "child")
    if(!is.null(lengths))
        v <- min(spl[node]*lengths[node], sum(sapply(children, function(child) findOptSum(tree, spl, child, lengths))))
    else
        v <- min(spl[node], sum(sapply(children, function(child) findOptSum(tree, spl, child))))
    globArr[node] <<- v
    return(v)
}

# optMeans <- rep(-100,nrow(yAll))
# counts <- rep(-100,nrow(yAll))
findOptMean <- function(tree, spl, node) {
    if(optMeans[node] != -100)
        return(c(optMeans[node]*counts[node], counts[node]))
    if(node <= length(tree$tip)) {
        optMeans[node] <<- spl[node]
        counts[node] <<- 1
        return(c(spl[node], 1))
    }
    
    s <- 0
    n <- 0
    children <- Descendants(tree, node, "child")
    for(child in children) {
        agg <- findOptMean(tree, spl, child)
        #print(agg)
        s <- s + agg[1]
        n <- n + agg[2]
    }
    if(s/n < spl[node]) {
        optMeans[node] <<- s/n
        counts[node] <<- n
        return(c(s,n))
    }
    optMeans[node] <<- spl[node]
    counts[node] <<- 1
    return(c(spl[node], 1))
}

findCuts <- function(tree, vals, spl, node, length = NULL) {
    if(!is.null(length)) {
        if(vals[node] == spl[node]*length[node])
            return(node)
    }
    else{
        if(vals[node] == spl[node])
            return(node)    
    }
    
    children <- Descendants(tree, node, "child")
    return(unlist(sapply(children, function(child) findCuts(tree, vals, spl, child, length))))
}


findMaxSum <- function(tree, mv, node, lengths = NULL) {
    if(globArr[node] != -100)
        return(globArr[node])
    if(node <= length(tree$tip)) {
        globArr[node] <<- mv[node]
        return(mv[node])
    }
    children <- Descendants(tree, node, "child")
    if(!is.null(lengths))
        v <- max(mv[node]*lengths[node], sum(sapply(children, function(child) findMaxSum(tree, mv, child, lengths))))
    else
        v <- max(mv[node], sum(sapply(children, function(child) findMaxSum(tree, mv, child))))
    globArr[node] <<- v
    return(v)
}

getLog2FC <- function(infRepsArray, condition, pc=5, array=FALSE) {
    dims <- dim(infRepsArray)
    cond1 <- condition == levels(condition)[1]
    cond2 <- condition == levels(condition)[2]
    diffs <- matrix(nrow=dims[1],ncol=dims[3])
    for (k in seq_len(dims[3])) {
        diffs[,k] <- log2(rowMeans(infRepsArray[,cond2,k]) + pc) -
            log2(rowMeans(infRepsArray[,cond1,k]) + pc)
    }
    if (array) {
        return(diffs)
    }
    # median over inferential replicates
    rowMedians(diffs)
}

