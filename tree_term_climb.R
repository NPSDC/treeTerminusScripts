getInfReps <- function(ys) {
    
    infRepError <- function(infRepIdx) {
        if (length(infRepIdx) == 0) {
            stop("there are no inferential replicates in the assays of 'y'")
        }
    }
    
    infRepIdx <- grep("infRep",assayNames(ys))
    infRepError(infRepIdx)
    infReps <- assays(ys)[infRepIdx]
    abind::abind(as.list(infReps), along=3)
}

computeSign <- function(y, x, pc = 5, minP = 0.7) {
    infRepsArray <- getInfReps(y)
    condition <- colData(y)[[x]]
    stopifnot(is.factor(condition))
    stopifnot(nlevels(condition) == 2)
    stopifnot(!anyNA(condition))
    
    dims <- dim(infRepsArray)
    cond1 <- condition == levels(condition)[1]
    cond2 <- condition == levels(condition)[2]
    diffs <- matrix(nrow=dims[1],ncol=dims[3])
    for (k in seq_len(dims[3])) {
        diffs[,k] <- log2(rowMeans(infRepsArray[,cond2,k]) + pc) -
            log2(rowMeans(infRepsArray[,cond1,k]) + pc)
    }
    signs <- rep(0, dims[1])
    pos <- diffs>0
    #print(head(pos))
    signs[rowMeans(pos) >= minP] = 1
    signs[rowMeans(!pos) >= minP] = -1
    return(signs)
}

checkGoUp <- function(parent, desc, children, signs, mInfRV, infDiff, nodeSig)
{
    if(parent > length(mInfRV) | any(desc > length(mInfRV)))
        stop("Invalid parent or child indexes")
    signDesc <- signs[desc]
    if(all(signDesc >= 0) | all(signDesc <= 0))  ##same sign of children
    {
        if(all(!nodeSig[desc])) ## All are non signficant
        {
            pIRV <- mInfRV[parent]
            cIRV <- mInfRV[children]
            diff <- pIRV - mean(cIRV)
            #if(diff <= infDiff)
            if(pIRV >= infDiff)
                return(T)
            #print("irv")
            return(F)
        }
        else
        {
            #print("significant")
            return(F)
        }
            
    }
    else
    {
        #print("sign")
        return(F)
    }
}

computeReps <- function(tree, y, cond, pCutOff = 0.05, infDiff)
{
    # nodeSig <- rep(T, nrow(y)) ### node signficant or not
    # nodeSig[mcols(y)[["pvalue"]] > pCutoff] = F
    
    tested <- rep(F, nrow(y)) ### tested or not (for IHW/qvalue/bonferroni)
    candNode <- rep(F, nrow(y)) ### Nodes that would be the output
    nodesLooked <- rep(F, nrow(y)) ### Nodes that we have already gone over and wont be going over any further
    
    nonSigLeaves <- intersect(seq_along(tree$tip.label), which(!nodeSig))
    sigLeaves <- intersect(seq_along(tree$tip.label), which(nodeSig))
    nodesLooked[sigLeaves] = T
    candNode[sigLeaves] = T
    #tested[tree$tip.label] = T
    #signs <- computeSign(y, cond)
    for(n in nonSigLeaves) {
        if(nodesLooked[n])
            next()
        curNode <- n
        foundCand <- F
        while(1)
        {
            p <- Ancestors(tree,curNode,"parent")
            if(nodesLooked[p])
                break()
            child <- unlist(Descendants(tree,p,"children"))
            desc <- unlist(Descendants(tree,p,"all"))
            #if(!checkGoUp(p, child, signs, mcols(y)[["meanInfRV"]], infDiff, nodeSig))
            if(!checkGoUp(p, desc, child, signs, mInfRV, infDiff, nodeSig))
                break()
            foundCand <- T
            curNode <- p
        }
        nodesLooked[c(p,curNode,unlist(Descendants(tree,curNode,"children")))] <- T ## not want to set right Descendants of parent if not a candidate
        if(foundCand)
        {
            if(nodeSig[curNode])
                candNode[curNode] = T
        }
        #print(nodesLooked)
    }
    return(candNode)
}