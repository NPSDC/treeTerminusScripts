source("tree_helper_function.R")
library(IHW)
library(foreach)
library(doParallel)

### Differences between infRV of a group and its descendants
computeInfRVDiff <- function(tree, y)
{
    mIRVDiff <- rep(0.0, dim(y)[1])
    infRVs <- mcols(y)[["meanInfRV"]]
    mIRVDiff[seq_along(tree$tip.label)] <- infRVs[seq_along(tree$tip.label)]
    innNodes <- length(tree$tip.label)+seq(tree$Nnode)
    children <- Descendants(tree, innNodes, "children")
    for(i in seq_along(innNodes))
    #for(i in seq(10))
    {
        # print(mean(infRVs[children[[i]]]))
        # print(infRVs[innNodes[i]])
        mIRVDiff[innNodes[i]] = infRVs[innNodes[i]] - mean(infRVs[children[[i]]])
    }
    return(mIRVDiff)
}

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
    signDesc <- rep(1, length(c(parent, desc)))
    #print(signs)
    if(!is.null(signs))
        signDesc <- signs[c(parent, desc)]
    #print(signDesc)
    if(all(signDesc >= 0) | all(signDesc <= 0))  ##same sign of children
    {
        if(all(!nodeSig[desc])) ## All are non signficant
        {
            pIRV <- mInfRV[parent]
            cIRV <- mInfRV[children]
            diff <- pIRV - mean(cIRV)
            if(diff <= infDiff)
            #if(pIRV >= infDiff)
                return(T)
            print("irv")
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

doIHW <- function(y, tree, alpha, nbins=40, inds = NULL)
{
    levels <- node.depth(tree, 2)
    levels <- ifelse(levels > 4, 5, levels)
    df <- data.frame(IRVCut = cut(mcols(y)[["meanInfRV"]], breaks = quantile(mcols(y)[["meanInfRV"]], 0:8/8), include.lowest = T),
                     levels, pvalue = mcols(y)[["pvalue"]], infRVs = mcols(y)[["meanInfRV"]])
    
    rChildNodes <- Descendants(tree, length(tree$tip.label)+1, "children") ### Root child nodes
    folds <- rep(1, nrow(y))
    rParNodes <- rChildNodes[rChildNodes %in% seq_along(tree$tip.label)] ## Txps that directly map to root
    remChildNodes <- setdiff(rChildNodes, rParNodes)
    #ch <- cut(remChildNodes, breaks = quantile(remChildNodes, 0:4/4))
    ch <- cut(remChildNodes, breaks = quantile(remChildNodes, 0:5/5))
    ch[1] <- ch[2]
    d <- split(remChildNodes, ch)
    
    ch2 <- cut(rParNodes, breaks = quantile(rParNodes, 0:5/5))
    ch2[1] <- ch2[2]
    d2 <- split(rParNodes, ch2)
    for(i in seq_along(d))
        #folds[c(d[[i]],unlist(Descendants(tree, d[[i]],"all")))] <- i
        folds[c(d[[i]],unlist(Descendants(tree, d[[i]],"all")),d2[[i]])] <- i
    print(table(folds))
    if(!is.null(inds))
    {
        if(!is.logical(inds) | length(inds) != nrow(y))
            stop("Inds either not logical or does not have the same length as y")
        df <- df[inds,]
        folds <- folds[inds]
        print(table(folds))    
    }
    resgroup <- ihw(pvalue ~ infRVs, data = df, nbins=nbins,
                    alpha = alpha, folds = folds,
                    covariate_type = "ordinal")
}

runTreeTermAlphas <- function(tree, y, cond, infDiff, pCutOff = 0.05, pChild = 0.05, ihwType = c("b"), alphas = c(0.01, 0.05, 0.1), cSign = T, cores = 1)
{
    if(!ihwType %in% c("a", "b"))
        stop(paste("Invalid IHW type entered", ihwType))
    
    sols <- vector("list", length(alphas))
    names(sols) <- as.character(alphas)
    if(ihwType == "b")
    {
        resAlphas <- lapply(alphas, function(alpha) doIHW(y, tree, alpha))
        resDfs <- mclapply(resAlphas, function(res) {
                df <- as.data.frame(res)
                df <- cbind(df, inds = seq(nrow(y)))
                df
            }, mc.cores = cores)
        nSol <- mclapply(seq_along(resDfs), function(i) {
            runTreeTermAlpha(tree, y, cond, infDiff, resDfs[[i]][["adj_pvalue"]], alphas[i], alphas[i], cSign = cSign)
        }, mc.cores = cores)
        for(i in seq_along(alphas))
        {
            sols[[i]][["candNodeO"]] <- which(nSol[[i]][["candNode"]])
            sols[[i]][["negNode"]] <- nSol[[i]][["negNode"]]
            sols[[i]][["negNodeO"]] <- which(nSol[[i]][["negNode"]])
        }
    }
    else
    {
        # nSol <- runTreeTermAlpha(tree, y, cond, infDiff, mcols(y)[["pvalue"]], alphas, alphas)
        nSol <- mclapply(alphas, function(alpha) {
            runTreeTermAlpha(tree, y, cond, infDiff, mcols(yAll)[["pvalue"]], alpha, alpha, cSign = cSign)
        }, mc.cores = cores)
        resAlphas <- lapply(seq_along(alphas), function(i) doIHW(y, tree, alphas[i], inds = nSol[["nodesLooked"]]))
        resDfs <- mclapply(resAlphas, function(res) 
            {
                df <- as.data.frame(res)
                df <- cbind(df, inds = which(nSol[["nodesLooked"]]))
                df
            }, mc.cores = cores)
        for(i in seq_along(alphas)) {
            sols[[i]][["candNodeO"]] <- intersect(which(nSol[["candNode"]]), resDfs[[i]][resDfs[[i]]$adj_pvalue <= alphas[i],"inds"])
           # remNegNodes <- setdiff(which(nSol[["candNode"]]), sols[[i]][["candNodeO"]])
            remNegNodes <- setdiff(which(nSol[[i]][["candNode"]]), sols[[i]][["candNodeO"]])
            if(sum(sols[[i]][["negNode"]][remNegNodes]) > 0)
                stop("Incorrect neg nodes ")
            # sols[[i]][["negNode"]] <- nSol[["negNode"]]
            sols[[i]][["negNode"]] <- nSol[[i]][["negNode"]]
            sols[[i]][["negNode"]][remNegNodes] <- T
            sols[[i]][["negNodeO"]] <- which(sols[[i]][["negNode"]])
        }
        # nSol <- list(nSol)
    }
    for(i in seq_along(alphas))
    {
        j = i
        sols[[i]][["resIHW"]] <- resAlphas[[i]]
        sols[[i]][["resDf"]] <- resDfs[[i]]
        # if(ihwType=="a")
        #     j=1
        sols[[i]][["nodesLooked"]] <- nSol[[j]][["nodesLooked"]]
        sols[[i]][["candNode"]] <- nSol[[j]][["candNode"]]
    }
    return(sols)
}

runTreeTermAlpha <- function(tree, y, cond, infDiff, pvalue, pCutOff = 0.05, pChild = 0.05, cSign = T)
{
    nodeSig <- rep(T, nrow(y)) ### node signficant or not
    nodeSig[pvalue > pChild] = F
    
    tested <- rep(F, nrow(y)) ### tested or not (for IHW/qvalue/bonferroni)
    candNode <- rep(F, nrow(y)) ### Nodes that would be the output
    nodesLooked <- rep(F, nrow(y)) ### Nodes that we have already gone over and wont be going over any further
    negNode <- rep(F, nrow(y)) ### Nodes that are output as negative
    nonSigLeaves <- intersect(seq_along(tree$tip.label), which(!nodeSig))
    
    sigLeaves <- intersect(seq_along(tree$tip.label), which(nodeSig))
    nodesLooked[sigLeaves] = T
    candNode[sigLeaves] = T
    #tested[tree$tip.label] = T
    if(cSign) signs <- computeSign(y, cond) else signs <- NULL
    #print(nonSigLeaves)
    for(n in nonSigLeaves) {
        if(nodesLooked[n])
            next()
        curNode <- n
        foundCand <- F
        while(1)
        {
            if(curNode==(length(tree$tip.label)+1))
                break()
            p <- Ancestors(tree,curNode,"parent")
            if(nodesLooked[p])
                break()
            child <- unlist(Descendants(tree,p,"children"))
            desc <- unlist(Descendants(tree,p,"all"))
            if(any(nodesLooked[desc]))
                break()
            if(!checkGoUp(p, desc, child, signs, mcols(y)[["meanInfRV"]], infDiff, nodeSig))
            #if(!checkGoUp(p, desc, child, signs, mInfRV, infDiff, nodeSig))
                break()
            foundCand <- T
            curNode <- p
        }
        #print(curNode)
        negNode[curNode] <- T
        nodesLooked[c(curNode,unlist(Descendants(tree,curNode,"all")))] <- T ## not want to set right descendants of parent if not a candidate
        
        if(foundCand)
        {
            if(pvalue[curNode] <= pCutOff)
            {
                candNode[curNode] <- T
                negNode[curNode] <- F
            }
        }
        #print(nodesLooked)
    }
    return(list("candNode" = candNode, "negNode" = negNode, "nodesLooked" = nodesLooked))
}