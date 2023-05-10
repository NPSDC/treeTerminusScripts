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

estimatePThresh <- function(y, adjPval = 0.05, type = "BH") {
    library(fdrtool)
    pvalues <- c()
    if(is(y, "SummarizedExperiment"))
        pvalues <- mcols(y)[["pvalue"]]
    else{
        if(class(y) == "numeric")
            pvalues <- y
        else
            stop("invalid class")
    }

    if(type=="BH") {
        adPval <- which.min(abs(p.adjust(pvalues, method = "BH") - adjPval))
        return(pvalues[adPval])
    }
    else{
        pvalues <- pvalues[!is.na(pvalues)]
        fdrIn <- which.min(abs(fdrtool(pvalues, statistic="pvalue")[["lfdr"]] - adjPval))
        return(pvalues[fdrIn])
    }


}

checkGoUp <- function(parent, desc, children, signs, mInfRV, minInfRV, nodeSig, temp = F)
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
        descPVal <- nodeSig[desc]
        descNA <- which(is.na(descPVal))
        descPVal <- descPVal[setdiff(seq_along(descPVal), descNA)]

        if(all(!descPVal)) ## All are non signficant
        #if(any(!descPVal) | length(descNA) > 0) ## Any node is non significant
        {
            pIRV <- mInfRV[parent]
            cIRV <- mInfRV[children]
            # diff <- pIRV - mean(cIRV)
            # if(diff <= infDiff)
            # if(pIRV >= minInfRV)
            #     return(T)\
            if(temp) {
                if(all(cIRV >= minInfRV))
                    return(T)
                #print("infRV")
                return(F)
            }
            return(T)
        }
        else
        {
         #   print("significant")
            return(F)
        }

    }
    else
    {
        #print("sign")
        return(F)
    }
}

doIHW <- function(y, tree, alpha, iRVBin = 4, mCountBin = 4, nbins=NULL, inds = NULL)
{
    group <- !(is.null(iRVBin) | is.null(mCountBin))
    levels <- node.depth(tree, 2)
    levels <- ifelse(levels > 4, 5, levels)
    mCount <- rowMeans(assays(y)[["counts"]])
    mCCut = tryCatch( {
            cut(mCount, breaks = quantile(mCount, 0:mCountBin/mCountBin), include.lowest = T)
    },
    error=function(cond) {
        message("Changing bins to 2 since 4 gave an error")
        return(cut(mCount, breaks = quantile(mCount, 0:2/2), include.lowest = T))
    })
    df <- data.frame(IRVCut = cut(mcols(y)[["meanInfRV"]], breaks = quantile(mcols(y)[["meanInfRV"]], 0:iRVBin/iRVBin), include.lowest = T),
                     mCCut = mCCut,
                     levels, pvalue = mcols(y)[["pvalue"]], infRVs = mcols(y)[["meanInfRV"]], mCount = mCount)
    if(group)
        df[["group"]] <- factor(paste(df[["IRVCut"]], df[["mCCut"]]))

    ## Splitting the nodes into folds, such that each fold contains equal amount of root child txps and inner nodes
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
    # print(table(folds))

    ### For running on subset of y
    if(!is.null(inds))
    {
        if(!is.logical(inds) | length(inds) != nrow(y))
            stop("Inds either not logical or does not have the same length as y")
        df <- df[inds,]
        folds <- folds[inds]
        print(table(folds))
    }

    if(group){
        #print(table(df[["group"]]))
        resgroup <- ihw(pvalue ~ group, data = df, alpha = alpha, folds = folds,
                        covariate_type = "nominal")
        return(resgroup)
    }
    if(!is.null(nbins)) {
        resgroup <- ihw(pvalue ~ infRVs, data = df, nbins=nbins,
                        alpha = alpha, folds = folds,
                        covariate_type = "ordinal")
        return(resgroup)
    }
    if(!is.null(iRVBin)) {
        resgroup <- ihw(pvalue ~ IRVCut, data = df,
                        alpha = alpha, folds = folds,
                        covariate_type = "ordinal")
        return(resgroup)
    }
    if(!is.null(mCountBin)) {
        resgroup <- ihw(pvalue ~ mCCut, data = df,
                        alpha = alpha, folds = folds,
                        covariate_type = "ordinal")
        return(resgroup)
    }
}


runTreeTermAlphas <- function(tree, y, cond, minInfRV, pCutOff = 0.05, pChild = 0.05, corr = "IHW", runType = c("a"), alphas = c(0.01, 0.05, 0.1),
                              compPThresh = T, cSign = T, cores = 1, temp = T, file=NULL, minP = 0.70)
{
    library(qvalue)
    library(fdrtool)
    if(!corr %in% c("qvalue", "IHW", "BH", "lfdr"))
        stop(paste("Invalid correction type entered, should be either qvalue or IHW", corr))
    if(!runType %in% c("a", "b"))
        stop(paste("Invalid IHW type entered, should be either a or b", ihwType))

    sols <- vector("list", length(alphas))
    names(sols) <- as.character(alphas)
    if(runType == "b")
    {
        if (corr == "IHW") {
            resAlphas <- mclapply(alphas, function(alpha) doIHW(y, tree, alpha), mc.cores = cores)
            gc()
            resDfs <- lapply(resAlphas, function(res) {
                    df <- as.data.frame(res)
                    df <- cbind(df, inds = seq(nrow(y)))
                    df
                })
        }
        else {
            if(corr == "qvalue")
                resAlphas <- list(qvalue(mcols(y)[["pvalue"]]))
            else if(corr == "BH")
                resAlphas <- list(list(qvalues = p.adjust(mcols(y)[["pvalue"]], method="BH")))
            else {
                qvals <- rep(1, nrow(y))
                pvals <- mcols(yAll)[["pvalue"]]
                naInds <- is.na(pvals)
                lfdr <- fdrtool(pvals[!naInds], statistic="pvalue")[["lfdr"]]
                qvals[!naInds] <- lfdr
                resAlphas <- list(list(qvalues = qvals))

            }

            resDfs <- list(data.frame(pvalue=mcols(y)[["pvalue"]], adj_pvalue=resAlphas[[1]][["qvalues"]], inds = seq(nrow(y))))
        }
        nSol <- mclapply(seq_along(alphas), function(i) {
            j=i
            if(corr != "IHW")
                j=1
            runTreeTermAlpha(tree, y, cond, minInfRV, resDfs[[j]][["adj_pvalue"]], alphas[i], alphas[i], cSign = cSign, temp=temp,minP=minP)
        }, mc.cores = cores)

        for(i in seq_along(alphas))
        {
            sols[[i]][["candNodeO"]] <- which(nSol[[i]][["candNode"]])
            sols[[i]][["negNode"]] <- nSol[[i]][["negNode"]]
            sols[[i]][["negNodeO"]] <- which(nSol[[i]][["negNode"]])
            sols[[i]][["naNodeO"]] <- which(nSol[[i]][["naNode"]])
        }
    }
    else
    {
        nSol <- mclapply(alphas, function(alpha) {
            pThresh <- alpha
            if(compPThresh)
                pThresh <- estimatePThresh(y[1:length(tree$tip),], alpha)
            else
                pThresh <- alpha
            print(pThresh)
            runTreeTermAlpha(tree, y, cond, minInfRV, mcols(y)[["pvalue"]], pCutOff = pThresh, pChild = pThresh, cSign = cSign, temp=temp,minP=minP)
        }, mc.cores = cores)
        if(!is.null(file))
            save(nSol, file = file)
        if (corr == "IHW") {
            resAlphas <- mclapply(seq_along(alphas), function(i) doIHW(y, tree, alphas[i], inds = nSol[[i]][["nodesLooked"]]), mc.cores = cores)
            gc()
            resDfs <- lapply(seq_along(resAlphas), function(i)
                {
                    df <- as.data.frame(resAlphas[[i]])
                    df <- cbind(df, inds = which(nSol[[i]][["nodesLooked"]]))
                    df
                })
        }
        else {
            resAlphas <- lapply(seq_along(alphas), function(i) {
                nLooked <- nSol[[i]][["nodesLooked"]]
                if(corr=="qvalue")
                    qvalue(mcols(y)[["pvalue"]][nLooked])
                else
                    list(qvalues=p.adjust(mcols(y)[["pvalue"]][nLooked], method="BH"))
            })
            resDfs <- lapply(seq_along(resAlphas), function(i) {
                    nLooked <- which(nSol[[i]][["nodesLooked"]])
                    data.frame(pvalue=mcols(y)[["pvalue"]][nLooked], adj_pvalue=resAlphas[[i]][["qvalues"]], inds = nLooked)
                })
        }
        for(i in seq_along(alphas)) {
            sols[[i]][["candNodeO"]] <- intersect(which(nSol[[i]][["candNode"]]), resDfs[[i]][resDfs[[i]]$adj_pvalue <= alphas[i],"inds"])

            # remNegNodes <- setdiff(which(nSol[["candNode"]]), sols[[i]][["candNodeO"]])
            remNegNodes <- setdiff(which(nSol[[i]][["candNode"]]), sols[[i]][["candNodeO"]])
            if(sum(sols[[i]][["negNode"]][remNegNodes]) > 0)
                stop("Incorrect neg nodes ")
            # sols[[i]][["negNode"]] <- nSol[["negNode"]]
            sols[[i]][["negNode"]] <- nSol[[i]][["negNode"]]
            sols[[i]][["negNode"]][remNegNodes] <- T
            sols[[i]][["negNodeO"]] <- which(sols[[i]][["negNode"]])
            sols[[i]][["naNodeO"]] <- which(nSol[[i]][["naNode"]])
        }
        # nSol <- list(nSol)
    }
    for(i in seq_along(alphas))
    {
        j=i
        if(runType == "b" & corr != "IHW")
            j=1
        sols[[i]][["resIHW"]] <- resAlphas[[j]]
        sols[[i]][["resDf"]] <- resDfs[[j]]
        # if(ihwType=="a")
        #     j=1
        j = i
        sols[[i]][["nodesLooked"]] <- nSol[[j]][["nodesLooked"]]
        sols[[i]][["candNode"]] <- nSol[[j]][["candNode"]]
        sols[[i]][["pCut"]] <- nSol[[j]][["pCut"]]
        sols[[i]][["pChild"]] <- nSol[[j]][["pChild"]]
    }
    return(sols)
}

runTreeTermAlpha <- function(tree, y, cond, minInfRV, pvalue, pCutOff = 0.05, pChild = 0.05, cSign = T, temp=T,minP=0.70)
{
    nodeSig <- rep(T, nrow(y)) ### node significant or not
    nodeSig[pvalue > pChild] = F
    nodeSig[is.na(pvalue)] = NA ### Setting pvalues for nodes to NA

    candNode <- rep(F, nrow(y)) ### Nodes that would be the output
    nodesLooked <- rep(F, nrow(y)) ### Nodes that we have already gone over and wont be going over any further
    negNode <- rep(F, nrow(y)) ### Nodes that are output as negative
    naNode <- rep(F, nrow(y)) ### Nodes that are output as NA
    naLooked <- rep(F, nrow(y)) ### Nodes that are looked are output as NA
    nonSigLeaves <- intersect(seq_along(tree$tip.label), union(which(!nodeSig), which(is.na(nodeSig)))) ### we want to aggregate leaves with NA pvalues as well
    print(length(nonSigLeaves))
    sigLeaves <- intersect(seq_along(tree$tip.label), which(nodeSig))
    nodesLooked[sigLeaves] = T
    candNode[sigLeaves] = T
    if(cSign) signs <- computeSign(y, cond, minP = minP) else signs <- NULL
    #print(nonSigLeaves)
    for(n in nonSigLeaves) {
        #print(paste("i is", n))
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
            if(!temp) {
                if(mcols(y)[["meanInfRV"]][curNode] <= minInfRV)
                    break()
            }

            if(!checkGoUp(p, desc, child, signs, mcols(y)[["meanInfRV"]], minInfRV, nodeSig, temp=temp))
            #if(!checkGoUp(p, desc, child, signs, mInfRV, infDiff, nodeSig))
                break()
            foundCand <- T
            curNode <- p
        }
        #print(curNode)
        if(is.na(pvalue[curNode])) {
            if(!naLooked[curNode])
            {
                naNode[curNode] <- T
                naLooked[c(curNode,unlist(Descendants(tree, curNode, "all")))] <- T
            }
            if(curNode > length(tree$tip))
                naNode[unlist(Descendants(tree, curNode, "all"))] <- F
        }
        else {
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
        }
    }
    nodesLooked[!mcols(y)[["keep"]]] <- F
    naNode[unlist(Descendants(tree, sort(c(which(candNode), which(negNode))), "all"))] <- F
    return(list("candNode" = candNode, "negNode" = negNode, "naNode" = naNode, "nodesLooked" = nodesLooked, "pCut" = pCutOff, "pChild" = pChild))
}


meetCriteria <- function(yAll, tree, signs, mIRVCut, pCut, nInd) {
    #print(nInd)
    if(mcols(yAll)[nInd, "pvalue"] > pCut | is.na(mcols(yAll)[nInd, "pvalue"]))
        return(F)
    if(nInd < length(tree$tip)) {  ## leaf node that is significant
       if(is.na(mcols(yAll)[nInd,"pvalue"]))
           return(F)
       if(mcols(yAll)[nInd,"pvalue"] < pCut)
           return(T)
    }

    desc <- unlist(Descendants(tree, nInd,  'all')) ## Inner nodes
    children <- Descendants(tree, nInd,  'child')
    #if((all(signs[desc] >= 0) | all(signs[desc] <= 0)) & all(mcols(yAll)[desc,"meanInfRV"] > mIRVCut)) { ##same sign and minInfRV
    if((all(signs[desc] >= 0) | all(signs[desc] <= 0)) & all(mcols(yAll)[children,"meanInfRV"] > mIRVCut)) { ##same sign and minInfRV
        ##IS.NA()
            if(sum(is.na(mcols(yAll)[desc,"pvalue"])) > 0 | !all(mcols(yAll)[desc,"pvalue"] < pCut, na.rm=T)) {
                if(sum(is.na(mcols(yAll)[children,"pvalue"])) > 0)
                    return(T)
                if(!all(mcols(yAll)[children,"pvalue"] < pCut))  ##all children cant be significant
                    return(T)
            }
    }
    return(F)
}

findCNodes <- function(yAll, tree, ind, pCut, mIRVCut, signs) {
    # dNodes <- Descendants(tree, nrow(y)+1, "child")
    # dNodes <- dNodes[sapply(Descendants(tree, dNodes),length) > 1]
    # dNodes <- dNodes[sapply(Descendants(tree, dNodes, "all"), function(nodes) sum(mcols(yAll)[nodes,"qvalue"] < 0.1, na.rm=T)>1)]
    # sums <- sapply(Descendants(tree, dNodes, "all"), function(nodes) sum(mcols(yAll)[nodes,"qvalue"] < 0.1, na.rm=T))
    #
    # pNodes <- Descendants(tree, nrow(y)+1, "child")
    # pNodes <- pNodes[sapply(Descendants(tree, pNodes),length) > 1]
    # pNodes <- pNodes[sapply(Descendants(tree, pNodes), function(nodes) sum(mcols(yAll)[nodes,"qvalue"] < 0.1, na.rm=T)==1)]
    desc <- unlist(Descendants(tree, ind, "all"))
    sigNodes <- desc[which(mcols(yAll)[desc,"pvalue"] < pCut)]
    if(!is.na(mcols(yAll)[ind,"pvalue"])) {
        if(mcols(yAll)[ind,"pvalue"] < pCut & ind > length(tree$tip))
            sigNodes <- c(ind, sigNodes)
    }

    if(length(sigNodes) == 0)
        return(c())
    # names(mSigNodes) <- sigNodes

    sigInn <- sort(sigNodes[sigNodes > length(tree$tip)])
    if(length(sigInn) == 0)
        return(sigNodes)
    sigCons <- c() ##sigNodes that are non overlapping
    sigLeft <- sigInn

    #descSig <- Descendants(tree, sigInn, "all")
    #print(paste("sigLeft", sigLeft))
    while(length(sigLeft) > 0) {
        sigCons <- c(sigCons, sigLeft[1])
        sigLeft <- setdiff(sigLeft, c(sigLeft[1],unlist(Descendants(tree,sigLeft[1], "all"))))
    }

    sigCondD <- unlist(Descendants(tree, sigCons,"all"))
    sigCons <- c(sigCons, setdiff(sigNodes,c(sigCons,sigCondD)))
    # print(sigNodes)
    # print(sigCons)
    nodes <- lapply(sigCons, function(ind) findTNode(yAll, tree, ind, sigNodes, pCut, mIRVCut, signs))
    if(sum(is.na(unlist(nodes))) > 0)
      print(ind)
    return(unlist(nodes))
}

findTNode <- function(yAll, tree, ind, sigNodes, pCut, mIRVCut, signs) {
    if(length(intersect(c(ind,unlist(Descendants(tree, ind, "all"))), sigNodes)) == 0) ##Existing node does not have any descendants in the known significant node
        return(c())

    #print(sigNodes)
    if(meetCriteria(yAll, tree, signs, mIRVCut, pCut, ind))
        return(ind)
    children <- Descendants(tree, ind, "child")
    nodes <- c()
    for(child in children)
        nodes <- c(nodes,findTNode(yAll, tree, child, sigNodes, pCut, mIRVCut, signs))
    if(sum(is.na(nodes))) {
      print(ind)
    }
    return(nodes)
}

climbMax <- function(y, tree, mIRVCut, alpha, signs, minP=0.80, cores=1) {
    l <- length(tree$tip)
    pThresh <- estimatePThresh(y, alpha)
    cNodes <- Descendants(tree, l+1, "child")
    leaves <- cNodes[cNodes <= l]
    sigNodes <- leaves[which(mcols(y)[leaves, "pvalue"] <= pThresh)]
    innNodes <- setdiff(cNodes, leaves)
    descN <- Descendants(tree, innNodes, "all")
    descN <- innNodes[sapply(descN, function(nodes) sum(mcols(y)[nodes,"pvalue"] < pThresh, na.rm=T)>0)]
    sigs <- unlist(mclapply(descN, function(node) findCNodes(y, tree, node, pThresh, mIRVCut, signs), mc.cores=cores))
    sigNodes <- c(sigNodes, sigs)
    return(sigNodes)
}
