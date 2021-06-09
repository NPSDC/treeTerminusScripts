extractBouthNodes <- function(bouth_ob, se, type = "driver") {
    if(type=="driver")
        detNodes <- bouth_ob$results.by.node[bouth_ob$results.by.node$is.driver,]
    else
        detNodes <- bouth_ob$results.by.node[bouth_ob$results.by.node$is.detected,]
    levels <- ncol(detNodes) - 3
    det <- c()
    for(i in seq(nrow(detNodes)))
    {
        for(j in seq(levels))
        {
            if(detNodes[i,j]!=" " & detNodes[i,j+1]==" ")
            {
                det <- c(det, as.numeric(detNodes[i,j]))
                break()
            }

            else if(j==levels)
                det <- c(det, match(detNodes[i,j], rownames(y)))
        }    
    }
    return(det)
}


getTruth <- function(n, trueInds)
{
    truth <- rep(0,n)
    truth[trueInds] <- 1
    as.factor(truth)
}

### type in c(all, leaves)
computeTPFP <- function(tSig, sSig, y, logFC, tree = NULL, type = "all", pTab = F)
{
    # nodeDf <- bouth_ob$tree@node
    # nodeDf$id[nodeDf$id=="Root"] <- as.character(nrow(y)+1)
    # inds <- match(nodeDf$id, as.character(seq_along(logFC)))
    # pvalues <- bouth_ob$tree@test$pvalue[inds]
    if(type == "all")
    {
        truth <- getTruth(length(logFC), tSig)
        simTruth <- getTruth(length(logFC), sSig)    
    }
    else 
    {
        if(is.null(tree))
            stop("Tree cannot be null")
        tSig <- unique(unlist(Descendants(tree, tSig, type = "tips")))
        sSig <- unique(unlist(Descendants(tree, sSig, type = "tips")))
        truth <- getTruth(nrow(y), tSig)
        simTruth <- getTruth(nrow(y), sSig)
    }
    
    tab <- table(simTruth, truth)
    tpr <- tab[2,2]/colSums(tab)[2]
    fdr <- tab[2,1]/rowSums(tab)[2]
    if(pTab)
        print(tab)
    list(fdr=fdr, tpr=tpr)
}