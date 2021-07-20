library(tximeta)
library(tximport)
library(ape)
library(phangorn)
library(fishpond)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(parallel)
library(DESeq2)

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
    ySwish <- computeInfRV(ySwish)
    return(ySwish)
}

## Given the seObject, create leaves for the txps in seObject missing in the tree and then arrange the rows of seObject in the order of leaves of the trees
## (for the index version refer previous commit on github)
mergeLeaves <- function(tree, ySwish) {
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
    ySwish <- ySwish[tree$tip.label,]
    return(list("tree" = tree, "ySwish" = ySwish))
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
## Will return matrix, with leaves first followed by inner nodes
computeAggNodesU <- function(tree, nodeID, se_counts, group_inds = NULL) {
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
    
    lInds <- Descendants(tree, nodeID[innNodes], type = "tips")
    names(lInds) <- as.character(nodeID[innNodes])
    ls <- sapply(lInds, length)
    
    if(length(leaves) > 0)
        mat <- rbind(mat, performColAgg(se_counts[nodeID[leaves],]))
    if(length(innNodes) > 0)
        mat <- rbind(mat, performColAgg(se_counts, lInds))
    mat <- performRowAgg(mat, group_inds)       
  
    return(mat)
}

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

prepSwish <- function(tree, seObLeaves) 
{
    innNodes <- nrow(seObLeaves)+1:tree$Nnode
    asList <- vector(mode = "list", length(assays(seObLeaves)))
    #asList <- vector(mode = "list", 5)
    names(asList) <- assayNames(seObLeaves)
    for(n in names(asList))
        asList[[n]] <- computeAggNodesU(tree, c(1:nrow(seObLeaves),innNodes), assays(seObLeaves)[[n]])
    #rData <- rowData(swishObLeaves)
    
    y <- SummarizedExperiment(assays = asList, colData = colData(seObLeaves), metadata = metadata(seObLeaves))
    metadata(y)$infRepsScaled=F
    y
}

library(data.table)
createDESeq2Ob <- function(indsIV, aggCounts, y, cores=6)
{
  if(length(unlist(indsIV)) != nrow(aggCounts) | any(duplicated(unlist(indsIV))))
    stop("all indexes not covered")
  ddsL <- mclapply(indsIV, function(inds) {
    d <- DESeqDataSetFromMatrix(countData = round(aggSimCountsNodes[inds,]),
                                colData = colData(y),
                                design = ~ condition)
    d <- DESeq(d, minReplicatesForReplace=Inf)
  }, mc.cores = cores)
  names(ddsL) <- names(indsIV)
  
  resL <- mclapply(seq_along(indsIV), function(i) {
    res <- results(ddsL[[i]], name="condition_2_vs_1", independentFilter = F, cooksCutoff = Inf)
    res$node <- indsIV[[i]]
    as.data.frame(res)
  }, mc.cores = cores)
  resL <- as.data.frame(rbindlist(resL))
  resL <- resL[order(resL$node),]
  rownames(resL) <- rownames(aggCounts)
  return(list("dds"=ddsL,"res"=resL))
}

createSwishOb <- function(indsIV, y, cores=6) {
  yAll <- prepSwish(tree, y)
  yAll <- scaleInfReps(yAll, lengthCorrect = T)
  yAll <- labelKeep(yAll)
  mcols(yAll)$keep=T
  yL <- mclapply(indsIV, function(inds)
    {
    swish(yAll[inds,], x = "condition")
  }, mc.cores = cores)
  names(yL) <- names(indsIV)
  
  sR <- mclapply(seq(indsIV), function(i)
    {
    swishRes <- data.frame(mcols(yL[[i]]))
    swishRes$node <- indsIV[[i]]
    swishRes <- swishRes[,-c(2)] ### Removing keep
  }, mc.cores = cores)
  gc()
  sR <- as.data.frame(rbindlist(sR))
  sR <- sR[order(sR$node),]
  rownames(sR) <- rownames(yAll)
  return(list("yL"=yL, "df"=sR))
}

createCand <- function(tree, resL, cores = cores, type = "deseq2")
{
  logFC <- "log2FoldChange"
  if(type=="swish")
    logFC <- "log2FC"
  candDeseqL <- mclapply(resL, function(res) {
                        getCand(tree,
                        score_data = res,
                        node_column = "node",
                        p_column = "pvalue",
                        sign_column = logFC,
                        message=T)
    }, mc.cores = cores)
  gc()
  return(candDeseqL)
}

getTruth <- function(n, trueInds)
{
  truth <- rep(0,n)
  truth[trueInds] <- 1
  as.factor(truth)
}

### type in c(all, leaves)
computeTPFP <- function(tSig, sSig, y, logFC = NULL, tree = NULL, type = "all", pTab = F)
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

findDriverGround <- function(detNodes, tree, logFCs)
{
    nLeaves <- length(tree$tip.label)
    detNodes <- sort(detNodes)
    innNodes <- detNodes[detNodes > nLeaves]
    notDriver <- c()
    driver <- c()
    detNodes <- c(innNodes, detNodes[detNodes <= nLeaves])
    desc <- Descendants(tree, detNodes, type = "tips")
    
    for(i in seq_along(detNodes))
    {
        if(detNodes[i] %in% notDriver)
          next()
        fcs <- logFCs[desc[[i]]]
        if(sum(fcs > 0) > 0 & sum(fcs < 0) > 0) {
          next()
        }
        notDriver <- c(notDriver, unlist(Descendants(tree, detNodes[i], "all")))
        driver <- c(driver, detNodes[i])
    }
    return(list(driver, notDriver))
}

findDistance <- function(tree, fpNodes, trueNodes)
{
    Anc <- Ancestors(tree, fpNodes, "all")
    Desc <- Descendants(tree, fpNodes, "children")
    inds <- c()
    for(i in seq_along(fpNodes))
    {
      #print(length(fpNodes))
      aInd <- which(Anc[[i]] %in% trueNodes)
      dInd <- sum(Desc[[i]] %in% trueNodes)
      if(length(aInd) == 1)
        inds <- c(inds,aInd)
    #  if()
      else{
        if(dInd == 1)
          inds <- c(inds, -1)
        else
        {
          lev = -1
          while(1) {
            lev = lev - 1
            desc <- unlist(Descendants(tree, Desc[[i]], "children"))
            if(sum(desc %in% trueNodes) == 1) {
              inds <- c(inds, lev)
              break()
            }
            if(sum(desc %in% seq(length(tree$tip))) == length(desc)) {
              print(paste(i,'1'))
              inds <- c(inds, NA)
              break()
            }
              
          }
        }
      }
      
    }
    return(inds)
}

###Figure out the list thing with ancestors/descendants based on the *input* being a number/vector and output being a *vector/number*
###Check the above w.r.t NA nodes as well
###For the above check whether my driver nodes are in actual driver nodes
###Come up with a metrc for error (include a way to include depth and the other nodes)

### Compute metrics w.r.t output predicted by a method as baseline
computeMetOut2 <- function(logFCTrue, methNodes, lfcThresh = 0.067)
{
  tNodes <- which(abs(logFCTrue) >= lfcThresh)
  tp = length(intersect(methNodes, tNodes))/length(methNodes)
  fp = length(setdiff(methNodes, tNodes))/length(methNodes)
  return(list("tp" = tp, "fp" = fp))
}

computeMetOut <- function(detNodes, logFCTrue, negNodes = NULL, nodesLooked = NULL, tree = NULL, onlyDet = F, lfcThresh = 0.067, res = T)
{
  if(onlyDet)
  {
    if(is.null(tree))
      stop("Tree cannot be null")
    if(!is.null(negNodes))
      stop("neg nodes should be null")
    descLeaves <- unlist(Descendants(tree, detNodes, "tips"))
    nodesLooked <- sort(c(setdiff(seq_along(tree$tip.label), descLeaves), detNodes))
  }
  if(!is.null(negNodes))
    nodesLooked <- sort(c(negNodes, detNodes))
  
  logFCTrue <- logFCTrue[nodesLooked]
  allNodes <- seq_along(nodesLooked)
  names(allNodes) <- as.character(nodesLooked)
  sNodes <- allNodes[as.character(detNodes)]
  tpNodes <- which(abs(logFCTrue) >= lfcThresh)
  tp = length(intersect(sNodes, tpNodes))/length(sNodes)
  fp = length(setdiff(sNodes, tpNodes))/length(sNodes)
  
  print(paste("tp", tp))
  met <- computeTPFP(tpNodes, sNodes, y = NULL, logFC = logFCTrue, type = "all")
  if(res)
    return(met)
  
  fps <- as.numeric(names(allNodes)[setdiff(sNodes, tpNodes)])
  fns <- as.numeric(names(allNodes)[setdiff(tpNodes, sNodes)])
  tps <- as.numeric(names(allNodes)[tpNodes])
  tns <- as.numeric(names(allNodes)[setdiff(allNodes, tpNodes)])
  return(list(fps = fps, fns = fns, tps = tps, tns = tns))
}

computeTDepth <- function(tree, methNodes)
{
  return(sum(node.depth(tree, 2)[methNodes])/length(methNodes))
}