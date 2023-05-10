library(tximeta)
library(tximport)
library(ape)
library(phangorn)
library(fishpond)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(parallel)
library(DESeq2)

### Correcting for batch effects
corrBESwish <- function(se, samples, fit, condition = "condition")
{
  y <- se[,!(colnames(se) %in% samples)]
  if(is(y[[condition]], "factor"))
          y[[condition]] <- droplevels(y[[condition]])
  y[[condition]] <- factor(y[[condition]])
  y <- scaleInfReps(y, saveMeanScaled=TRUE)
  
  infRepIdx <- grep("infRep",assayNames(y),value=TRUE)
  nreps <- length(infRepIdx)
  
  mm <- model.matrix(formula(paste("~", condition)), colData(y))
  
  pc <- .1
  for (k in seq_len(nreps)) {
    logInfRep <- log(assay(y, infRepIdx[k]) + pc)
    logInfRep <- limma::removeBatchEffect(
      logInfRep,
      covariates=fit[["sv"]],
      design=mm)
    assay(y, infRepIdx[k]) <- exp(logInfRep)
  }
  gc()
  
  y <- swish(y, x=condition)
  y
}

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
runSwishTree <- function(tree, ySwish, type) {
    if(!type %in% c("union", "tree"))
        stop(paste("invalid type", type))
    if(type == "union")
        tnames <- union(tree$tip.label, rownames(ySwish)[mcols(ySwish)$keep])
    if(type == "tree")
        tnames <- tree$tip.label
    ySwish <- ySwish[tnames,]
    #mcols(ySwish)[,'keep'] <- TRUE
    ySwish <- swish(ySwish, x="condition")
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
        # if(is.null(dim(counts)))
        #   stop("counts has to be matrix/dataframe")
        if(is.null(row_inds)) {
          if(is.null(dim(counts)))
            return(matrix(counts, nrow=length(counts), ncol=1))
          return(counts)
        }
            
        if(is.null(names(row_inds)))
            stop("row indexes must be named")
        
        sInd <- F ###if not matrix then column sums
        
        if(dim(counts)[2] == 1)
          sInd <- T
        df <- matrix(0, nrow = length(row_inds), ncol = ncol(counts))
        for(i in seq_along(row_inds)){
     #     print(sInd)
          if(!sInd)
            df[i,] <- colSums(counts[row_inds[[i]],])
          else
            df[i,] <- sum(counts[row_inds[[i]]])
        }
          
        rownames(df) <- names(row_inds)
        colnames(df) <- colnames(counts)
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
    #print(dim(se_counts[nodeID[leaves],]))
    
    if(length(leaves) > 0) {
      mat <- rbind(mat, performColAgg(se_counts[nodeID[leaves],]))
    }
    
    if(length(innNodes) > 0) {
      mat <- rbind(mat, performColAgg(se_counts, lInds))
      #print("ss")
    }
        
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

### Associate a group with every gene
getGeneGroup <- function(tree, t2GDf)
{
  getGenes <- function(node)
  {
    txps <- tree$tip[unlist(Descendants(tree, node, "tip"))]
    genes <- t2GDf[txps,"gene_id"]
    unique(genes)
  }
  getNode <- function(g, cur_root, inds)
  {
    genes <- getGenes(cur_root)
    ug = unique(genes)
    if(length(ug) == 1)
    {
      if(ug == g)
        return(c(inds, cur_root))
    }
    children <- Descendants(tree, cur_root, "child")
    genes <- lapply(children, getGenes)
    nodeInds <- c()
    gInd <- which(sapply(genes, function(gene) g %in% gene))
    if(length(gInd) == 0)
      return(-1)
    #print(children)
    children <- children[gInd]
    for(ch in children)
      inds <- getNode(g, ch, inds)
    return(inds)
  }

  innNodes <- Descendants(tree, length(tree$tip.label)+1, "child")
  # innNodesAll <- unlist(Descendants(tree, innNodesAll, "all"))
  # innNodesM <- rep(F, length(innNodesAll))
  # names(innNodesM) <- as.character(innNodesAll)
  #innNodes <- innNodes[innNodes >= length(tree$tip.label)] ### All children of root that is not a leaf (aka highest)
  txpAll <- Descendants(tree, innNodes, "tips")
  genesAll <- unique(t2GDf[tree$tip[unique(unlist(txpAll))],"gene_id"])
  geneTGroup <- vector(mode = "list", length(genesAll))
  names(geneTGroup) <- genesAll
  nGenes <- c()
  for(i in seq_along(innNodes))
  {
    # if(innNodesM[as.character(innNodes[i])])
    #   next()
    print(i)
    if(innNodes[i] <= length(tree$tip))
    {
      gene <- t2GDf[tree$tip[innNodes[i]], "gene_id"]
      geneTGroup[[gene]] <- c(geneTGroup[[gene]], innNodes[i])
   #   innNodesM[names(innNodes)[i]] <- T
    }
    else
    {
      genes <- getGenes(innNodes[i])
      for(g in genes)
      {
        inds <- getNode(g, innNodes[i], geneTGroup[[g]])
        if(-1 %in% inds)
          print(paste("could not find a node for", g))
        desc <- Descendants(tree, inds, "all")
      #  innNodesM[as.character(unlist(desc))] <- T
        geneTGroup[[g]] <- c(geneTGroup[[g]], inds)
      }
    }
  }
  geneTGroup <- lapply(geneTGroup, unique)
  return(geneTGroup)
}

plotGene <- function(gene, y, gy, gg, tree, type = c("txp"), inds_consid)
{
  if(!type %in% c("txp", "sub_gene", "gene"))
    stop("Invalid name")
  library(fishpond)
  nodeInds <- gg[[gene]]
  nodeCons <- intersect(nodeInds, inds_consid)
  nodesNotCons <- setdiff(nodeInds, nodeCons)
  nodesO <- c()
  if(length(nodesNotCons) != 0)
  {
      ch <- unlist(Descendants(tree, nodesNotCons, "child"))
      
      while(sum(ch %in% inds_consid) != length(ch))
      {
        nodesO <- c(nodesO, ch[ch %in% inds_consid])
        rem <- ch[!(ch %in% inds_consid)]
        ch <- unlist(Descendants(tree, rem, "child"))
      }
      nodesO <- unique(c(nodesO, ch[ch %in% inds_consid]))
  }
  nodeInds <- c(nodesO, nodeCons)
  # print(nodeInds)
  # if(length(nodeInds) == 1)
  # {
  #   if(nodeInds > length(tree$tip))
  #     nodeInds <- Descendants(tree, nodeInds, "child")
  # }
  #plot all inner nodes
  #plot txp with max qvalue
  
  txps <- nodeInds[nodeInds <= length(tree$tip)]
  innNodes <- nodeInds[nodeInds > length(tree$tip)]
  nodeInds <- setdiff(nodeInds, txps)
  

  if(type == "txp")
  {

    if(length(innNodes) > 0)
      txps <- c(txps,unlist(Descendants(tree, innNodes, "tip")))
    row <- ifelse(length(txps) == 1, 1, 2)
    x <- ifelse(length(txps) <= 2, 1, 2)
    if(length(txps) > 3)
      txps <- match(rownames(y)[txps][order(mcols(y)[txps,"qvalue"])[1:4]], tree$tip)
      
    par(mfrow=c(row,x))
    print(txps)
    for(i in seq_along(txps))
      plotInfReps(y, txps[i], x="condition", main = paste("Transcript ", tree$tip[txps[i]]))
    #dev.off()
  }
  
  if(type == "sub_gene")
  {
    row <- ifelse(length(nodeInds) == 1, 1, 2)
    x <- ifelse(length(nodeInds) <= 2, 1, 2)
    if(length(nodeInds) > 4)
      nodeInds <- rownames(y)[nodeInds][order(mcols(y)[nodeInds,"qvalue"])[1:4]]
    par(mfrow=c(row,x))
    print(nodeInds)
    for(i in seq_along(nodeInds))
      plotInfReps(y, nodeInds[i], x="condition", main = paste("Sub gene group ", nodeInds[i], " with ", length(Descendants(tree, nodeInds[i], "tip")[[1]]), 
                                                              "transcripts"))
    #dev.off()
  }
  if(type == "gene")
  {
    par(mfrow=c(1,1))
    plotInfReps(gy, gene, x = "condition", main = paste("Gene ", gene))
  }
}


### Get a list of genes and their associated groups given by tree climbing algorithm
getGeneGroup <- function(tree, nodes, mapDf) {
  txpInds <- Descendants(tree, nodes, "tip")
  genes <- lapply(txpInds, function(x) {
      txps <- tree$tip[x]
      unique(mapDf[txps, "GENEID"])
    })
  allG <- unique(unlist(genes))
  geneGroups <- vector(mode="list", length=length(allG))
  names(geneGroups) <- allG
  for(i in seq_along(genes)) {
    print(i)
    for(g in genes[[i]])
      geneGroups[[g]] <- c(geneGroups[[g]], nodes[i])
  }
  return(geneGroups)
}
  
getGeneCombPvalue <- function(geneGroups, y, method = c("harmonic"), correct = c("BH")) {
  library(metap)
  library(harmonicmeanp)
  library(qvalue)
  
  if(!method %in% c("harmonic", "fisher"))
    stop("incorrect method entered")
  if(!is.null(correct))
  {
    if(!correct %in% c("BH","qvalue"))
      stop("incorrect correction enter")  
  }
  L <- length(unlist(geneGroups))
  print(length(geneGroups))
  pvalues <- sapply(geneGroups, function(g) {
    pvals <- mcols(y)[g,"pvalue"]
    if(sum(is.na(pvals)) == length(pvals))
      return(NA)
    pvals <- pvals[!is.na(pvals)]
    if(method=="harmonic") {
      if(is.null(correct))
        p.hmp(pvals, L=L)
      else
        p.hmp(pvals, L=length(pvals))
    }
      
    else {
      if(length(pvals) > 1)
        sumlog(pvals, log.p=F)[["p"]]
      else
        pvals
    }
      
  })
  pvalues <- pvalues[!is.na(pvalues)]
  print(length(pvalues))
  if(!is.null(correct)) {
    if(correct=="BH")
      pvalues <- p.adjust(pvalues, method = "BH")
    else
      pvalues <- qvalue(pvalues)[["qvalues"]]
  }
   return(pvalues) 
}

fixDf <- function(df) {
  library(tidyverse)
  df <- data.frame(lapply(df, unlist))
  df_gather <- df %>% gather(FDR, Value, 3:5) # the columns are not variables, but levels, so gather them
  df_metrics <- df_gather %>% spread(Metric, Value) # we want TPR and FDR as cols
  # can't do everything in one go bc of the following line, where we add back e.g. FDR_0.01
  df_metrics$level <- sub("FDR_","",df_gather$FDR[2 * 1:nrow(df_metrics)])
  return(df_metrics)
}

### root of the mean (over features) of squared differences across the mean (over samples) log2 abundance per condition
computeDeseqInf <- function(y, tNodes = c(), retTPM = F) {
  tpm <- assays(y)[["abundance"]]
  sf <- DESeq2::estimateSizeFactorsForMatrix(tpm)
  logScaledTpm <- log2(t(t(tpm)/sf)+1)
  colData(y)[["condition"]] <- as.factor(colData(y)[["condition"]])
  cond1 <- colData(y)[["condition"]] == levels(colData(y)[["condition"]])[1]
  cond2 <- colData(y)[["condition"]] == levels(colData(y)[["condition"]])[2]
  
  mC1 <- rowMeans(logScaledTpm[,cond1]) ##mean condition 1
  mC2 <- rowMeans(logScaledTpm[,cond2]) ##mean condition 2
  diff <- (mC2 - mC1)^2
  if(class(tNodes) == "list")
    tNodes <- c(tNodes[["candNodeO"]], tNodes[["negNodeO"]], tNodes[["naNodeO"]])
  else
    tNodes <- tNodes
  if(retTPM)
    return(diff)
  return(sqrt(mean(diff[tNodes])))
}

### root of the sum (over features) of weighted squared log2 fold change of fishpond-scaled counts, weighting by the InfVar(LFC)
computeWeightedLFC <- function(y, tNodes, type = "sum") {
  compLFC <- function(reps, conds) {
    mC1 <- rowMeans(reps[,conds==1]) ##mean condition 1
    mC2 <- rowMeans(reps[,conds==2]) ##mean condition 2
    lfc <- log2(mC2+1) - log2(mC1+1)
    return(lfc)
  }
  if(class(y) == 'list') {
    lfcCounts <- y[[1]]
    varThresh <- y[[2]]
    lfcVar <- y[[3]]  
    lfcVar[lfcVar < varThresh] = varThresh[lfcVar < varThresh] ## setting threshold for the feature to expected variance  
    if(class(tNodes) == "list")
      nodes <- c(tNodes[["candNodeO"]], tNodes[["negNodeO"]], tNodes[["naNodeO"]])
    else
      nodes <- tNodes
    if(type == "mean")
      return(sqrt(mean(lfcCounts[nodes]^2/lfcVar[nodes])))
    return(sqrt(sum(lfcCounts[nodes]^2/lfcVar[nodes])))
  }  
  lfcCounts <- compLFC(assays(y)[["counts"]], colData(y)[["condition"]]) ## count LFC
  #varThresh <- log2(exp(1))^2*(1/(rowMeans(assays(y)[["counts"]][,colData(y)[["condition"]]==1])+1) + 1/(rowMeans(assays(y)[["counts"]][,colData(y)[["condition"]]==2])+1)) ## threshold for variance
  varThresh <- log2(exp(1))^2*(1/(rowSums(assays(y)[["counts"]][,colData(y)[["condition"]]==1])+1) + 1/(rowSums(assays(y)[["counts"]][,colData(y)[["condition"]]==2])+1)) ## threshold for variance
  
  lfcReps <- sapply(assayNames(y)[grepl("infRep", assayNames(y))], function(rep) compLFC(assays(y)[[rep]], colData(y)[["condition"]])) ## lfc for inf replicates
  lfcVar <- rowVars(lfcReps) ## variance over LFC inferential replicates
  print(sum(lfcVar < varThresh))
  
  if(is.null(tNodes)) {
    return(list(lfcCounts, varThresh, lfcVar))
  }
}

compLFCDiff <- function(tLFC, cLFC, tNode) {
  return(sqrt(mean((tLFC[tNode] - cLFC[tNode])^2)))
}

compMDist <- function(y, tNodes) {
  mIRV <- mcols(y)[["meanInfRV"]]
  tpm <- assays(y)[["abundance"]]
  sf <- DESeq2::estimateSizeFactorsForMatrix(tpm)
  scaledTpm <- t(t(tpm)/sf)
  mC1 <- rowMeans(scaledTpm[,colData(y)[["condition"]]==1]) ##mean condition 1
  mC2 <- rowMeans(scaledTpm[,colData(y)[["condition"]]==2]) ##mean condition 2
  dist=(mC1-mC2)^2/(mIRV)
  return(sqrt(mean(dist[tNodes])))
}

convAllTreeTxp <- function(clusFile, txpInd=28288) {
  trees <- read.tree(clusFile)
  trees <- lapply(trees, function(tr) {
    if(length(tr$tip)>2){
      tips <- as.numeric(tr$tip)
      tips <- tips[tips > txpInd]
      return(drop.tip(tr, tip=as.character(tips)))
    }
  })
  trees <- trees[!sapply(trees, is.null)]
}

compEffLenTxp <- function(txp_len, max_len=1000, mean = 200, sd = 25)
{
  library(truncnorm)
  if(txp_len > max_len)
    return(txp_len-mean)
  new_mean = mean(truncnorm::rtruncnorm(n=10000, b=txp_len, mean=mean, sd=sd))
  return(ifelse(txp_len-new_mean > 0,txp_len-new_mean, new_mean-txp_len))
}

compTPM <- function(counts, effLen) {
  if(nrow(counts) != length(effLen))
    stop("Lenghts not same as counts")
    vals <- counts/effLen
    tpms <- 1e6*t(t(vals)/(colSums(vals)))
    return(tpms)
}

##meanV is for adding variance+mean in denominator
compSPL <- function(y, i, pc = 5, mv=T, meanV=F) {
  if(!mv) {
    infReps <- assays(y)[grep("infRep", assayNames(y))]
    infReps <- abind::abind(as.list(infReps), along = 3)
    infMean <- apply(infReps, 1:2, mean)
    infVar <- apply(infReps, 1:2, var)
  }
  else {
    infMean <- assays(y)[["mean"]]
    infVar <- assays(y)[["variance"]]
  }
  spl <- abs(infVar[,i]-infMean[,i])/(infMean[,i]+pc)
  if(meanV)
    spl <- abs(infVar[,i]-infMean[,i])/(infVar[,i] + infMean[,i]+pc)
  # if(log)
  #   spl <- abs(log2(spl))
  return(spl)
}