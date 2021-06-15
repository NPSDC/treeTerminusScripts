#### Loading data
source("tree_helper_function.R")
load("environment/ygroup.RData")
load("environment/tree_group.RData")
library(treeclimbR)
library(ggplot2)

#### Creating aggregated counts
tree$node.label <- as.character(c(nrow(y)+1:tree$Nnode))
aggSimCountsNodes <- computeAggNodesU(tree, c(1:(nrow(y)+tree$Nnode)), assays(y)[["counts"]], NULL)

### Preparing for DESeq2
# dds <- DESeqDataSetFromMatrix(countData = round(aggSimCountsNodes),
#                               colData = colData(y),
#                               design = ~ condition)
# dds <- DESeq(dds, minReplicatesForReplace=Inf)
# save(dds, file = "environment/innNodes/dds.RData")
# 
# res <- results(dds, name="condition_2_vs_1", independentFilter = F, cooksCutoff = Inf)
# res$node <- seq(1:nrow(aggSimCountsNodes))
# save(res, file = "environment/innNodes/res.RData")
# load("environment/innNodes/res.RData")
# 
# #### Running treeClimbR on pvalues obtained using DESeq2
# candDeseq2 <- getCand(tree,
#                      score_data = res,
#                      node_column = "node",
#                      p_column = "pvalue",
#                      sign_column = "log2FoldChange")
# save(candDeseq2, file = "environment/innNodes/candDeseq2.RData")
# bestDeseq2 <- lapply(c(0.01, 0.05, 0.1), function(x) evalCand(tree = tree, levels = candDeseq2$candidate_list,
#                       score_data = data.frame(res), node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FoldChange", limit_rej=x))
# save(bestDeseq2, file = "environment/innNodes/bestDESeq2.RData")
# 
# #### Preparing for Swish
# yAll <- prepSwish(tree, y)
# yAll <- scaleInfReps(yAll, lengthCorrect = T)
# yAll <- labelKeep(yAll)
# mcols(yAll)$keep=T
# yAll <- swish(yAll, x = "condition")
# yAll <- computeInfR
# save(yAll, file = "environment/innNodes/ySwish.RData")
# load("environment/innNodes/ySwish.RData")

# swishRes <- data.frame(mcols(yAll))
# swishRes$node <- seq(1:nrow(yAll))
# swishRes <- swishRes[,-c(2)] ### Removing keep
# 
# #### Running treeClimbR on pvalues obtained using Swish
# candSwish <- getCand(tree,
#                      score_data = swishRes,
#                      node_column = "node",
#                      p_column = "pvalue",
#                      sign_column = "log2FC")
# save(candSwish, file = "environment/innNodes/candSwish.RData")
# bestSwish <- lapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = candSwish$candidate_list,
#                       score_data = swishRes, node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FC", limit_rej=x))
# save(bestSwish, file = "environment/innNodes/bestSwish.RData")


#### Stratifying by inferential variance
# infRVs <- mcols(yAll)[,"meanInfRV"]
# infCuts <- cut_number(infRVs, 6)
# indsIV <- lapply(unique(infCuts), function(cut) which(infCuts == cut))
# names(indsIV) <- unique(infCuts)

# de <- createDESeq2Ob(indsIV, aggSimCountsNodes, y, cores=6)
# save(de, file = "environment/innNodes/dd_res_infRV.RData")
# load("environment/innNodes/dd_res_infRV.RData")
# candDeInf <- getCand(tree, score_data = de$res,
#                      node_column = "node",
#                      p_column = "pvalue",
#                      sign_column = "log2FoldChange")
# save(candDeInf, file = "environment/innNodes/cand_deseq2_infRV.RData")
# bestDeInf <- sapply(c(0.01, 0.05, 0.1), function(x) evalCand(tree = tree, levels = candDeInf$candidate_list,
#                       score_data = data.frame(de[["res"]]), node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FoldChange", limit_rej=x))
# save(bestDeInf, file = "environment/innNodes/best_deseq2_infRV.RData")

# swishL <- createSwishOb(indsIV, y, cores=6)
# save(swishL, file = "environment/innNodes/swish_infRV.RData")
# load("environment/innNodes/swish_infRV.RData")
# candSwInf <- getCand(tree,
#                      score_data = swishL[["df"]],
#                      node_column = "node",
#                      p_column = "pvalue",
#                      sign_column = "log2FC")
# save(candSwInf, file = "environment/innNodes/cand_swish_infRV.RData")
# bestSwInf <- lapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = candSwInf$candidate_list,
#                       score_data = swishL[["df"]], node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FC", limit_rej=x))
# save(bestSwInf, file = "environment/innNodes/bestSwInf.RData")
