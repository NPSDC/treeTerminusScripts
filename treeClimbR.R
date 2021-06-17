#### Loading data
source("tree_helper_function.R")
source("getCandU.R")
load("environment/ygroup.RData")
load("environment/tree_group.RData")
library(treeclimbR)
library(ggplot2)
library(tictoc)
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
load("environment/innNodes/res.RData")
# 
# #### Running treeClimbR on pvalues obtained using DESeq2
# tic()
# cDeseq2 <- getCandU(tree,
#                       score_data = res,
#                       node_column = "node",
#                       p_column = "pvalue",
#                       sign_column = "log2FoldChange", cores = 8)
# toc()
# save(cDeseq2, file = "environment/parTClimbR/cDeseq2.RData")
# bDeseq2 <- lapply(c(0.01, 0.05, 0.1), function(x) evalCand(tree = tree, levels = cDeseq2$candidate_list,
#                       score_data = data.frame(res), node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FoldChange", limit_rej=x))
# save(bDeseq2, file = "environment/parTClimbR/bDESeq2.RData")

#### Preparing for Swish
# yAll <- prepSwish(tree, y)
# yAll <- scaleInfReps(yAll, lengthCorrect = T)
# yAll <- labelKeep(yAll)
# mcols(yAll)$keep=T
# yAll <- swish(yAll, x = "condition")
# yAll <- computeInfR
# save(yAll, file = "environment/innNodes/ySwish.RData")
load("environment/innNodes/ySwish.RData")

swishRes <- data.frame(mcols(yAll))
swishRes$node <- seq(1:nrow(yAll))
swishRes <- swishRes[,-c(2)] ### Removing keep

#### Running treeClimbR on pvalues obtained using Swish
# tic()
# cSwish <- getCandU(tree,
#                  score_data = swishRes,
#                  node_column = "node",
#                  p_column = "pvalue",
#                  sign_column = "log2FC",
#                  cores = 6)
# toc()
# save(cSwish, file = "environment/parTClimbR/cSwish.RData")
# bSwish <- mclapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = cSwish$candidate_list,
#                       score_data = swishRes, node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FC", limit_rej=x), mc.cores = 3)
# save(bSwish, file = "environment/parTClimbR/bSwish.RData")
# gc()

#### Stratifying by inferential variance
infRVs <- mcols(yAll)[,"meanInfRV"]
infCuts <- cut_number(infRVs, 6)
indsIV <- lapply(unique(infCuts), function(cut) which(infCuts == cut))
names(indsIV) <- unique(infCuts)

# de <- createDESeq2Ob(indsIV, aggSimCountsNodes, y, cores=6)
# save(de, file = "environment/innNodes/dd_res_infRV.RData")
# load("environment/innNodes/dd_res_infRV.RData")
# tic()
# cDeseq2Inf <- getCandU(tree, score_data = de$res,
#                      node_column = "node",
#                      p_column = "pvalue",
#                      sign_column = "log2FoldChange",
#                      cores = 6)
# toc()
# save(cDeseq2Inf, file = "environment/innNodes/c_deseq2_infRV.RData")
# bDeseq2Inf <- mclapply(c(0.01, 0.05, 0.1), function(x) evalCand(tree = tree, levels = cDeseq2Inf$candidate_list,
#                       score_data = data.frame(de[["res"]]), node_column = "node",
#                       p_column = "pvalue", sign_column = "log2FoldChange", limit_rej=x), mc.cores = 3)
# save(bDeseq2Inf, file = "environment/parTClimbR/b_deseq2_infRV.RData")
# gc()

# swishL <- createSwishOb(indsIV, y, cores=6)
# save(swishL, file = "environment/innNodes/swish_infRV.RData")
load("environment/innNodes/swish_infRV.RData")
tic()
cSwishInf <- getCandU(tree,
                     score_data = swishL[["df"]],
                     node_column = "node",
                     p_column = "pvalue",
                     sign_column = "log2FC",
                     cores = 2)
toc()
save(cSwishInf, file = "environment/parTClimbR/c_swish_infRV.RData")
bSwishInf <- mclapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = cSwishInf$candidate_list,
                      score_data = swishL[["df"]], node_column = "node",
                      p_column = "pvalue", sign_column = "log2FC", limit_rej=x), mc.cores = 3)
save(bSwishInf, file = "environment/parTClimbR/bSwishInfRV.RData")
gc()