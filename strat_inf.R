source("tree_helper_function.R")
source("getCandU.R")
library(treeclimbR)
library(ggplot2)
library(tictoc)

cuts <- c(4,6,8,12,16)
cuts <- c(2,3,5)

## Whole transcriptome
load("environment/whole_txp/swish_y.RData")
load("environment/whole_txp/tree.RData")
load("environment/whole_txp/ySwishAllNodes.RData")

infRVs <- mcols(yAll)[,"meanInfRV"]
cSwishInf <- vector("list", length(cuts))
bSwishInf <- vector("list", length(cuts))
for(i in seq_along(cuts))
{
    infCuts <- cut_number(infRVs, cuts[i])
    indsIV <- lapply(unique(infCuts), function(cut) which(infCuts == cut))
    names(indsIV) <- unique(infCuts)

    swishL <- createSwishOb(indsIV, y, cores=2)
    print(swishL)
    gc()
    cSwishInf[[i]] <- getCandU(tree,
                          score_data = swishL[["df"]],
                          node_column = "node",
                          p_column = "pvalue",
                          sign_column = "log2FC",
                          cores = 1)
    save(cSwishInf, file = "environment/whole_txp/c_swish_infRV_diff_rem.RData")
    bSwishInf[[i]] <- mclapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = cSwishInf[[i]]$candidate_list,
                                                    score_data = swishL[["df"]], node_column = "node",
                                                    p_column = "pvalue", sign_column = "log2FC", limit_rej=x), mc.cores = 2)
    save(bSwishInf, file = "environment/whole_txp/b_swish_infRV_diff_rem.RData")
    gc()
}
# 
### Only on tree txps
# load("environment/innNodes/ySwish.RData")
# load("environment/ygroup.RData")
# load("environment/tree_group.RData")
# 
# infRVs <- mcols(yAll)[,"meanInfRV"]
# cSwishInf <- vector("list", length(cuts))
# bSwishInf <- vector("list", length(cuts))
# for(i in seq_along(cuts))
# {
#     infCuts <- cut_number(infRVs, cuts[i])
#     indsIV <- lapply(unique(infCuts), function(cut) which(infCuts == cut))
#     names(indsIV) <- unique(infCuts)
#     
#     swishL <- createSwishOb(indsIV, y, cores=2)
#     gc()
#     cSwishInf[[i]] <- getCandU(tree,
#                                score_data = swishL[["df"]],
#                                node_column = "node",
#                                p_column = "pvalue",
#                                sign_column = "log2FC",
#                                cores = 1)
#     gc()
#     save(cSwishInf, file = "environment/innNodes/c_swish_infRV_diff.RData")
#     bSwishInf[[i]] <- mclapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = cSwishInf[[i]]$candidate_list,
#                                                                       score_data = swishL[["df"]], node_column = "node",
#                                                                       p_column = "pvalue", sign_column = "log2FC", limit_rej=x), mc.cores = 2)
#     save(bSwishInf, file = "environment/innNodes/b_swish_infRV_diff.RData")
#     gc()
# }