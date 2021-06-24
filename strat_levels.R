source("tree_helper_function.R")
source("getCandU.R")
library(treeclimbR)
load("environment/whole_txp/tree.RData")
load("environment/whole_txp/swish_y.RData")

levels <- node.depth(tree, 2)
levels <- sapply(levels, function(x) {
    if(x <= 4)
        as.character(x)
    else
        "5"
})
levInds <- lapply(unique(levels), function(lev) which(levels == lev))
swishL <- createSwishOb(levInds, y, cores=4)
gc()
df <- swishL[["df"]]
rm(swishL)
gc()
cSwishLev <- getCandU(tree,
                           score_data = df,
                           node_column = "node",
                           p_column = "pvalue",
                           sign_column = "log2FC",
                           cores = 2)
gc()
save(cSwishLev, file = "environment/whole_txp/c_swish_levels.RData")
bSwishLev <- mclapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = cSwishLev$candidate_list,
                                                                  score_data = df, node_column = "node",
                                                                  p_column = "pvalue", sign_column = "log2FC", limit_rej=x), mc.cores = 2)
save(bSwishLev, file = "environment/whole_txp/b_swish_levels.RData")
gc()

### Only tree txps
load("environment/tree_group.RData")
load("environment/ygroup.RData")

levels <- node.depth(tree, 2)
levels <- sapply(levels, function(x) {
    if(x <= 4)
        as.character(x)
    else
        "5"
})
levInds <- lapply(unique(levels), function(lev) which(levels == lev))
swishL <- createSwishOb(levInds, y, cores=4)
gc()
df <- swishL[["df"]]
rm(swishL)
gc()
cSwishLev <- getCandU(tree,
                      score_data = df,
                      node_column = "node",
                      p_column = "pvalue",
                      sign_column = "log2FC",
                      cores = 2)
gc()
save(cSwishLev, file = "environment/innNodes/c_swish_levels.RData")
bSwishLev <- mclapply(c(0.01,0.05,0.1), function(x) evalCand(tree = tree, levels = cSwishLev$candidate_list,
                                                             score_data = df, node_column = "node",
                                                             p_column = "pvalue", sign_column = "log2FC", limit_rej=x), mc.cores = 2)
save(bSwishLev, file = "environment/innNodes/b_swish_levels.RData")
gc()