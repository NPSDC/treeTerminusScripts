#### Creating and modifying the data for whole transcriptome level FDR
source("tree_helper_function.R")
load("environment/se_swim.RData")

### Reading trees
trees <- read.tree("../terminus/data/term/cluster_nwk.txt")
tLengths <- sapply(trees, function(x) length(x$tip.label))
tree <- mergeTree(trees, se = seSwim)

### Running Swish
y <- scaleInfReps(seSwim)
y <- labelKeep(y)
y <- runSwishtree(tree, y, type = "union")
save(y, file = "environment/whole_txp/swish_y.RData")
# load("environment/swish_y.RData")

### Creating final merged tree
modOb <- mergeLeaves(tree, y)
tree <- modOb[["tree"]]
y <- modOb[["ySwish"]]
save(tree, file = "environment/whole_txp/tree.RData")
save(y, file = "environment/whole_txp/swish_y.RData")
# load("environment/tree.RData")

