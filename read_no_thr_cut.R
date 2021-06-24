#### Reading input of terminus run without no threshold cut
source("tree_helper_function.R")
# quantDir <- "../boot_gibbs/quant_output/human/swim_gencodev28"
# samples <- as.vector(outer(c(1:6), c("01","02"), function(x,y) paste("out",x,"sample",y,sep="_")))
# files <- file.path(quantDir, samples, "quant.sf")
# colData <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
# seSwim <- tximeta::tximeta(colData)
# save(seSwim, file = "environment/se_swim.RData")
load("environment/se_swim.RData")

### Reading trees
trees <- read.tree("../terminus/data/term_thr/cluster_nwk.txt")
tLengths <- sapply(trees, function(x) length(x$tip.label))
tree <- mergeTree(trees, se = seSwim)
save("environment/no_threshold/tree.RData")

print(length(tree$tip.label))

### Running Swish
y <- scaleInfReps(seSwim)
y <- labelKeep(y)
y <- runSwishtree(tree, y, type = "union")

### Creating final merged tree
modOb <- mergeLeaves(tree, y)
tree <- modOb[["tree"]]
y <- modOb[["ySwish"]]
save(y, file = "environment/no_threshold/swish_y.RData")
tree$node.label <- as.character(c(nrow(y)+1:tree$Nnode))
save(tree, file = "environment/no_threshold/tree.RData")
