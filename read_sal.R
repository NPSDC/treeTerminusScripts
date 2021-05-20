library(tximeta)
library(tximport)
library(ape)
library(phangorn)
library(fishpond)
library(TreeSummarizedExperiment)

## Rading data
quantDir <- "../boot_gibbs/quant_output/human/swim_gencodev28"
samples <- as.vector(outer(c(1:6), c("01","02"), function(x,y) paste("out",x,"sample",y,sep="_")))
files <- file.path(quantDir, samples, "quant.sf")
colData <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
seSwim <- tximeta::tximeta(colData)
save(seSwim, file = "environment/se_swim.RData")

### Reading trees
load("environment/se_swim.RData")
trees <- read.tree("../terminus/data/term/cluster_nwk.txt")
tLengths <- sapply(trees, function(x) length(x$tip.label))
tree <- mergeTree(trees)

### Running Swish
y <- scaleInfReps(seSwim)
y <- labelKeep(y)
y <- runSwishtree(tree, y)
save(y, file = "environment/swish_y.RData")
load("environment/swish_y.RData")

### Creating final merged tree
length(setdiff(as.numeric(tree$tip.label), match(rownames(y), rownames(seSwim)))) == 0
tree <- mergeLeaves(tree, y, seSwim)
save(tree, file = "environment/tree.RData")
load("environment/tree.RData")

### Creating data frame needed by BOUTH
treeDf <- getTreeDf(tree)

### Bouth Code

### TreeSummarizedExperiment
### Other checks needed, find the tree with max depth and see if it matches with our depth of tree






