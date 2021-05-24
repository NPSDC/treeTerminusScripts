library(BOUTH)
source("tree_helper_function.R")
## Reading data
# quantDir <- "../boot_gibbs/quant_output/human/swim_gencodev28"
# samples <- as.vector(outer(c(1:6), c("01","02"), function(x,y) paste("out",x,"sample",y,sep="_")))
# files <- file.path(quantDir, samples, "quant.sf")
# colData <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
# seSwim <- tximeta::tximeta(colData)
# save(seSwim, file = "environment/se_swim.RData")
load("environment/se_swim.RData")

### Reading trees
trees <- read.tree("../terminus/data/term/cluster_nwk.txt")
tLengths <- sapply(trees, function(x) length(x$tip.label))
tree <- mergeTree(trees)

### Running Swish
# y <- scaleInfReps(seSwim)
# y <- labelKeep(y)
# y <- runSwishtree(tree, y, type = "union")
# save(y, file = "environment/swish_y.RData")
load("environment/swish_y.RData")

### Creating final merged tree
# length(setdiff(as.numeric(tree$tip.label), match(rownames(y), rownames(seSwim)))) == 0
# tree <- mergeLeaves(tree, y, seSwim)
# save(tree, file = "environment/tree.RData")
load("environment/tree.RData")

### Creating data frame needed by BOUTH
treeDf <- getTreeDf(tree)
save(treeDf, file="environment/treeDf.RData")

### Other checks needed, find the tree with max depth and see if it matches with our depth of tree
## Depth given by tLength-1
maxDepth = max(tLengths) - 1 #91
## Depth of new tree is 93 which is what is expected


### BOUTH Code
pvalues <- mcols(y)[["pvalue"]][as.numeric(rownames(treeDf))]
plot(hist(pvalues))

bouth_test <- list()
bouth_test[['fdr_0.1']] <- bouth(anno.table = treeDf, pvalue.leaves = pvalues, na.symbol = "Undefined", far = 0.1, is.weighted = TRUE)
bouth_test[['fdr_0.05']] <- bouth(anno.table = treeDf, pvalue.leaves = pvalues, na.symbol = "Undefined", far = 0.05, is.weighted = TRUE)
save(bouth_test, file="environment/bouth_test_all.RData")
### Keeping the full trees that includes all txps seem to be making root as driver node


### Limiting the set
tree <- mergeTree(trees)
y <- scaleInfReps(seSwim)
y <- labelKeep(y)
y <- runSwishtree(tree, y, type = "tree")
tree <- convTree(tree, y, seSwim)
treeDfGroup <- getTreeDf(tree)
save(tree, file="environment/tree_group.RData")
save(treeDfGroup, file="environment/treeDfGroup.RData")
save(y, file="environment/ygroup.RData")

pvalues <- mcols(y)[["pvalue"]][as.numeric(rownames(treeDfGroup))]
plot(hist(pvalues))
bouth_test <- list()
bouth_test[['fdr_0.1']] <- bouth(anno.table = treeDfGroup, pvalue.leaves = pvalues, na.symbol = "Undefined", far = 0.1, is.weighted = TRUE)
bouth_test[['fdr_0.05']] <- bouth(anno.table = treeDfGroup, pvalue.leaves = pvalues, na.symbol = "Undefined", far = 0.05, is.weighted = TRUE)
save(bouth_test, file = "environment/bouth_test_group.RData")

### DESeq2 on simulated
txi <- tximeta::tximeta(colData, type = "salmon", txOut = TRUE, countsFromAbundance="lengthScaledTPM")
ddsSim <- DESeqDataSet(txi[rownames(y),], design=~condition)
ddsSim <- DESeq(ddsSim)
res <- results(ddsSim, name="condition_2_vs_1")
res$pvalue[is.na(res$pvalue)] <- 1
res$padj[is.na(res$padj)] <- 1
save(res, file="environment/resDeseq2_sim.RData")
pvalues <- res$pvalue[match(rownames(res),rownames(y)[as.numeric(tree$tip.label)])]
bouth_sim <- bouth(anno.table = treeDfGroup, pvalue.leaves = pvalues, na.symbol = "Undefined", far = 0.05, is.weighted = TRUE)
save(bouth_sim, file="environment/bouth_sim.RData")
tseSim <- createTreeSummarizedObject(y, treeDfGroup, tree)

### Issues with true counts
# 1. Missing txps
# 2. No normalization
### Reading true counts data
dir <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/swim_data"
samples <- as.vector(outer(c(1:6), c("1","2"), function(x,y) paste(paste("out",x, sep = "_"), paste(paste("out",y, sep=""), "counts.tsv", sep = "_"), sep="/")))
files <- file.path(dir, samples)
dd <- readTrueCounts(files, samples)
mTxps <- setdiff(rownames(y),rownames(dd)) ## txps that are missing
print(length(mTxps))
mDf <- matrix(0,nrow=length(mTxps),ncol=ncol(dd), dimnames=list(mTxps, colnames(dd)))
dd <- rbind(dd, mDf)
ddsTrue <- DESeqDataSetFromMatrix(countData = dd[rownames(y),], colData = colData, design=~condition)
ddsTrue <- DESeq(ddsTrue)
resTrue <- results(ddsTrue)
resTrue$pvalue[is.na(resTrue$pvalue)] <- 1
resTrue$padj[is.na(resTrue$padj)] <- 1
save(resTrue, file="environment/resDeseq2_true.RData")
pvalues <- resTrue$pvalue[match(rownames(resTrue),rownames(y)[as.numeric(tree$tip.label)])]
bouth_true <- bouth(anno.table = treeDfGroup, pvalue.leaves = pvalues, na.symbol = "Undefined", far = 0.05, is.weighted = TRUE)
save(bouth_true, file="environment/bouth_true.RData")