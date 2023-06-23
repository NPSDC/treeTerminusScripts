treeTermPath <- '/fs/cbcb-lab/rob/students/noor/Uncertainity/treeTerminusScripts'
suppressPackageStartupMessages(source(file.path(treeTermPath, "tree_helper_function.R")))
suppressPackageStartupMessages(source(file.path(treeTermPath, "tree_term_climb.R")))
suppressPackageStartupMessages(library(beaveR))

args <- commandArgs(trailingOnly = TRUE)

seed <- as.nuemric(args[1])
mainDir <- args[2]
clustFile <- file.path(mainDir, "terminus", paste("seed",seed,sep="="), "no_threshold0", "cluster_nwk.txt")
quantDir <- file.path(mainDir, "out_sal", paste("seed",seed,sep="="))
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep='_')))
quantFiles <- file.path(quantDir, samples, 'quant.sf')
coldata <- data.frame(files=quantFiles, names=samples, condition=factor(rep(1:2, each=6)))

tse <- buildTSE(treeTermFile = clustFile, coldata = coldata)
tree <- rowTree(tse)
nl <- length(tree$tip)

y <- tse[1:nl,]
y <- fishpond::scaleInfReps(y)
hist(mcols(pvalue))

yAll <- computeSizeFactors(tse)
yAll <- scaleInfReps(yAll)
yAll <- labelKeep(yAll)
set.seed(1)
yAll <- swish(yAll, x = "condition")
