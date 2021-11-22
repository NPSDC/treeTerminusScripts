source("tree_helper_function.R")
source("tree_term_climb.R")
set.seed(10)
l <- 10
tree <- rtree(l)
n <- 10
means1 <- c(10,10,5,6,20,20,10,10,20,5)
means2 <- c(10.8,10.9,6,8,20,20,10,10,20,10)

set.seed(5)
counts <- matrix(0, nrow = l, ncol = 2*n)
for(i in seq(l))
    counts[i,] <- rnorm(2*n, rep(c(means1[i], means2[i]), each = n))

aggCounts <- computeAggNodesU(tree, 1:(length(tree$tip) + tree$Nnode), counts, group_inds = NULL)
pvalues <- apply(aggCounts, 1, function(x) t.test(x[1:n], x[n+1:n])[["p.value"]])

lfcs <- log2(rowMeans(aggCounts[,1:n])) - log2(rowMeans(aggCounts[,n+1:n]))
signs <- rep(-1,nrow(aggCounts))
mInfRV <- rnorm(nrow(aggCounts),20)


pvalues <- c(0.4, 0.3, 0.1, 0.5, 0.3, 0.4, 0.2, 0.3, 0.1, 0.02,  0.4, 0.3, 0.01, 0.1, 0.04, 0.01, 0.5, 0.3, 0.1)
#pvalues[c(18,19)]=c(0.1,0.2)
nodeSig <- rep(T, 19) ### node signficant or not
nodeSig[pvalues > 0.05]=F
#checkGoUp(13, 1,1, signs, mInfRV, 0.5, nodeSig)

signs[c(6)] = 1
cands <- computeReps(tree, aggCounts, "cond", pCutOff = 0.05, 0.5)

load("environment/se_swim.RData")
load("environment/whole_txp/tree.RData")
load("environment/whole_txp/ySwishAllNodes.RData")
load("environment/whole_txp/swish_y.RData")
load("environment/treeTerm/tBefore3.RData")
trees <- read.tree("../terminus/data/term/cluster_nwk.txt")
tLengths <- sapply(trees, function(x) length(x$tip.label))
which(tLengths==20)
tree <- trees[[7419]]
plot(tree)
y <- seSwim[as.numeric(tree$tip.label)+1,]

y <- scaleInfReps(y)
y <- labelKeep(y)
y <- runSwishtree(tree, y, type = "union")
resDf <- tBefore3[[3]][["resDf"]]
yAll <- prepSwish(tree, y)
yAll <- scaleInfReps(yAll, lengthCorrect = T)
yAll <- labelKeep(yAll)
mcols(yAll)$keep=T
yAll <- swish(yAll, x = "condition")
yAll <- computeInfRV(yAll)
d <- runTreeTermAlpha(tree, yAll, "condition", 100, mcols(yAll)[["pvalue"]], pCutOff = 0.05)
termSwimThresh <- mean(c(-3.8223631318627858,-3.6433994917165045,-3.6568587953449123,-3.7615266573162742,-3.74218233312281,
                         -3.7825084537263063,-3.681508024209264,-3.6888354232562355,-3.811711665323895,-3.7328311655820956,-3.6988339545709255,-3.892522972438748))
tAfter <- runTreeTermAlphas(tree, yAll, "condition", termSwimThresh, pCutOff = 0.05, ihwType = c("a"), alphas = c(0.01, 0.05, 0.10), cores = 1)

signs <- computeSign(yAll, "condition")
tree2 <- trees[[51]]
tree2$tip.label <- rownames(seSwim)[as.numeric(tree2$tip.label)+1]
tInds <- sort(c(553303,unlist(Descendants(tree,55330,"all"))))
yA <- yAll[tInds,]
yR <- resDf[tInds,]

nodeSig <- rep(T, nrow(yR))
nodeSig[yR[["adj_pvalue"]] > 0.05] = F
mInfRV <- mcols(yAll)[["meanInfRV"]][tInds]
signsR <- signs[tInds]

p <- 9
child <- unlist(Descendants(tree2,p,"children"))
desc <- unlist(Descendants(tree2,p,"all"))
checkGoUp(p, desc, child, signsR, mInfRV, termSwimThresh, nodeSig)

dN <- tBefore3[[3]][["negNodeO"]]
dNI <- dN[dN > nrow(y)]
which(sapply(dNI,function(node) any(Ancestors(tree, node) %in% dN)))

d <- runTreeTermAlpha(tree2, yA, "condition", 100, yR[["adj_pvalue"]], pCutOff = 0.05)