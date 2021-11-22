suppressPackageStartupMessages(source("tree_helper_function.R"))
library(phangorn)
load("../mikelove-swimdown-216a1dd/simulate/data/simulate.rda")
load("environment/ygroup.RData")
load("environment/tree_group.RData")
load("environment/parTClimbR/bSwish.RData")
load("environment/parTClimbR/bSwishInfRV.RData")
load("environment/parTClimbR/bDESeq2.RData")
load("environment/parTClimbR/b_deseq2_infRV.RData")

missingTxps <- setdiff(rownames(y), rownames(sim.counts.mat))
sim.counts.mat <- rbind(sim.counts.mat, matrix(0, nrow = length(missingTxps), ncol = ncol(sim.counts.mat),
                                               dimnames = list(missingTxps, colnames(sim.counts.mat))))
sim.counts.mat <- sim.counts.mat[rownames(y),] ###Arranging in the manner of y

innNodes <- nrow(y)+1:tree$Nnode

aggCountsNodes <- computeAggNodes(tree, c(1:nrow(y),innNodes), sim.counts.mat)
fold_changes <- rbind(fold_changes, matrix(0, nrow = length(missingTxps), ncol = ncol(fold_changes),
                                           dimnames = list(missingTxps, colnames(fold_changes))))
logFCNodes <- ifelse(rowSums(aggCountsNodes)==0, 0, log2(aggCountsNodes[,2]/aggCountsNodes[,1]))
fold_changes <- fold_changes[rownames(y),]

trueCountNodes <- which(abs(logFCNodes) >= 0.066) ##Making sure root is non differentiated
trueDriverNodes <-  findDriverGround(trueCountNodes, tree, logFCNodes)[[1]] ### All underlying children have same sign
trueDriverNodes <- sort(trueDriverNodes)

detNodes <- list()
detNodes[["OnlySwishL"]] <- lapply(c(0.01,0.05,0.1), function(x) which(mcols(y)[,"qvalue"] <= x)) ##Leaves
detNodes[["tcrSw"]] <- sapply(bSwish, function(sw) sw$output[sw$output$signal.node, "node"]) ##Treeclimbr Swish on all nodes
detNodes[["tcrSwInf"]] <- sapply(bSwishInf, function(sw) sw$output[sw$output$signal.node, "node"]) ##Treeclimbr Swish on nodes stratified by 

fpNodes <- setdiff(detNodes[[2]][[3]],trueDriverNodes)
fpNodesLeaves <- fpNodes[fpNodes <= nrow(y)]
fpInnNodes <- fpNodes[fpNodes > nrow(y)]
print(length(fpNodesLeaves))
print(length(fpInnNodes))
l <- findDistance(tree, fpNodes, trueDriverNodes)
print(table(l))
print(sum(is.na(l)))
### Many txps have same ratio or both are 0

load("../mikelove-swimdown-216a1dd/simulate/data/simulate.rda")
load("environment/whole_txp/swish_y.RData")
load("environment/whole_txp/tree.RData")
missingTxps <- setdiff(rownames(y), rownames(sim.counts.mat))
sim.counts.mat <- rbind(sim.counts.mat, matrix(0, nrow = length(missingTxps), ncol = ncol(sim.counts.mat),
                                               dimnames = list(missingTxps, colnames(sim.counts.mat))))
sim.counts.mat <- sim.counts.mat[rownames(y),] ###Arranging in the manner of y

innNodes <- nrow(y)+1:tree$Nnode

aggCountsNodes <- computeAggNodes(tree, c(1:nrow(y),innNodes), sim.counts.mat)
fold_changes <- rbind(fold_changes, matrix(0, nrow = length(missingTxps), ncol = ncol(fold_changes),
                                           dimnames = list(missingTxps, colnames(fold_changes))))
logFCNodes <- ifelse(rowSums(aggCountsNodes)==0, 0, log2(aggCountsNodes[,2]/aggCountsNodes[,1]))
fold_changes <- fold_changes[rownames(y),]

trueCountNodes <- which(abs(logFCNodes) >= 0.066) ##Making sure root is non differentiated
trueDriverNodes <-  findDriverGround(trueCountNodes, tree, logFCNodes)[[1]] ### All underlying children have same sign
trueDriverNodes <- sort(trueDriverNodes)

load("environment/whole_txp/bSwish.RData")
load("environment/whole_txp/bSwishInfRV.RData") 
detNodes <- list()
detNodes[["OnlySwishL"]] <- lapply(c(0.01,0.05,0.1), function(x) which(mcols(y)[,"qvalue"] <= x)) ##Leaves
detNodes[["tcrSw"]] <- sapply(bSwish, function(sw) sw$output[sw$output$signal.node, "node"]) ##Treeclimbr Swish on all nodes
detNodes[["tcrSwInf"]] <- sapply(bSwishInf, function(sw) sw$output[sw$output$signal.node, "node"]) ##Treeclimbr Swish on nodes stratified by 
