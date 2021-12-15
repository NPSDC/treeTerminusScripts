library(ape)
library(phangorn)
merged_group_file <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/brain_sim/out_term/merged_groups_length.txt"
mgroupFiles <- paste("/fs/cbcb-lab/rob/students/noor/Uncertainity/brain_sim/out_term/", sort(sapply(as.vector(outer(c(1:6),c(1:2), function(x,y) paste(x,y, sep="_"))), function(x) paste(x,"mgroup_nwk.txt", sep="/"))), sep="")

trees <- lapply(mgroupFiles, function(f) read.tree(f))
inpPhyl <- rep("a", length(trees)*length(trees[[1]])+length(trees[[1]]))
j=1
for(i in seq_along(trees[[1]])) {
   inpPhyl[j] <- paste(">Group", rownames(rDf)[i], sep = " ")
   j=j+1
   for(k in seq_along(trees)) {
      # inpPhyl[j] <- gsub(";", ":0;", write.tree(trees[[k]][[i]]))
       inpPhyl[j] <- write.tree(trees[[k]][[i]])
       j=j+1
   }
}
write(inpPhyl, file = "../brain_sim/out_term/trees_phyl_inp.nwk")

rDf <- read.delim(merged_group_file, sep="\t", header=F, row.names=1)
dfScore <- data.frame(row.names = rownames(rDf)) #number of trees per sample corresponding to the group
for(i in seq_along(mgroupFiles)) {
    scores <- sapply(strsplit(rDf[,2*i], ",", fixed = T), length)
    dfScore <- cbind(dfScore, scores)
}
colnames(dfScore)[1:ncol(dfScore)] <- paste("nTrees", sapply(strsplit(mgroupFiles, "/", fixed=T), function(x) x[1]), sep = "_")

dfNtxps <- data.frame(row.names = rownames(rDf)) ##number of txps covered in a sample corresponding to the group
for(i in seq_along(mgroupFiles)) {
    scores <- sapply(strsplit(rDf[,2*i], ",", fixed = T), function(x) {
        tGs <- unlist(x)
        length(unlist(strsplit(tGs,"_", fixed=T)))
        })
    dfNtxps <- cbind(dfNtxps, scores)
}
colnames(dfNtxps)[1:ncol(dfNtxps)] <- paste("nTrees", sapply(strsplit(mgroupFiles, "/", fixed=T), function(x) x[1]), sep = "_")

## type - r remove all leaves 
## type - NULL keep it
getSepTrees <- function(tree, score, ntxps, type = NULL) {
    nL <- length(tree$tip)
    if(score == 0)
        return(NULL)
    if(score==1 & ntxps == nL) ##tree already correct
        return(write.tree(tree))
    if(!is.null(type)) {  
        if(ntxps == 2 & type == "r") {### Taking care of the 2txp case
            if(score!=1)
                stop("Score should be 1 if ntxps = 2")
            return(write.tree(tree))
        }
    }
    
    child <- Descendants(tree, nL+1, type="children")
    if(all(child <= nL)) { ### atleast 1 internal node since it does not meet nL==ntxps
        print("all leaves")
        return(-1)
    }
    
    remChild <- child[child > nL] ## Nodes left after removing txps - these are the ones not added
    nwks <- sapply(remChild, function(node) write.tree(extract.clade(tree, node)))
    nwksL <- sapply(strsplit(nwks, ",", fixed=T), length)
    if(sum(nwksL) != ntxps)
        stop("Something wrong in length calculation")
    if(any(nwksL == 2) & score > 1) {
        nwks <- gsub(";", "", nwks) 
        nwks <- paste(nwks,collapse=",")
        nwks <- paste("(", nwks, ");", sep = "")
    }
    return(nwks)
}

### 2 options - remove all doublets and keep all
getSepTrees <- function(tree, score, ntxps, groupLength) {
    
    
}

getNTrees <- function(trees, score, ntxps, groupLength) {
    txps <- table(unlist(lapply(trees, function(tree) tree$tip)))
}
### Write old
groupLength <- sapply(strsplit(rownames(rDf), split="_", fixed=T), lengtxh)


vals <- rep("a", nrow(rDf) + sum(dfScore))
k=1
for(i in seq(nrow(rDf)))
{
    vals[k] <- paste(">Group", rownames(rDf)[i], sep = " ")
    k=k+1
    for(j in seq_along(trees)) {
        treeNwks <- getSepTrees(trees[[j]][[i]], dfScore[i,j], dfNtxps[i,j])
        lTree <- length(treeNwks)
        if(!is.null(treeNwks)) {
            # if(lTree != dfScore[i,j]){
            #     stop("not enough trees")
            # }
                
            vals[k:(k+lTree-1)] <- treeNwks
            k <- k + lTree
        }
    }   
}

vals <- vals[vals!="a"]
write(vals, file = "trees_rem.nwk")

vals <- rep("a", nrow(rDf) + sum(dfScore))
k=1
for(i in seq(nrow(rDf)))
{
    vals[k] <- paste(">Group", rownames(rDf)[i], sep = " ")
    k=k+1
    for(j in seq_along(trees)) {
        treeNwks <- getSepTrees(trees[[j]][[i]], dfScore[i,j], dfNtxps[i,j], "r")
        lTree <- length(treeNwks)
        if(!is.null(treeNwks)) {
            # if(lTree != dfScore[i,j]){
            #     stop("not enough trees")
            # }
            
            vals[k:(k+lTree-1)] <- treeNwks
            k <- k + lTree
        }
    }   
}

vals <- vals[vals!="a"]
write(vals, file = "trees_not_rem.nwk")