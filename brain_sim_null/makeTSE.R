library(beaveR)
args <- commandArgs(trailingOnly = TRUE)

treeFile <- args[1]
quantDir <- args[2]
savePath <- args[3]
seed <- as.numeric(args[4])

if(!(file.exists(treeFile))) {
    stop("invalid tree file")
}
if(!(dir.exists(quantDir))) {
    stop("invalid quant dir")
}
if(!(dir.exists(savePath))) {
    stop("invalid save dor")
}
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep='_')))
quantFiles <- file.path(quantDir, samples, 'quant.sf')
coldata <- data.frame(files=quantFiles, names=samples, condition=factor(rep(1:2, each=6)))

tse <- buildTSE(treeTermFile = treeFile, coldata = coldata)
save(tse, file = file.path(savePath, paste("tse", "seed.RData", sep ="_")))
