#quants <- read.delim("ERR188297/quant.sf", string=FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
libSize <- as.numeric(args[1])
nsampgenes <- as.numeric(args[2])
cortex_tpm <- read.delim("/fs/cbcb-lab/rob/students/noor/Uncertainity/treeTerminusScripts/brain_simulation/frontal_cortex.tsv.gz", skip=2)
seed <- as.numeric(args[3])
saveDir <- args[4]
dataDir <- args[5]

# average over the 10 samples
quants <- data.frame(row.names = cortex_tpm$transcript_id,
                     TPM = rowMeans(cortex_tpm[,-c(1,2)]))
# rounding difference
quants$TPM <- quants$TPM * 1e6 / sum(quants$TPM)

suppressPackageStartupMessages(library(GenomicFeatures))
ann_dir <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/brain_sim/annotation"
gtf <- file.path(ann_dir, "gencode.v26.annotation.gtf.gz")
txdb.filename <- file.path(ann_dir, "gencode.v26.annotation.sqlite")
#gtf <- "gencode.v26.annotation.gtf.gz"
#txdb <- makeTxDbFromGFF(gtf)
#saveDb(txdb, txdb.filename)
txdb <- loadDb(txdb.filename)

txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
all(rownames(quants) %in% txdf$TXNAME)
quants$GENEID <- txdf$GENEID[match(rownames(quants), txdf$TXNAME)]
quants <- quants[order(quants$GENEID),]

ebt <- exonsBy(txdb, by="tx", use.names=TRUE)
txp.len <- sum(width(ebt))
quants$Length <- txp.len[rownames(quants)]
quants$NumReads <- quants$TPM * quants$Length
quants$NumReads <- round(quants$NumReads * 50e6 / sum(quants$NumReads))
table(quants$NumReads > 10)

quantsUp <- quants ##Updated quants sim
quantsUp[quantsUp$NumReads < 10, c("TPM", "NumReads")] <- 0

# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz
fasta0 <- file.path(ann_dir, "gencode.v26.transcripts.fa.gz")
# fasta <- file.path(save_dir, "gencode.v26.transcripts.short.names.fa")
suppressPackageStartupMessages(library(Biostrings))
txseq <- readDNAStringSet(fasta0)
names(txseq) <- sub("\\|.*","",names(txseq))
# writeXStringSet(txseq, file=fasta)
#gc <- as.numeric(letterFrequency(txseq, "GC", as.prob=TRUE))

### simulate DTU and DTE
set.seed(seed)
# ratio of expressed genes with DGE, DTE or DTU

exprs.txs <- rownames(quantsUp)[quantsUp$TPM > 0]
exprs.genes <- unique(txdf$GENEID[match(exprs.txs, txdf$TXNAME)])
nsampgenes <- ifelse(nsampgenes > 15000, length(exprs.genes), nsampgenes)
samp.genes <- sample(exprs.genes, nsampgenes)
samp.txps <- txdf[txdf$GENEID %in% samp.genes,"TXNAME"]
quantsUp[setdiff(rownames(quantsUp), samp.txps), c("TPM", "NumReads")]=0
# re-order txdf for speed
txdf <- txdf[match(rownames(quantsUp),txdf$TXNAME),]
all.equal(rownames(quantsUp), txdf$TXNAME)

# this will keep the isoform TPMs for the two groups
tpms <- matrix(0,nrow=nrow(quantsUp), ncol=2)
rownames(tpms) <- rownames(quantsUp)
tpms[,1] <- quantsUp[["TPM"]]
tpms[,2] <- quantsUp[["TPM"]]

# make counts from TPMs
fraglen <- 200
per.nuc <- tpms * pmax(quantsUp[["Length"]] - fraglen, 1)
sim.counts.mat <- t(t(per.nuc) / colSums(per.nuc) * libSize)

# only simulate when expected count 5 or more in one group
keep <- sim.counts.mat[,1] >= 5
sim.counts.mat <- sim.counts.mat[keep,]

# subset the transcripts
all(rownames(sim.counts.mat) %in% names(txseq))
txseq <- txseq[ rownames(sim.counts.mat) ]
stopifnot(all(names(txseq) == rownames(sim.counts.mat)))
writeXStringSet(txseq, file.path(saveDir, paste("transcripts_", seed, ".fa", sep="")))


sim.means <- rowMeans(sim.counts.mat)
load("/fs/cbcb-lab/rob/students/noor/Uncertainity/treeTerminusScripts/brain_simulation/meanDispPairs/meanDispPairs.rda")
match.idx <- sapply(sim.means, function(mu) {
 which.min(abs(mu - meanDispPairs$baseMean))
})
disps <- meanDispPairs$dispGeneEst[match.idx]
plot(sim.means, disps, log="xy", cex=.1) # looks right

# log10 disp ~ N(-1,0.5)
#set.seed(1)
#disps <- 10^rnorm(nrow(sim.counts.mat),-1,0.5)

summary(disps)
summary(disps[sim.means > 1000])

library(alpine)
load("/fs/cbcb-lab/rob/students/noor/Uncertainity/treeTerminusScripts/brain_simulation/data/fitpar_all.rda")
fitpar[[1]][["model.params"]]$gc.knots <- seq(from=.4, to=.6, length=3)
fitpar[[1]][["model.params"]]$gc.bk <- c(0,1)
# plotGC(fitpar, "all", col=rep(1:2,each=15))
frag_GC_bias <- plotGC(fitpar, "all", return.type=2)
frag_GC_bias <- frag_GC_bias[,16:30]
#plot(0:100/100, frag_GC_bias[,1], type="l", ylim=c(0,1))
#for (i in 2:15) {
#  lines(0:100/100, frag_GC_bias[,i])
#}
sim_counts <- sim.counts.mat[,1]
fold_changes <- rep(1, nrow(sim.counts.mat))
save(fold_changes, tpms,
     txdf, sim.counts.mat, sim_counts,
     disps, frag_GC_bias, samp.txps, samp.genes,
     file=file.path(dataDir,paste("simulate_seed=", seed, ".rda", sep="")))

library(devtools)
si <- session_info()$packages
write.table(si, file="session_info.txt", quote=FALSE, sep="\t")

