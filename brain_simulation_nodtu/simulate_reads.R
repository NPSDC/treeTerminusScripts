args <- commandArgs(trailingOnly = TRUE)
print(args[1])
i <- as.numeric(args[1])
fastaPath <- args[2]
outPath <- args[3]
print(paste(i, fastaPath, outPath))
library(polyester)
load("../environment/brain_sim_nodtu/simulate.rda")
se <- simulate_experiment(fasta=fastaPath,
                          outdir=file.path(outPath,  i),
                          num_reps=c(1,1),
                          reads_per_transcript=sim_counts,
                          size=1/disps,
                          fold_changes=fold_changes,
                          fraglen=200,
                          fragsd=25,
                          frag_GC_bias=frag_GC_bias[,c(i,i)],
                          readlen=100,
                          seed=i)
