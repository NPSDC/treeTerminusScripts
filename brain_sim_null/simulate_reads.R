args <- commandArgs(trailingOnly = TRUE)
print(args)
i <- as.numeric(args[1])
fastaPath <- args[2]
outPath <- args[3]
seed <- as.numeric(args[4])
sim_path <- args[5]

library(polyester)
set.seed(seed+i)
load(sim_path)
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
