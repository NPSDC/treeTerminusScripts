suppressPackageStartupMessages(source("tree_helper_function.R"))
load("environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/yAll.RData")
load("environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/y.RData")
load("environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/seBrainSim.RData")
load("environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/tree.RData")

yNS <- seBrainSim[tree$tip,] ##Not scaled
yNS <- computeInfRV(yNS, meanVariance=F)
yAggNS <- prepSwish(tree, yNS)
yAggNS <- computeInfRV(yAgg, meanVariance=F)

yAgg <- prepSwish(tree, y)
yS <- y
mcols(yAgg)[["meanInfRV"]] <- mcols(yAggNS)[["meanInfRV"]]
mcols(yS)[["meanInfRV"]] <- mcols(yNS)[["meanInfRV"]]
save(yS, file="environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/updated/yS.RData")
save(yAgg, file="environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/updated/yAgg.RData")
save(yNS, file="environment/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100/updated/yNS.RData")