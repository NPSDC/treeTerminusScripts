library(TreeSummarizedExperiment)
library(ggplot2)
library(ggtree)
library(ggpubr)

p1 <- ggtree(ape::read.tree(text="((4,5), ((1,2),3));")) + geom_tiplab()

p2 <- ggtree(ape::read.tree(text="((1,2), (5,3));")) + geom_tiplab()

p13 <- ggtree(ape::read.tree(text="(1,5);")) + geom_tiplab()
p23 <- ggtree(ape::read.tree(text="(2,3);")) + geom_tiplab()
p3 <- ggarrange(p13,p23,nrow=2,ncol=1)

p14 <- ggtree(ape::read.tree(text="((1,2),3);")) + geom_tiplab()
p24 <- ggtree(ape::read.tree(text="(4,5);")) + geom_tiplab()
p4 <- ggarrange(p14,p24, nrow=2,ncol=1)

p5 <- ggtree(ape::read.tree(text="(((1,3),2),5);")) + geom_tiplab()

pOld <- ggarrange(p1,p2,p3,p4,p5, nrow=1, ncol=5)

###New
p1 <- ggtree(ape::read.tree(text="((4,5), ((1,2),3));")) + geom_tiplab()

p2 <- ggtree(ape::read.tree(text="(((1,2), (5,3)),4);")) + geom_tiplab()

p3 <- ggtree(ape::read.tree(text="((2,3), (1,5), 4);")) + geom_tiplab()

p4 <- ggtree(ape::read.tree(text="((3,(1,2)), (4,5));")) + geom_tiplab()

p5 <- ggtree(ape::read.tree(text="((((1,3),2),5),4);")) + geom_tiplab()

pNew <- ggarrange(p1, p2, p3, p4, p5, nrow = 1)
ggarrange(pOld, pNew, nrow = 2)
ggsave(pOld, filename = "~/cbcb_rob/Uncertainity/treeTerminusScripts/TreeTerminusImages/Trees/tree1.png", dpi = 600)
ggsave(pNew, filename = "~/cbcb_rob/Uncertainity/treeTerminusScripts/TreeTerminusImages/Trees/tree2.png", dpi = 600)
