# how do corems perform on DREAM5 data
# which corems so signs of activity in new data set
# are the number of active conditions greater than random
# how does distribution of active conditions for corems in DREAM5 data set compare to Distiller data set

# starting dir: /home/abrooks/Documents/EGRIN2/Eco_ensemble_2/validation/DREAM5

# source scripts
source("~/Documents/git/Corems/processEGRIN.R")

# load data
load("DREAM5_normed_ratios.RData")

# make a random condition filehash for DREAM5 data 
resampleHash <- resampleRandomConditions(geneSetSize=seq(3,20,1),ratios,resamples=20000,method="cvar",mode="none",filehash=T)

# DREAM5 gBg compared to Distiller gBg
# load data
load("./validation/DREAM5/DREAM5_gene_gene_for_aaron.RData")
# compute common genes
cg <- intersect(rownames(e$ratios[[1]]),rownames(gene.gene)) 
# compute backbone
# normalize gene.gene
R.m <- gene.gene
# make diag = 0
diag(R.m) <- 0
R.m.2 <- t(apply(R.m,1,function(x){if (sum(x)>0) {return(x/sum(x))} else {return(x/Inf)} }))
gene.gene <- R.m.2
gene.gene.backbone <- multiscaleBackbone(gene.gene, pval=0.05,multicore=F,makeFullyConnected=F)
# correlation 
gBg.cor <- cor(as.vector(o$gg$gBg[cg,cg]),as.vector(gene.gene[cg,cg]))
gBg.backbone.cor <- cor(as.vector(o$gg$gBg.backbone[cg,cg]),as.vector(gene.gene.backbone[cg,cg]))

save(gene.gene,gene.gene.backbone,gBg.cor,gBg.backbone.cor,file="./validation/DREAM5/DREAM5_gBg.RData")

##
# make corems for DREAM5 network
##
# parameters
BACKBONE.PVAL=0.05
FULLY.CONNECTED=F 
MULTICORE=F
CORES=1
RESID.FILTER=1
MINMOTIFSCORE=Inf
OUTDIR="./out/"
LINKCOMM.SCORE=0
LINKCOMM.SIMTHRESINC=.1
LINKCOMM.SIMSCORE=3
COREMSIZETHRESH = 3
COREMBYSIZE = F
RDATANAME = "corems.RData"
CONDITIONRESAMPLES = 20000
CONDITIONMETHOD = "cvar"
CONDITIONFILEHASH = T
COREMMETHOD = c("all","clean_density","clean_size")[2]
HAL = F
# run it
source("~/Documents/git/Corems/processEGRIN.R")
source("~/Documents/git/Corems/main.R")
load("DREAM5_gBg.RData")
load("DREAM5_normed_ratios.RData")
params <- setdiff(ls()[grep("[[:upper:]]",ls())],ls()[grep("[[:lower:]]",ls())])
o <- runCorems(gBg=gene.gene,ratios=ratios)