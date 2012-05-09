####################################################################################
# Author: Aaron Brooks
# Affiliation: Institute for Systems Biology, Seattle, WA
# Date of creation: 04/19/2012
# Last update: 04/19/2012
####################################################################################
# DESCRIPTION
####################################################################################
# Main function. Calls accessory functions to post-process EGRIN2
####################################################################################
# REQUIRED PARAMETER DESCRIPTIONS
####################################################################################
# [e] reference to EGRIN environment
####################################################################################
# CUSTOMIZABLE PARAMETER DESCRIPTIONS####################################################################################
# NOTE: set these in the file main.R
# BACKBONE.PVAL [0.05] sig required by backbone extraction
# FULLY.CONNECTED [F] assume fully connected graph for backbone extraction. 
# MULTICORE [T] run on multiple cores
# CORES [4] how many cores to run on
# RESID.FILTER [1] filter out biclusters above given residual value
# MINMOTIFSCORE [Inf] filter out biclusters above given motif e value
# OUTDIR [./out/] where files should be written
# LINKCOMM.SCORE [0] use link similarity def of (0) Ahn or (1) Kalinka
# LINKCOMM.SIMTHRESINC [.1] amount to increment community detection threshold.
# LINKCOMM.SIMSCORE [5] score used to evaluate global density of communities (1,2,3,4,5)
####################################################################################


# PARAMS
HAL = F
BACKBONE.PVAL=0.05
FULLY.CONNECTED=F 
MULTICORE=T
CORES=4
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

params <- setdiff(ls()[grep("[[:upper:]]",ls())],ls()[grep("[[:lower:]]",ls())])

require(multicore)
options(cores=CORES)

# load EGRIN env
runCorems <- function(gBg=NULL) {
  #source("processEGRIN.R")
  o <- new.env(parent = baseenv())
  o$parameters <- lapply(params,function(i){eval(as.symbol(i))})
  names(o$parameters) <- params
  system(paste("mkdir",OUTDIR,sep=" "))
  system(paste("mkdir filehash"))
  if (is.null(gBg)){
    cat("Making gene-gene co-occurence matrix from cMonkey data\n")
    gBg <- make.r.gBg()
  } else {
    cat("Using user supplied gene-gene co-occurence matrix\n")
  }
  cat("Extracting backbone\n")
  gBg.backbone <- multiscaleBackbone(gBg)
  cat("Writing edge list\n")
  writeEdgeList(gBg.backbone)
  cat("Running corem detection\n")
  runCoremDetection(numGenes = dim(gBg)[1])
  o$link.community.threshold <- chooseCutoff()
  cat("Reading in corems\n")
  # unload filehashRO
  unload("filehashRO")
  require(filehash)
  dbCreate("./filehash/gg_filehash.dump")
  o$gg <- dbInit("./filehash/gg_filehash.dump")
  o$gg$gBg <- gBg; rm(gBg)
  o$gg$gBg.backbone <- gBg.backbone; rm(gBg.backbone)
  dbCreate("./filehash/corem_filehash.dump")
  o$corems <- dbInit("./filehash/corem_filehash.dump")
  o$corems$all <- loadCorems(o$link.community.threshold)
  o$corems$clean_density <- o$corems$all[o$corems$all[,Community.Weighted.Density]>0,]
  if (COREMBYSIZE) {
    o$corems$clean_size <- cleanCoremsBySize(o$corems$all)
  }
  o$gg$gBg.backbone.corems <- coremsTOgbg(o$corems[[COREMMETHOD]])
  # reload filehashRO
  unload("filehash")
  require(filehashRO)
  o$ratios <- normalizeRatios(e$ratios[[1]])
  save(o,file=RDATANAME)
  return(o)
}

loadEnv <- function() {
  load(RDATANAME)
  # re.init filehash
  unload("filehash")
  require(filehashRO)
  require(data.table)
  if (file.exists("./filehash/gg_filehash.dump")) {
    cat("Loading gene-gene co-occurence matrices into env$gg\n")
    o$gg <- dbInit("./filehash/gg_filehash.dump")
  }
  if (file.exists("./filehash/corem_filehash.dump")) {
    cat("Loading corem data.table into env$corems\n")
    o$corems <- dbInit("./filehash/corem_filehash.dump")
  }
  return(o)
}

processCorems <- function() {
  o<-loadEnv()
  o$corem_list <- list()
  o$corem_list$corems <- unique(o$corems[[COREMMETHOD]][,Community.ID])
  o$corem_list$genes <- lapply(o$corem_list$corems,function(i) getGenes(i,o$corems[[COREMMETHOD]]))
  names(o$corem_list$genes) <- o$corem_list$corems
  o$corem_list$conditions <- findCoremConditions.group(o$corem_list,o$ratios,ratios.normalized=T,
                                                       method=CONDITIONMETHOD,resamples=CONDITIONRESAMPLES,
                                                       all=F,padjust=F,pval=0.05,enforce.diff=F,
                                                       diff.cutoff=2,filehash=CONDITIONFILEHASH,lookup.table=lookup.table)
  save(o,file=RDATANAME)
}

analyzeCorems <- function() {
  
}

if (F) {
  
}