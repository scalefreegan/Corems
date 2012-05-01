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
BACKBONE.PVAL=0.05
FULLY.CONNECTED=F 
MULTICORE=T
CORES=4
RESID.FILTER=1
MINMOTIFSCORE=Inf
OUTDIR="./out/"
LINKCOMM.SCORE=0
LINKCOMM.SIMTHRESINC=.1
LINKCOMM.SIMSCORE=5
COREMSIZETHRESH = 3

params <- setdiff(ls()[grep("[[:upper:]]",ls())],ls()[grep("[[:lower:]]",ls())])

# load EGRIN env
runCorems <- function() {
  #source("processEGRIN.R")
  o <- new.env(parent = baseenv())
  o$parameters <- lapply(params,function(i){eval(as.symbol(i))})
  names(o$parameters) <- params
  system(paste("mkdir",OUTDIR,sep=" "))
  system(paste("mkdir filehash"))
  cat("Making gene-gene co-occurence matrix from cMonkey data\n")
  o$gBg <- make.r.gBg()
  cat("Extracting backbone\n")
  o$gBg.backbone <- multiscaleBackbone(o$gBg)
  cat("Writing edge list\n")
  writeEdgeList(o$gBg.backbone)
  cat("Running corem detection\n")
  o$link.community.threshold <- runCoremDetection()
  cat("Reading in corems")
  # unload filehashRO
  unload("filehashRO")
  require(filehash)
  dbCreate("./filehash/corem_filehash.dump")
  o$corems <- dbInit("./filehash/corem_filehash.dump")
  o$corems$all <- loadCorems(o$link.community.threshold)
  o$corems$clean_density <- o$corems$all[o$corems$all[,Community.Weighted.Density]>0,]
  o$corems$clean_size <- cleanCoremsBySize(o$corems$all)
  # reload filehashRO
  unload("filehash")
  require(filehashRO)
  o$ratios <- normalizeRatios(e$ratios[[1]])
}

processCorems <- function() {
  
}

analyzeCorems <- function() {
  
}

if (F) {
  
}