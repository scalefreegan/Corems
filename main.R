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
LINKCOMM.SIMSCORE=3
COREMSIZETHRESH = 3
RDATANAME = "corems.RData"
CONDITIONRESAMPLES = 20000
CONDITIONMETHOD = "cvar"
CONDITIONFILEHASH = T

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
  gBg <- make.r.gBg()
  cat("Extracting backbone\n")
  gBg.backbone <- multiscaleBackbone(gBg)
  cat("Writing edge list\n")
  writeEdgeList(gBg.backbone)
  cat("Running corem detection\n")
  runCoremDetection()
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
  o$corems$clean_size <- cleanCoremsBySize(o$corems$all)
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
  # unload filehashRO
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

processCorems <- function(method=c("all","clean_density","clean_size")[2],filehash=CONDITIONFILEHASH) {
  o<-loadCorems()
  o$corem_list <- list()
  o$corem_list$corems <- unique(o$corems[[method]][,Community.ID])
  o$corem_list$genes <- lapply(o$corem_list$corems,function(i) getGenes(i,o$corems[[method]]))
  if (filehash) {
    # store random resamples in filehash
    # WARNING: This may be VERY large file. >50GB
    resampleRandomConditions(geneSetSize=sort(unique(sapply(o$corem_list$genes,length))),
                       o$ratios,resamples=CONDITIONRESAMPLES,method=CONDITIONMETHOD,mode="none")
    o$corem_list$conditions <- mclapply(o$corem_list$corems,function(g) {
      findCoremConditions(o$corem_list$genes[[g]],o$ratios,ratios.normalized=T,method=CONDITIONMETHOD,resamples=CONDITIONRESAMPLES,
                          all=F,padjust=F,pval=0.05,enforce.diff=F,diff.cutoff=2,filehash=T,lookup.table=NULL...)
  } else {
    # compute enrichments w/o storing in filehash
    # find lengths of corems
    c.len <- sapply(o$corem_list$genes,len)
    c.len.unique <- unique(c.len)
    o$corem_list$conditions <- mclapply(seq(1,length(c.len.unique),1),function(i){
      if (i%%4 == 0) {
        cat(paste(signif((i/length(c.len.unique))*100,2),"% complete\n",sep=""))
      }
      len = c.len.unique[i]
      print(len)
      # find corems with this length
      c.tmp <- o$corem_list$corems[[which(c.len)==len]]
      # resample genes 
      lookup.table <- resampleRandomConditions
    })
  }
  
}

analyzeCorems <- function() {
  
}

if (F) {
  
}