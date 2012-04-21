####################################################################################
# Author: Aaron Brooks
# Affiliation: Institute for Systems Biology, Seattle, WA
# Date of creation: 04/19/2012
# Last update: 04/19/2012
####################################################################################
# DESCRIPTION
####################################################################################
# Functions that take ensemble cMonkey run (from Dave),
# transform into gene-by-gene co-occurence matrix,
# remove noisy edges by backbone extraction,
# write to edge list file,
# run linkcommunity detection to find corems,
# process corems,
# computes basic statistics
####################################################################################
# REQUIRED PARAMETER DESCRIPTIONS
####################################################################################
# [e] reference to EGRIN environment
####################################################################################
# CUSTOMIZABLE PARAMETER DESCRIPTIONS
####################################################################################
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
# COREMSIZETHRESH [3] minimum size of corem, # edges
####################################################################################

make.r.gBg <- function(clusterStack = e$clusterStack,filt.resid=RESID.FILTER,minMotifScore=MINMOTIFSCORE) {
  require(multicore)
  R.m = matrix(0,nrow=length(rownames(e$ratios[[1]])),ncol=length(rownames(e$ratios[[1]])))
  rownames(R.m) = rownames(e$ratios[[1]]); colnames(R.m) = rownames(e$ratios[[1]]) 
  l = length(e$clusterStack)
  for(i in seq(1:length(e$clusterStack))) {
    #print(i)
    if (i%%5000==0) {
      cat(paste(signif((i/l)*100,2),"% complete\n",sep=""))
    }
    if (is.na(e$clusterStack[[i]]$resid)) {
        } else if (e$clusterStack[[i]]$resid<filt.resid) {
          if (is.na(e$clusterStack[[i]]$e.val)){
          } else if (e$clusterStack[[i]]$e.val < minMotifScore) {
            g<-t(combn(e$clusterStack[[i]]$rows,2))
            R.m[g] = R.m[g]+1
            R.m[cbind(g[,2],g[,1])] = R.m[cbind(g[,2],g[,1])] + 1
          }
        }
  }
  # make diag = 0
  diag(R.m) <- 0
  # normalize rows
  R.m.2 <- t(apply(R.m,1,function(x){if (sum(x)>0) {return(x/sum(x))} else {return(x/Inf)} }))
  return(R.m.2)
}

multiscaleBackbone <- function(gBg, pval=BACKBONE.PVAL,multicore=MULTICORE,makeFullyConnected=FULLY.CONNECTED) {
  # gBg is any arbitrary matrix
  require(multicore)
  if (makeFullyConnected) {
    k <- sum(apply(gBg,1,sum)>0)
  }
  # define uniform random integrand
  integrand <- function(x,k) {(1-x)^(k-2)}
  # fcn to calculate degree,k
  calc_k <- function(i) {sum(i>0)}
  # only bother with non-zero entries (also to determine k)
  index <- which(gBg>0,arr.ind=T)
  # score each non-zero entry, computes p that defines "significance of edge"
  if (multicore) {
    a_ij <- unlist( mclapply( seq( 1:dim(index)[1] ), function(i) { 
      # Calc degree
      if (!makeFullyConnected) {
        k<-calc_k( gBg[ index[i,1], ] ) ;
      } 
      # Integration to determine p
      o <- 1-(k-1)*integrate(integrand,0,gBg[ index[i,1], index[i,2] ],k=k)$value; 
      return(o) 
    } ) )
  } else {
    a_ij <- unlist( lapply( seq( 1:dim(index)[1] ), function(i) { 
      # Calc degree
      if (!makeFullyConnected) {
        k<-calc_k( gBg[ index[i,1], ] ) ;
      } 
      # Integration to determine p
      o <- 1-(k-1)*integrate(integrand,0,gBg[ index[i,1], index[i,2] ],k=k)$value; 
      return(o) 
    } ) )
  }
  # Select edges with p less than pval cutoff
  a_ij_sig <- a_ij <= pval
  index_sig <- index[a_ij <= pval, ]
  o <- gBg; o[] <- 0
  o[index_sig] <- gBg[index_sig]
  # Make symmetric. If two undirected edges are both sig, keep highest weight. 
  upper.index <- which(upper.tri(o),arr.ind=T)
  upper.index.sym <- cbind(upper.index[,2],upper.index[,1])
  max.weight <- apply(cbind(o[upper.index],o[upper.index.sym]),1,max)
  o[upper.index] = max.weight
  o[upper.index.sym] = max.weight
  return(o)
}

writeEdgeList <- function(matrix,file=paste(OUTDIR,"out",sep=""), weighted=T) {
  index <- which(upper.tri(matrix),arr.ind=T)
  index <- index[which(matrix[index]>0),]
  namesI <- rownames(matrix)
  namesJ <- colnames(matrix)
  for (i in 1:dim(index)[1]) {
    if (weighted) {
      write(x=paste(namesI[index[i,1]],namesJ[index[i,2]],matrix[index[i,1],index[i,2]]),file=file,append=T)
    } else {
      write(x=paste(namesI[index[i,1]],namesJ[index[i,2]]),file=file,append=T)
    }
  }
}

runCoremDetection <- function(numGenes = dim(e$ratios[[1]])[1], dir = OUTDIR,s = LINKCOMM.SCORE) {
  require(multicore)
  cwd <- getwd()
  if (length(system("which adjmat2wpairs",inter=T))==0) {
    print("Please compile adjmat2wpairs")
    return(NULL)
  }
  setwd(OUTDIR)
  system("adjmat2wpairs out 0 0")
  if (length(system("which compute_tanimoto",inter=T))==0) {
    print("Please compile compute_tanimoto")
    return(NULL)
  }
  empty<-mclapply(seq(1,CORES),function(i){
    ref <- round(seq(0,numGenes+1,length.out = CORES+1))
    system(paste("compute_tanimoto out",s,ref[i],ref[i+1],sep=" "))
    invisible(NULL)
  })
  system("cat out.tanimoto_* > out.tanimoto")
  empty<-mclapply(seq(0,1,LINKCOMM.SIMTHRESINC),function(i){
    system(paste("cluster_communities out",i,sep=" "))
    invisible(NULL)
  })
  system("cat out.density_* > out.density")
  density <- read.table("out.density",header=F)
  def_map = c(3,5,7,8,9); names(def_map) <- c(1,2,5,3,4)
  ind <- as.integer(def_map[as.character(LINKCOMM.SIMSCORE)])
  maxT <- density[,1][which(density[,ind]==max(density[,ind]))]
  cat("Summary plotted. See density_stats.pdf\n")
  ######
  # Plot
  ######
  pdf("density_stats.pdf")
  plot(density[,1],density[,3],col="black",type="l",lty=1,
       ylim=c(0, 1.1*max(c(density[,3],density[,5],density[,8]),na.rm=T)),ylab="Similarity Score",
       xlab="Threshold",yaxt="n")
  lines(density[,1],density[,5],type="l",lty=2)
  lines(density[,1],density[,8],type="l",lty=3)
  axis(2, pretty(c(0, 1.1*max(c(density[,3],density[,5],density[,8]),na.rm=T))), col='black')
  par(new=T)
  plot(density[,1],density[,9],col="blue",type="l",lty=1,axes=F,
       ylim=c(0, 1.1*max(c(density[,7],density[,9]),na.rm=T)),ylab="",xaxt="n",xlab="")
  lines(density[,1],density[,7],col="blue",type="l",lty=2)
  axis(4, pretty(c(0, 1.1*max(c(density[,7],density[,9]),na.rm=T))), col='blue')
  lines(density[,1],density[,ind],type="l",lty=1,lwd=2,col="red")
  lines(c(maxT,maxT),c(-1,1),type="l",lty=2,col="red")
  legend("topleft",c("(1) Unweighted (all)","(2) Unweighted (n>2)","(3) Weighted (local)",
                     "(4) Weighted (global)","(5) Weighted (both)",
                     paste("Your choice",paste("Max =",maxT,sep=" "),sep=": ")),
         col=c("black","black","black","blue","blue","red"),lty = c(1,2,3,1,2,1),lwd=1)
  dev.off()
  setwd(cwd)
  cat(paste("Threshold is", maxT,"\n",sep=" "))
  return(maxT)
}

loadRegulons <- function(fCutoff,fRoot="out",dataTable=T) {
  require(data.table)
  # Load/Process/Save file
  # Call getting_communities as follows:
  # getting_communities ROOT_OF_FILENAME THRESHOLD
  fName <- paste(OUTDIR,fRoot,".communities_",fCutoff,sep="")
  # Process Raw file
  cwd <- getwd()
  setwd(OUTDIR)
  system(paste("getting_communities",fRoot,fCutoff))
  setwd(cwd)
  out <- read.table(fName)
  colnames(out) <- 
    c("Gene1","Gene2","Community.ID","Community.Density","Community.Weighted.Density")
  out[,3] <- as.character(out[,3])
  if (dataTable) {
    require(data.table)
    out <- data.table(out)
    setkey(out,"Community.ID")
  }
  return(out)
}

cleanCoremsBySize <- function(regulons.table, threshold = COREMSIZETHRESH,gene=F) {
  # Remove regulons with fewer than x genes
  require(multicore)
  require(data.table)
  setkey(regulons.table,"Community.ID")
  if (gene) {
    # By gene size
    r <- as.character(unique(regulons.table[,Community.ID]))
    r.len <- unlist(lapply(r,function(i){i<-length(getGenes(i,regulons.table))}))
    names(r.len) <- r
    toRemove <- names(r.len)[which(r.len <= threshold)]
    o <- regulons.table
    for (i in toRemove) {
      print(i)
      o <- o[-which(o[,Community.ID]==i),]  
    }
  } else {
    # By edges
    r <- table(regulons.table[,Community.ID])
    r <- names(r[r<=threshold])
    toRemove <- unlist(mclapply(r,function(i){i<-which(regulons.table[,Community.ID]==i)}))
    o <- regulons.table
    o <- o[-toRemove,]
  }
  return(o)
}

