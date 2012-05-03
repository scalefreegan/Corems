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

####################################################################################
# Discover corems given ensemble cMonkey
####################################################################################

# 
# Transform cMonkey output into weighted adjacency matrix
# co-occurence of genes in biclusters
#
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

#
# Extract the backbone from the graph
# removes noise
# From PMID: 19357301
#
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

#
# Write clean gBg matrix in format that can be read by Antoine's
# corem detection code
#
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

#
# Run Antoine's C++ corem detection code
#
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
  setwd(cwd)
}

chooseCutoff <- function() {
  cwd <- getwd()
  setwd(OUTDIR)
  density <- read.table("out.density",header=F)
  def_map = c(3,5,7,8,9); names(def_map) <- c(1,2,5,3,4)
  ind <- as.integer(def_map[as.character(LINKCOMM.SIMSCORE)])
  maxT <- density[,1][which(density[,ind]==max(density[,ind]))]
  cat("Summary plotted. See density_stats.pdf\n")
  ######################################################################
  # Plot corem density distributions across cut height given alternative
  # scoring approaches
  ######################################################################
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

# 
# Loads corems with max density according to scoring metric used
#
loadCorems <- function(fCutoff,fRoot="out",dataTable=T) {
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

####################################################################################
# Postprocess corems
####################################################################################

cleanCoremsBySize <- function(corems.table, threshold = COREMSIZETHRESH,gene=F) {
  # Remove corems with fewer than x genes
  require(multicore)
  require(data.table)
  setkey(corems.table,"Community.ID")
  if (gene) {
    # By gene size
    r <- as.character(unique(corems.table[,Community.ID]))
    r.len <- unlist(lapply(r,function(i){i<-length(getGenes(i,corems.table))}))
    names(r.len) <- r
    toRemove <- names(r.len)[which(r.len <= threshold)]
    o <- corems.table
    for (i in toRemove) {
      print(i)
      o <- o[-which(o[,Community.ID]==i),]  
    }
  } else {
    # By edges
    r <- table(corems.table[,Community.ID])
    r <- names(r[r<=threshold])
    toRemove <- unlist(lapply(r,function(i){i<-which(corems.table[,Community.ID]==i)}))
    o <- corems.table
    o <- o[-toRemove,]
  }
  return(o)
}

coremsTOgbg <- function(corem.table) {
  # Make corem data.table into gBg co-occurence matrix for analysis
  g<-unique(c(as.character(corem.table[,Gene1]),as.character(corem.table[,Gene2])))
  m <- matrix(0,nrow=length(g),ncol=length(g),dimnames=list(g,g))
  m[cbind(as.character(corem.table[,Gene1]),as.character(corem.table[,Gene2]))] <- corem.table[,Community.Density]
  m[cbind(as.character(corem.table[,Gene2]),as.character(corem.table[,Gene1]))] <- corem.table[,Community.Density]
  return(m)
}

#Get genes from corems function
getGenes <- function(coremID = "1", corems.table = corems) {
  if (class(corems.table)[1]=="data.table") {
    require(data.table)
    setkey(corems.table,"Community.ID")
    g <- corems.table[coremID,mult="all"]
    g <- unique(c(as.character(g[,Gene1]),as.character(g[,Gene2])))
  } else {
    g<-unique(c(as.character(corems.table[which(corems.table[,3]==coremID),1]),
                as.character(corems.table[which(corems.table[,3]==coremID),2])))
  }
  g <- g[!is.na(g)]
  return(g) 
}

getcorems <- function(geneName = "VNG0700G", corems.table = corems) {
  if (class(corems.table)[1]=="data.table") {
    require(data.table)
    setkey(corems.table,Gene1)
    g1 <- unique(as.character(corems.table[geneName,mult="all"][,Community.ID]))
    setkey(corems.table,Gene2)
    g2 <- unique(as.character(corems.table[geneName,mult="all"][,Community.ID]))
    g <- unique(c(g1,g2))
    setkey(corems.table,Community.ID)
  } else {
    g<-unique(c(as.numeric(corems.table[which(corems.table[,1]==geneName),3]),
                as.numeric(corems.table[which(corems.table[,2]==geneName),3])))
  }
  g <- g[!is.na(g)]
  return(g) 
}

resampleRandomConditions <- function(geneSetSize=seq(3,200,1),ratios,resamples=20000,
                                     method=c("sd","cvar")[2],mode="none",filehash=T) {
  require(multicore)
  require(Matrix)
  if (filehash) {
    unload("filehashRO")
    require(filehash)
    fn <- paste("./filehash/corem_",paste(method,resamples,"filehash",sep="_"),".dump",sep="")
    dbCreate(fn)
    o <- dbInit(fn)
  } else {
    o <- list()
  }
  genePool <- rownames(ratios)
  geneSetSize <- as.character(geneSetSize)
  if (method == "sd") {
    o$sd<-mclapply(seq(1,length(geneSetSize)),function(i) {
      len = as.integer(geneSetSize[i])
      #print(len)
      if (i%%10==0) {
        cat(paste(signif((i/length(geneSetSize))*100,2),"% complete\n",sep=""))
      }
      i <- do.call(rbind,lapply(seq(1:resamples),function(j){apply(ratios[sample(genePool,len),colnames(ratios)],2,sd,na.rm=T)}))
      i.2 <- do.call(cbind,lapply(seq(1:dim(i)[2]),function(j){
        j <- i[,j]
        j.ecdf <- ecdf(j)
        # make everything above lowest 7.5% quantile = 0
        j[which(j>quantile(j.ecdf,probs=seq(0,1,.075))[2])] = 0
        return(j)
      }))
      colnames(i.2) <- colnames(i)
      i.2 <- as(i.2,"sparseMatrix")    
    }) 
    names(o$sd) <- as.character(geneSetSize)
  } else if (method == "cvar") {
    o$cvar<-mclapply(seq(1,length(geneSetSize)),function(i) {
      len = as.integer(geneSetSize[i])
      #print(len)
      if (i%%10==0) {
        cat(paste(signif((i/length(geneSetSize))*100,2),"% complete\n",sep=""))
      }
      i <- do.call(rbind,lapply(seq(1:resamples),function(j){cvar(genes=sample(genePool,len),conditions=colnames(ratios),ratios=ratios,mode=mode)}))
      i.2 <- do.call(cbind,lapply(seq(1:dim(i)[2]),function(j){
        j <- i[,j]
        j.ecdf <- ecdf(j)
        # make everything above lowest 7.5% quantile = 0
        j[which(j>quantile(j.ecdf,probs=seq(0,1,.075))[2])] = 0
        return(j)
      }))
      colnames(i.2) <- colnames(i)
      i.2 <- as(i.2,"sparseMatrix")
    })
    names(o$cvar) <- as.character(geneSetSize)
  }
  unload(filehash)
  require(filehashRO)
  invisible(o)
}

findCoremConditions.ind <- function(genes,ratios,ratios.normalized=F,method=c("sd","cvar")[2],resamples=20000,
                                       all=F,padjust=F,pval=0.05,enforce.diff=F,diff.cutoff=2,filehash=T,lookup.table=NULL...) {
  require(multicore)
  len = as.character(length(genes))
  if (filehash&&is.null(lookup.table)) {
    fn <- paste("./filehash/corem_",paste(method,resamples,"filehash",sep="_"),".dump",sep="")
    lookup.table <- dbInit(fn)
  } else if (is.null(lookup.table)) {
    lookup.table <-resampleRandomConditions(geneSetSize=len,ratios,resamples,method,"none",F)
  }
  if (method == "sd") {
    # Make 0 values Inf -- from sparsification above
    tmp.table = lookup.table[[method]][[len]]
    tmp.table[which(tmp.table==0)] = Inf
    lookup.ecdf <- apply(tmp.table,2,ecdf)
    val <- apply(ratios[genes,],2,sd,na.rm=T)
    o <- sapply(seq(1,length(val)),function(i){lookup.ecdf[[names(val)[i]]](val[i])})
  } else if (method == "cvar") {
    tmp.table = lookup.table[[method]][[len]]
    tmp.table[which(tmp.table==0)] = Inf
    lookup.ecdf <- apply(tmp.table,2,ecdf)
    val <- cvar(genes,conditions=colnames(ratios),ratios=ratios,mode="none")
    o <- sapply(seq(1,length(val)),function(i){lookup.ecdf[[names(val)[i]]](val[i])})
  }
  if (padjust) {
    o <- p.adjust(o,method="BH")
  }
  if (!all) {
    # Only report conds with p<=pval
    o <- o[o<=pval]
  }
  if (enforce.diff) {
    # normalize ratios (just in case they're not)
    if (ratios.normalized == F) {
      ratios <- normalizeRatios(ratios)
    }
    meanExp <- abs(colMeans(ratios[genes,names(o)]))
    meanExp <- which(meanExp>=diff.cutoff)
    o <- o[meanExp]
  }
  return(o)
}

findCoremConditions.group <- function(coremStruct,ratios,ratios.normalized=F,method=c("sd","cvar")[2],resamples=20000,
                                    all=F,padjust=F,pval=0.05,enforce.diff=F,diff.cutoff=2,filehash=T,lookup.table=NULL) {
  # Corem struct is:
  # env$corem_list
  require(multicore)
  if (filehash) {
    # store random resamples in filehash
    # WARNING: This may be VERY large file. >50GB
    if (is.null(lookup.table)) {
      cat("Couldn't find precomputed resamples. Computing now. This might take awhile. \n")
      lookup.table<-resampleRandomConditions(geneSetSize=sort(unique(sapply(coremStruct$genes,length))),
                             o$ratios,resamples=resamples,method=method,mode="none")
    } 
    cat("Using user supplied precomputed resamples\n")
    o <- mclapply(seq(1,length(coremStruct$corems)),function(i) {
      if (i%%10==0) {
        cat(paste(signif((i/length(coremStruct$corems))*100,2),"% complete\n",sep=""))
      }
      out<-findCoremConditions.ind(coremStruct$genes[[coremStruct$corems[[i]]]],ratios,T,method,resamples,
                          all,padjust,pval,enforce.diff,diff.cutoff,filehash,lookup.table=lookup.table)
      return(out)
      })
    names(o) <- coremStruct$corems
  } else {
    ###
    # UNFINISHED!!
    ###
    # compute enrichments w/o storing in filehash
    # find lengths of corems
    c.len <- sapply(coremStruct$genes,len)
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
  return(o)
}


####################################################################################
# General functions
####################################################################################

# Unloads packages. Mostly to unload filehashRO
unload<-function(pkg) {if (!is.na(match(paste("package", pkg, sep=":"), search()))) detach(pos = match(paste("package", pkg, sep=":"), search()),unload=T)}

# Normalize ratios
normalizeRatios <- function(ratios){
  ratios.norm <- t(scale(t(ratios), center = apply(ratios, 1, median, na.rm = T), scale = apply(ratios, 1, sd, na.rm = T)))
  return(ratios.norm)
}

###############################
# Coefficient of variation
###############################
cvar <- function(genes,conditions,ratios,mode=c("median","mean","none")[1]) {
  m <- abs(apply(ratios[genes,conditions,drop=F],2,mean,na.rm=T))
  m[which(m==0)] = 1e-6
  if (mode=="mean") { 
    var.m <- mean(apply(ratios[genes,conditions,drop=F],2,sd,na.rm=T)/m,na.rm=T)
  } else if (mode=="median"){
    var.m <- median(apply(ratios[genes,conditions,drop=F],2,sd,na.rm=T)/m,na.rm=T)
  } else {
    var.m <- apply(ratios[genes,conditions,drop=F],2,sd,na.rm=T)/m
  }
  return(var.m)
}

###############################
# Residuals                
###############################
residual <- function( rats, allRatios, maxRowVar = mean( apply( allRatios[,], 1, var, use="pair" ), na.rm=T )) {
  d.rows <- rowMeans( rats, na.rm=T )
  d.cols <- colMeans( rats, na.rm=T )
  d.all <- mean( d.rows, na.rm=T )
  rij <- rats + d.all
  rij <- rij - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
  rij <- rij - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
  average.r <- mean( abs( rij ), na.rm = TRUE )
  row.var <- mean( apply( rats, 1, var, use = "pairwise.complete.obs" ), na.rm=T )
  if ( is.na( row.var ) || row.var > maxRowVar ) row.var <- maxRowVar
  average.r <- average.r / row.var
  return(average.r)
}

###############################
# General plotting              
###############################

# Microarray color scheme
require(gplots)
blue2yellow <- colorpanel(200,"blue","black","yellow")

plotHeatmap <- function(matrix,file="",cor=F,RowSideColors=rep("red",dim(matrix)[1]),scale="none",
                        dendrogram="row",cexRow=0.75,cexCol=0.75,...) {
  require(gplots)
  if (scale != "none") {
    mi.col <- colorRampPalette(c("white","red"))
    heatmap.2(matrix,trace="none",cexRow=cexRow,cexCol=cexCol,
              col=bluered(256),density.info="none",RowSideColors=RowSideColors,scale=scale,dendrogram=dendrogram,...)
  } else if (cor == F) {
    mi.col <- colorRampPalette(c("white","red"))
    heatmap.2(matrix,trace="none",cexRow=cexRow,cexCol=cexCol,
              breaks=seq(0,1,.01),col=mi.col(100),density.info="none",RowSideColors=RowSideColors,scale=scale,
              dendrogram=dendrogram,...)
  } else {
    heatmap.2(matrix,trace="none",cexRow=cexRow,cexCol=cexCol,
              breaks=seq(-1,1,.02),col=bluered(100),density.info="none",RowSideColors=RowSideColors,scale=scale,
              dendrogram=dendrogram,...)
  }
}