#########################################################################################
#
# Core functions
#
#########################################################################################

require(gplots)
blue2yellow <- colorpanel(200,"blue","black","yellow")




#Get genes from corem function
getCoremGenes <- function(coremID = "1", corem.table = corem) {
  if (class(corem.table)[1]=="data.table") {
    require(data.table)
    setkey(corem.table,"Community.ID")
    g <- corem.table[regulonID,mult="all"]
    g <- unique(c(as.character(g[,Gene1]),as.character(g[,Gene2])))
  } else {
    g<-unique(c(as.character(corem.table[which(corem.table[,3]==regulonID),1]),
                as.character(corem.table[which(corem.table[,3]==regulonID),2])))
  }
  g <- g[!is.na(g)]
  return(g) 
}

# Remove corems with fewer than x genes
cleanCoremsBySize <- function(corem.table = corem, threshold = 3,gene=F) {
  require(multicore)
  require(data.table)
  setkey(corem.table,"Community.ID")
  if (gene) {
    # By gene size
    r <- as.character(unique(corem.table[,Community.ID]))
    r.len <- unlist(lapply(r,function(i){i<-length(getGenes(i,corem.table))}))
    names(r.len) <- r
    toRemove <- names(r.len)[which(r.len <= threshold)]
    o <- corem.table
    for (i in toRemove) {
      print(i)
      o <- o[-which(o[,Community.ID]==i),]  
    }
  } else {
    # By edges
    r <- table(corem.table[,Community.ID])
    r <- names(r[r<=threshold])
    toRemove <- unlist(mclapply(r,function(i){i<-which(corem.table[,Community.ID]==i)}))
    o <- corem.table
    o <- o[-toRemove,]
  }
  return(o)
}



getCorems <- function(geneName = "VNG0700G", corem.table = corem) {
  if (class(corem.table)[1]=="data.table") {
    require(data.table)
    setkey(corem.table,Gene1)
    g1 <- unique(as.character(corem.table[geneName,mult="all"][,Community.ID]))
    setkey(corem.table,Gene2)
    g2 <- unique(as.character(corem.table[geneName,mult="all"][,Community.ID]))
    g <- unique(c(g1,g2))
    setkey(corem.table,Community.ID)
  } else {
    g<-unique(c(as.numeric(corem.table[which(corem.table[,1]==geneName),3]),
                as.numeric(corem.table[which(corem.table[,2]==geneName),3])))
  }
  g <- g[!is.na(g)]
  return(g) 
}

getCoremEdges <- function(coremID = "1", corem.table = corem) {
  if (class(corem.table)[1]=="data.table") {
    require(data.table)
    setkey(corem.table,"Community.ID")
    g <- corem.table[as.character(coremID),mult="all"]
    g <- cbind(as.character(g$"Gene1"),as.character(g$"Gene2"))
  } else {
    g<-as.character(corem.table[which(corem.table[,3]==regulonID),c("Gene1","Gene2")])
  }
  return(g) 
}

getComembers <- function(geneName = "VNG0700G", corem.table = corem) {
  g <- getCorems(geneName,corem.table)
  o <- sapply(g,function(i){getGenes(i,corem.table)})
  return(o)
}

compareCorems <- function(regulon1=NULL,regulon2=NULL,corem.table,byGene=T,runall=F,p.val=F,jaccard=T) {
  if (byGene) {
    # Get genes
    g1 <- getGenes(regulon1,corem.table)
    g2 <- getGenes(regulon2,corem.table)
    # Calc jaccard similarity
    # o <- length(intersect(g1,g2))/length(unique(c(g1,g2)))
    # Using min -- to correct for different gene set sizes
    if (p.val) {
      o <- phyper(length(intersect(g1,g2)),length(g2),2400-length(g2),length(g1),F)
    } else if (jaccard) {
      o <- length(intersect(g1,g2))/length(unique(c(g1,g2)))
      } else {
      o <- length(intersect(g1,g2))/min(length(g1),length(g2))
    }
    return(o)
  } 
}

makeSubRegulonTable <- function(corem,corem.table) {
  if (class(corem.table)[1]=="data.table") {
    o <- corem.table[corem,mult="all"]
  } else {
    index <- lapply(corem,function(i){which(corem.table[,3]==i)})
    o <- corem.table[as.integer(unlist(index)),]
  }
  return(o)
}

compareAllCorems <- function(corem.table,p.val=T) {
  require(multicore)
  require(data.table)
  print("Comparing all corems")
  if (class(corem.table)[1]=="data.table") {
    setkey(corem.table,Community.ID)
    r.names <- unique(as.character(corem.table[,Community.ID]))
  } else {
    r.names <- unique(corem.table[,3])
  }
  corem.mat <- matrix(0,nrow=length(r.names),ncol=length(r.names),dimnames=list(r.names,r.names))
  arr.i <- which(upper.tri(corem.mat),arr.ind =T)
  arr.i2 <- which(upper.tri(corem.mat))
  val <- unlist(mclapply(seq(1,dim(arr.i)[1]), function(i){
    if(i%%100000==0){print(i)}
    i <- compareCorems(r.names[arr.i[i,1]],r.names[arr.i[i,2]],corem.table,p.val=p.val)
  }))
  corem.mat[arr.i2] <- val
  corem.mat[lower.tri(corem.mat)] <- t(corem.mat)[lower.tri(corem.mat)]
  diag(corem.mat) <- 1
  return(corem.mat)
}

clusterCorems <- function(corem.table, corem.mat=NULL,method=c("merge","cluster")[1],pval=NULL,upper.tail=F) {
  # Get gene overlap
  if (is.null(corem.mat)) {
    print("Computing regulon similarity from scratch")
    corem.mat <- compareAllCorems(corem.table)
  }
  if (method == "cluster") {
    print("Clustering corems")
    # Clusters highly similar corems
    # Will reduce total number of corems
    require(dynamicTreeCut)
    corems.hclust <- hclust(as.dist(1-corem.mat))
    corems.hclust.cut <- cutreeDynamic(corems.hclust)
    o <- list()
    o$corems <- unique(corems.hclust.cut)
    o$genes <- lapply(o$corems,function(i){
     index <- corems.hclust$labels[which(corems.hclust.cut==i)]
     tmp.genes <- lapply(index,function(j){getGenes(j,corem.table)})
     tmp.genes <- unique(unlist(tmp.genes))
     i <- tmp.genes
    })
    names(o$genes) <- o$corems
  } else if (method == "merge") {
    # Merge corems with pval less than bonferonni corrected pval
    # Get sub corems.mat matrix given corem.table
    r <- as.character(unique(corem.table[,Community.ID]))
    corem.mat <- corem.mat[r,r]
    print("Merging corems")
    if (is.null(pval)) {
      pval <- 0.05/dim(corem.mat)[1]
    } else {
      pval <- pval
    }
    if (upper.tail) {
      index <- which(corem.mat>=pval,arr.ind=T)
    } else {
      index <- which(corem.mat<=pval,arr.ind=T)
    }
    o <- list()
    tmp.m <- vector(mode="integer",length=length(unique(as.vector(index))))
    names(tmp.m) <- unique(as.vector(index))
    int <- 1
    for (i in 1:dim(index)[1]) {
      print(paste(i,"out of",dim(index)[1]))
      if ( tmp.m[as.character(index[i,1])] == 0 && tmp.m[as.character(index[i,2])] == 0 ) {
        tmp.m[as.character(index[i,1])] = int
        tmp.m[as.character(index[i,2])] = int
        int = int + 1
      } else if ( tmp.m[as.character(index[i,1])] != 0 && tmp.m[as.character(index[i,2])] == 0 ) {
        tmp.m[as.character(index[i,2])] = tmp.m[as.character(index[i,1])]
      } else if ( tmp.m[as.character(index[i,1])] == 0 && tmp.m[as.character(index[i,2])] != 0 ) {
        tmp.m[as.character(index[i,1])] = tmp.m[as.character(index[i,2])]
      }
    }
    names(tmp.m) <- unique(rownames(corem.mat)[as.vector(index)])
    o$corems <- lapply(sort(unique(tmp.m)),function(i){i <- names(tmp.m)[which(tmp.m==i)]})
    names(o$corems) <- seq(1,length(o$corems))
    o$genes <- lapply(o$corems,function(i){getGenes(i,corem.table)})
  } 
  
  return(o)
}

clusterCorems.2 <- function(geneList) {
  # Assumes that give gene list
  # if there is overlap in the geneList
  # these genes should be merged into
  # a single cluster
  # Works if you've already filtered for 
  # conditional activity
  # Returns a length n list
  o <- vector("integer",length(unique(unlist(geneList))))
  names(o) <- unique(unlist(geneList))
  int <- 1
  pool <- c()
  for (i in 1:length(geneList)) {
    vals <- unique(o[geneList[[i]]])
    if (sum(vals > 0) > 0) {
      # One or more of the genes has already been assigned
      vals <- vals[vals>0]
      if (length(vals) == 1) {
        cID <- vals
        o[geneList[[i]]] = cID
      } else {
        cID <- min(vals)
        o[geneList[[i]]] = cID
        # Reassign other values as well
        for (m in vals) {
          o[which(o==m)] = cID
        }
        pool <- sort(unique(append(pool,vals[-which(vals==cID)])))
      }
    } else {
      if (length(pool) > 0) {
        int <- sort(pool)[1]
        pool <- pool[-1]
      }
      o[geneList[[i]]] = int
      if (length(pool) > 0) {
        int <- sort(pool)[1]
      } else {
        int = int+1
      }
    }
  }
  return(o)
}

getRegulonHits <- function(geneList,corem.table) {
  corems <- lapply(geneList,function(i){getCorems(i,corem.table)})
  # tabulate corems
  corems.t <- table(unlist(corems))
  corems.size <- unlist(lapply(names(corems.t),function(i){
    g <- getGenes(i,corem.table)
    tmp.grep <- grep("VNG*",g)
    if (length(grep("_",g[tmp.grep]))>0) {
      tmp.grep <-tmp.grep[-grep("_",g[tmp.grep])]
    }
    i <- length(tmp.grep)
                          }))
  names(corems.size) <- names(corems.t)
  out <- sort(corems.t[names(corems.t)]/corems.size[names(corems.t)],decreasing=T)
}

getMI <- function(genes = c("VNG2126C", "VNG2624G", "VNG0700G"),MI.index = egrin2.mi, n = 100) {
  o <- list()
  o[["mi"]]<-egrin2.mi[genes,genes]
  o[["mi.mean"]] <- mean(o$mi[upper.tri(o$mi)])
  p.val <- sapply(seq(1,n),function(i){
    i<-sample(rownames(MI.index),length(genes))
    mi<-egrin2.mi[i,i]
    i<- mean(mi[upper.tri(mi)])
    })
  o[["p.val"]] <- sum(o$mi.mean<=p.val)/n
  o[["resamples"]] <- n
  return(o) 
}

getCor <- function(genes = c("VNG2126C", "VNG2624G", "VNG0700G"),
                   exp.ratios = ratios, n = 100, conditions = c("none","condition1,etc")[1],
                   plot = F) {
  o <- list()
  if (length(conditions)<=1) {
    o[["cor"]]<-cor(t(ratios[genes,]))
  } else {
    o[["cor"]]<-cor(t(ratios[genes,conditions]))
    o[["conditions"]]<-conditions
  }
  o[["cor.mean"]] <- mean(o$cor[upper.tri(o$cor)])
  p.val <- sapply(seq(1,n),function(i){
    i<-sample(rownames(ratios),length(genes))
    if (length(conditions)<=1) {
      mi<-cor(t(ratios[i,]))
    } else {
      mi<-cor(t(ratios[i,conditions]))
    }
    i<- mean(mi[upper.tri(mi)])
    })
  o[["p.val"]] <- sum(o$cor.mean<=p.val)/n
  o[["resamples"]] <- n
  return(o) 
}

  
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

getMotifs <- function(genes=NULL,regulon=NULL,corem.table,cutoff=NULL) {
  # Requires egrin2 to be loaded
  if (!is.null(regulon)) {
    g<-getGenes(regulon,regulon.table)
  } else if (!is.null(genes)) {
    g<-genes
  } else {
    cat("you need to provide either a vector of genes or a regulon number and regulon table\n")
    return(invisible(NULL))
  }
  out<-list()
  out$motifs <- lapply(g,function(i){
    #print(i)
    i <- try(egrin2.agglom(src=i,srcType="gene",targetType="motif.cluster",path="motif"),silent=T)
    if (class(i)=="try-error") {
      return(i)
    } else {
      if (!is.null(cutoff)) {
        i<- i[i[,2]<=cutoff,]
      } else {
       i <- i[i[,3]<=0.05,]
     }
     return(i) 
    }
  })
  # Remove try-errors 
  toRemove <- lapply(out$motifs,function(i){class(i)})
  toRemove <- which(unlist(toRemove)=="try-error")
  if (length(toRemove)>0) {
    out$motifs <- out$motifs[-toRemove]
    g <- g[-toRemove]
  }
  # Remove motifs above MOTC_301
  out$motifs <- lapply(seq(1:length(out$motifs)),function(i){
    n <- rownames(out$motifs[[i]])
    n.split <- sapply(n,function(j){strsplit(j,split="_")[[1]][2]})
    i <- out$motifs[[i]][!as.integer(n.split)>301,]
  })
  # Remove genes without motifs
  toRemove <- lapply(seq(1:length(out$motifs)),function(i){dim(out$motifs[[i]])[1]})
  toRemove <- which(unlist(toRemove)==0)
  if (length(toRemove)>0) {
    out$motifs <- out$motifs[-toRemove]
    g <- g[-toRemove]
  }
  if (length(out$motifs)==0) {
    # no motifs
    return(NULL)
  }
  #names(out$motifs2) <- g2
  all <- c()
  mots <- list()
  for (i in 1:length(out$motifs)) {
    mots[[i]]<-out$motifs[[i]][,1]
    names(mots[[i]]) <- rownames(out$motifs[[i]])
    # Remove bad clusters
    toRemove <- which(names(mots[[i]])%in%bad.clusts==TRUE)
    if (length(toRemove)>0) {
      mots[[i]]<-mots[[i]][-toRemove]
    }
    all <- append(all,names(mots[[i]]))
  }
  names(mots) <- g
  out$motifs <- mots
  all.n <- unique(all)
  all <- vector(mode = "integer", length=length(all.n))
  names(all) <- all.n 
  for (i in 1:length(mots)) {
    all[names(mots[[i]])] <- all[names(mots[[i]])] + mots[[i]] 
  }
  out$all.motifs <- sort(all,decreasing=T)
  # Combine motifs into common score
  mot.score <- matrix(nrow=length(g),ncol=length(g),dimnames=list(g,g))
  for (i in g) {
    for (j in g) {
      # jaccard
      mot.score[i,j] <- length(intersect(names(mots[[i]]),names(mots[[j]])))/length(union(names(mots[[i]]),names(mots[[j]])))
    }
  }
  out$pairwise.motif.score <- mot.score
  # Moitf Gene Counts
  mot.score2 <- matrix(0,nrow=length(out$all.motifs),ncol=length(g),
                       dimnames=list(names(out$all.motifs),g))
  for (i in g) {
    mot.score2[names(mots[[i]]),i] <- mots[[i]]/sum(mots[[i]])
  }
  out$motif.score <- mot.score2
  return(out)
}

getConditions <- function(genes=NULL,regulon=NULL,corem.table,cutoff=NULL,pval=0.5,
                          enforce.diff=F,diff.cutoff=2,ratios=NULL,ratios.normalized=F) {
  # Requires egrin2 to be loaded
  if (!is.null(regulon)) {
    g<-getGenes(regulon,regulon.table)
  } else if (!is.null(genes)) {
    g<-genes
  } else {
    cat("you need to provide either a vector of genes or a regulon number and regulon table\n")
    return(invisible(NULL))
  }
  out <- egrin2.agglom(src=g,srcType="gene",targetType="condition",path="bicluster")
  if (!is.null(cutoff)) {
    out <- out[out[,2]<=cutoff,]
  } else {
    out <- out[out[,3]<=0.05,]
  }
#   conds.score <- matrix(0,nrow=length(out$conditions),ncol=length(g),
#                        dimnames=list(names(out$conds),rownames(out$pairwise.conds.score)))
#   for (i in rownames(out$conds.score)) {
#     conds.score[rownames(conds.i[[i]]),i] <- (conds.i[[i]][,1])/sum(conds.i[[i]][,1])
#   }
#   out$conds.score <- conds.score
  out <- out[order(out[,3]),]
  if (enforce.diff) {
    if (is.null(ratios)) {
      print("Please provide ratios.")
      return(NULL)
    }
    normRatios <- function(ratios) {
      o <- t(scale(t(ratios), 
                center = apply(ratios, 1, median, 
                  na.rm = T), scale = apply(ratios, 
                  1, sd, na.rm = T)))
      return(o)
    }
    # normalize ratios (just in case they're not)
    if (ratios.normalized == F) {
      ratios <- normRatios(ratios)
    }
    meanExp <- abs(colMeans(ratios[g,rownames(out)]))
    meanExp <- which(meanExp>=diff.cutoff)
    out <- out[meanExp,]
  }
  return(out)
}

compareConditions <- function(regulonConditionList,method=c("jaccard","cosine")[1]) {
  out<-list()
  if (method == "cosine") {
    require(lsa)
    index <- colnames(ratios)
    pairwise.conds.score <- matrix(nrow=length(index),ncol=length(regulonConditionList),
                                 dimnames=list(index,names(regulonConditionList)))
   for (i in 1:length(regulonConditionList)) {
      tmp <- regulonConditionList[[i]][,1]
      names(tmp)<-(regulonConditionList[[i]])
      pairwise.conds.score[index,names(regulonConditionList)[i]] <- tmp[index]
   }
   pairwise.conds.score[is.na(pairwise.conds.score)] <- 0
   pairwise.conds.score <- cosine(pairwise.conds.score)
  } else if (method == "jaccard") {
    pairwise.conds.score <- matrix(nrow=length(regulonConditionList),
                                   ncol=length(regulonConditionList),
                                   dimnames=list(names(regulonConditionList),
                                                 names(regulonConditionList)))
    for (i in 1:length(regulonConditionList)) {
      for (j in 1:length(regulonConditionList)) {
        pairwise.conds.score[i,j] <- length(intersect((regulonConditionList[[i]]),
                                            (regulonConditionList[[j]])))/length(union((regulonConditionList[[i]]),
                                            (regulonConditionList[[j]])))
      }
   }
  }
  return(pairwise.conds.score)
}  

compareMotifs <- function(regulonMotifList,method=c("jaccard","cosine")[1]) {
  out<-list()
  if (method == "cosine") {
    require(lsa)
    index <- c()
    for (i in 1:length(regulonMotifList)) {
         index <- append(index,names(regulonMotifList[[i]]))
    }
    index <- unique(index)
    pairwise.motif.score <- matrix(nrow=length(index),ncol=length(regulonMotifList),
                                 dimnames=list(index,names(regulonMotifList)))
   for (i in 1:length(regulonMotifList)) {
      tmp <- regulonMotifList[[i]]
      pairwise.motif.score[index,names(regulonMotifList)[i]] <- tmp[index]
   }
   pairwise.motif.score[is.na(pairwise.motif.score)] <- 0
   pairwise.motif.score <- cosine(pairwise.motif.score)
  } else if (method == "jaccard") {
    pairwise.motif.score <- matrix(nrow=length(regulonMotifList),
                                   ncol=length(regulonMotifList),
                                   dimnames=list(names(regulonMotifList),
                                                 names(regulonMotifList)))
    for (i in 1:length(regulonMotifList)) {
      for (j in 1:length(regulonMotifList)) {
        pairwise.motif.score[i,j] <- length(intersect(names(regulonMotifList[[i]]),
                                            names(regulonMotifList[[j]])))/length(union(names(regulonMotifList[[i]]),
                                            names(regulonMotifList[[j]])))
      }
   }
  }
  return(pairwise.motif.score)
}

compareMotifs.gene <- function(genes,pval=1e-6,count.cutoff=5,alt.annotation=F,...) {
  # assumes egrin2 env loaded
  require(multicore)
  g.m <- mclapply(genes,function(i) {
    to.r <- try(out$plot.promoter.architecture(i,dont.plot=T,verbose=F,p.value=pval,...))
    #to.r <- try(out$plot.promoter.architecture(i,dont.plot=T,verbose=F,p.value=pval))
    if (class(to.r)=="try-error") {
      return(character())
    } else {
      to.r <- to.r$mot.tab[to.r$mot.tab>=count.cutoff]
      return(to.r)
    }
  })
  names(g.m) <- genes
  # currently doesn't integrate counts - can think of a way to integrate in future
  if (alt.annotation) {
    # return annotation data frame
    u.m <- unique(unlist(sapply(g.m,function(i)names(i))))
    ord <- order(as.numeric(sapply(u.m,function(i)strsplit(i,split="_")[[1]][2])),decreasing=F)
    u.m <- u.m[ord]
    motif.annot <- as.data.frame(do.call(rbind,mclapply(genes,function(i){
      to.r <- rep(0,length(u.m))
      to.r[u.m%in%names(g.m[[i]])]=1
      names(to.r)<-u.m
      return(to.r)
        })))
    rownames(motif.annot) <- genes
    return(motif.annot)
  } else {
    pairwise.motif.score <- do.call(cbind,mclapply(genes,function(i){
      sapply(genes,function(j){
        j<-length(intersect(names(g.m[[i]]),names(g.m[[j]])))/length(names(g.m[[i]]))
        if (is.nan(j)) {
          return(0)
        } else {
          return(j)
        }
      })
    }))
    colnames(pairwise.motif.score) <- genes
    return(pairwise.motif.score)
  }
}

combineMotifScore <- function(regulonMotifScoreList) {
  all.motifs <- unique(unlist(sapply(regulonMotifScoreList,rownames)))
  all.genes <- unique(unlist(sapply(regulonMotifScoreList,colnames)))
  out <- matrix(data=0,nrow=length(all.motifs),ncol=length(all.genes),dimnames=list(all.motifs,all.genes))
  for (i in 1:length(regulonMotifScoreList)) {
    out[rownames(regulonMotifScoreList[[i]]),colnames(regulonMotifScoreList[[i]])] <- regulonMotifScoreList[[i]]
  }
  #renormalize
  out <- apply(out,2,function(i){i/sum(i)})
  return(out)
}
  
getResid <- function (rowNames,colNames, rats,rats.inds = "COMBINED", varNorm = F, in.cols = T, 
                      row.score.func="orig") {
    residual.submatrix <- function(rats, rows, cols, varNorm = F, 
        ...) {
        rows <- rows[rows %in% rownames(rats)]
        cols <- cols[cols %in% colnames(rats)]
        if (length(rows) <= 1 || length(cols) <= 1) 
            return(1)
        rats <- rats[rows, cols]
        if (is.vector(rats) || any(dim(rats) <= 1) || mean(is.na(rats)) > 
            0.95) 
            return(1)
        d.rows <- rowMeans(rats, na.rm = T)
        d.cols <- colMeans(rats, na.rm = T)
        d.all <- mean(d.rows, na.rm = T)
        rij <- rats + d.all
        rij[, ] <- rij[, ] - matrix(d.cols, nrow = nrow(rij), 
            ncol = ncol(rij), byrow = T)
        rij[, ] <- rij[, ] - matrix(d.rows, nrow = nrow(rij), 
            ncol = ncol(rij), byrow = F)
        average.r <- mean(abs(rij), na.rm = TRUE)
        if (varNorm) {
            maxRowVar <- attr(rats, "maxRowVar")
            row.var <- mean(apply(rats, 1, var, use = "pairwise.complete.obs"), 
                na.rm = T)
            if (is.na(row.var) || row.var > maxRowVar) 
                row.var <- maxRowVar
            average.r <- average.r/row.var
        }
        average.r
    }
    inds <- rats.inds
    if (rats.inds[1] == "COMBINED") 
        inds <- names(get("row.weights",envir=e))
    resids <- sapply(e$ratios[inds], function(rn) {
        if (row.score.func == "orig") {
            if (in.cols) 
                residual.submatrix(rn, rowNames, colNames, 
                  varNorm = varNorm)
            else residual.submatrix(rn, rowNames, colnames(rn)[!colnames(rn) %in% 
                colNames], varNorm = varNorm)
        }
        else {
            print("This is not supported at the moment")
#             if (in.cols) 
#                 mean(get.row.scores(k, for.rows = get.rows(k), 
#                   ratios = rn, method = row.score.func))
#             else mean(get.row.scores(k, cols = get.cols(k), for.rows = get.rows(k), 
#                 ratios = rn, method = row.score.func))
        }
    })
    if (rats.inds[1] == "COMBINED") 
        resids <- weighted.mean(resids, e$row.weights[inds], na.rm = T)
    if (rats.inds[1] != "COMBINED" && length(resids) < length(inds) && 
        all(is.na(resids))) {
        resids <- rep(NA, length(inds))
        names(resids) <- inds
    }
    return(resids)
}

plotExpression <- function(x,sort=NULL,mode=c(1)[1],toHighlight=NULL,colHighlight=NULL,...) {
  # Input, x, should be matrix or list of matrices
  if (!exists("sortExp")) {
    source("~/R/sortExp.R")
  }
  if (!exists("plotCI")) {
    source("~/R/plotCI.R")
  }
  # Matrix = gene expression, row = gene, col = condition
  # To highlight is a list of vectors of genes (or just vector)
  if (class(x) == "list") {
    m <- list()
    m$m <- list()
    m$s <- list()
    m$r <- list()
    for (i in 1:length(x)) {
      if (class(x[[i]]) !="matrix") {
        print("you need to supply a matrix or list of matrices")
        return(NULL)
      }
      if (!is.null(sort)) {
        x[[i]] <- sortExp(x[[i]],sort)
      }
      if (mode==1) {
        m$m[[i]] <- colMeans(x[[i]])
        m$s[[i]] <- sd(x[[i]])
        m$r[[i]] <- rep(getResid(rownames(x[[i]]),colnames(x[[i]])),length(m$m[[i]]))
        m$h[[i]] <- x[[i]][toHighlight[[i]],]
        m$ch[[i]] <- colHighlight[[i]]
      }
    }
    plotCI(m$m,m$s,m$r,add.data=m$h,add.colors=m$ch,...)
  } else if (class(x) == "matrix") {
    if (!is.null(sort)) {
      x = sortExp(x,sort)
   }
    if (mode == 1) {
     m <- colMeans(x)
     m.s <- sd(x)
     m.r <- rep(getResid(rownames(x),colnames(x)),length(m))
     m.h <- x[toHighlight,]
     m.ch <- colHighlight
     plotCI(m,m.s,m.r,add.data=m.h,add.colors=m.ch,...)
   }
  }
  #legend("topleft",legend = rownames(toHighlight),col = seq(1:length(toHighlight)),
  #       lty=seq(1,1,length.out=length(toHighlight)),lwd=4)
}

getSigCoexpression <- function(regulonNumbers,refObj=NULL,findRowsList=NULL, n = 20,
                               ratios=ratios,samples=NULL,test=c("residual","mi")[1]) {
  if (is.null(refObj)) {
      print("You must provide ref object containing [...Corems]$conditions and [...Corems]$genes")
      invisible(NULL)
    }
  # If unspecified, run findRows2 to get conditions most sig associated
  # With all combinations of 3 corems
  if (is.null(findRowsList)) {
      findRowsList <- findRows2(regulonConditionList=refObj$conditions[regulonNumbers],method="all",n=n)
  }
  # Make and fill a matrix with p.vals from resampling
  # Row = set of conditions, e.g. r1.r2.r3 is intersection of conditions assigned to all 3 corems
  # Col = set of genes, e.g. g1.g2.g3 is union of genes from all three corems
   # Set # resamples
  if (is.null(samples)) {
    samples <- (1/(0.05/(length(findRowsList)^2)))
  }
  m <- matrix(0,nrow=length(findRowsList),ncol=length(findRowsList),
              dimnames=list(names(findRowsList),names(findRowsList)))
  for (i in 1:length(findRowsList)) {
    conds <- findRowsList[[i]]
    if (test == "mi") {
      if (!exists("aracneOnTheFly")) {
          require("rARACNE.R")
        }
      rats.sub <- ratios[,conds]
      mi.m <- aracneOnTheFly(dataMatrix=rats.sub,clean=T)
    }
    #tmpL <- list()
    for (j in 1:length(findRowsList)) {
      genes <- unique(unlist(lapply(unlist(strsplit(names(findRowsList)[[j]],split=".",fixed=T)),
                                    function(x){x<-refObj$genes[[x]]})))
      if (test=="residual") {
        if (!exists("e")) {
          print("Please load EGRIN env, e")
        }
        src <- getResid(rowNames=genes,colNames=conds,rats=ratios)
        comp <- c()
        for (s in 1:samples) {
          gs <- sample(rownames(ratios),size=length(genes))
          comp <- append(comp,getResid(rowNames=gs,colNames=conds,rats=ratios))
        }
        m[i,j] <- sum(src>=comp)/length(comp)
        m2 <- p.adjust(m,method="bonferroni")
        m3<-matrix(m2,nrow=dim(m)[1],ncol=dim(m)[2],dimnames=list(rownames(m),colnames(m)))
      } else if (test == "mi") {
        src <- mean(mi.m$Mutual.Information[genes,genes])
        comp <- c()
        for (s in 1:samples) {
          gs <- sample(rownames(ratios),size=length(genes))
          comp <- append(comp,mean(mi.m$Mutual.Information[gs,gs]))
        }
        m[i,j] <- sum(src<=comp)/length(comp)
        m2 <- p.adjust(m,method="bonferroni")
        m3<-matrix(m2,nrow=dim(m)[1],ncol=dim(m)[2],dimnames=list(rownames(m),colnames(m)))
        #tmpL[[names(findRowsList)[j]]]$src <- src
        #tmpL[[names(findRowsList)[j]]]$comp <- comp
      }
    }
  }
 return(m3) 
 #return(tmpL)
}

getSigCoexpression.2 <- function(regulon,ratios,gBg_list = gBg_backbone_0.59_clean_list,
                                 mode=c("cvar","egrin2")[1],resamples=14340) {
  # 14340 is number of resamples to hit bonferroni corrected pval given 717 corems
  randGene <- function(size.g) {
    o <- sample(rownames(ratios),size.g)
  }
  n.genes <- length(gBg_list$genes[[regulon]])
  if (mode == "cvar") {
    # Compute mean cor of this regulon
    ref <- mean(cor(t(ratios[gBg_list$genes[[regulon]],names(gBg_list$conditions.cvar[[regulon]])])))
    # Calc # of conditions
    n.conds <- length(names(gBg_list$conditions.cvar[[regulon]]))
    o <- unlist(mclapply(seq(1:resamples),function(i){
      # Compute cvar
      g <- randGene(n.genes)
      # Choose n.conds conditions with lowest cvar for these genes
      best.conds <- (sort(apply(ratios[g,],2,sd)/abs(colMeans(ratios[g,]))))[1:n.conds]
      # Compute correlation
      o <- mean(cor(t(ratios[g,names(best.conds)])))
      return(o)
    }))
  } else if (mode == "egrin2") {
    # Compute mean cor of this regulon
    ref <- mean(cor(t(ratios[gBg_list$genes[[regulon]],rownames(gBg_list$conditions.egrin2[[regulon]])])))
    # Calc # of conditions
    n.conds <- length(rownames(gBg_list$conditions.egrin2))
    o <- unlist(mclapply(seq(1:resamples),function(i){
      # Compute egrin2 score
      g <- randGene(n.genes)
      # Choose n.conds conditions with lowest egrin2 pval for these genes
      conds <- egrin2.agglom(g,'gene','condition','bicluster')
      best.conds <- conds[order(conds[,2]),][1:n.conds,]
      # Compute correlation
      o <- mean(cor(t(ratios[g,rownames(best.conds)])))
      return(o)
    }))
  }
  out <- sum(ref<o)/length(o)
  return(out)
}

getMemeHitLocations <- function(regulon,refObj,distance=c(-20,500)) {
  require(stringr)
  # Must work in egrin local environment
  wd <- getwd()
  #setwd("/proj/baliga/dreiss/egrin2_v0.4/")
  genes <- refObj$genes[[regulon]]
  gene.sequences <- e$get.sequences(genes,distance=distance)
  starts.stops <- attr(gene.sequences,"start.stops")
  motif.clusters <- names(refObj$motifs[[regulon]]$all.motifs)
  motifs <- lapply(motif.clusters,function(i){unlist(get.motifs(motif.clusters=i))})
  names(motifs) <- motif.clusters
  motif.cluster.info <- lapply(motifs,function(i){
    out<-lapply(i,function(j){
      i<-get.motif.info(j)
    })
    names(out) <- i
    return(out)
  })
  # Remove genes with no sequence
  toRemove <- which(!genes%in%rownames(starts.stops))
  if (length(toRemove > 0)) {
    genes <- genes[-toRemove]
  }
  o <- list()
  for (g in genes) {
    print(g)
    o[[g]] <- list()
    s <- starts.stops[g,"start"]
    es <- starts.stops[g,"end"]
    strand <- starts.stops[g,"strand"]
    chr <- starts.stops[g,"contig"]
    if (strand == "D") {
      Gseq <- seq(s-1,es,1)
    } else if (strand == "R") {
      Gseq <- seq(es-1,s,-1)
    }
    o[[g]]$counts <- matrix(0,nrow=length(motif.clusters),ncol=length(Gseq),
                         dimnames=list(motif.clusters,Gseq))
    o[[g]]$strands <- list()
    o[[g]]$consensus <- list()
    for (mc in names(motif.cluster.info)) {
      motifs.to.consider <- lapply(motif.cluster.info[[mc]],
                                   function(i){i<-grep(g,
                                                    i[[1]]$posns[,"gene"])})
      motifs.to.consider <- names(unlist(motifs.to.consider))
      for (m in motifs.to.consider) {
        g.m<-which(motif.cluster.info[[mc]][[m]][[1]]$posns[,"gene"]==g)                           
        s.match <- as.character(motif.cluster.info[[mc]][[m]][[1]]$posns[g.m,"site"])
        pm <- motif.cluster.info[[mc]][[m]][[1]]$posns[g.m,"strand"]
        if (pm == "+") {
          pm <- 1
        } else if (pm == "-") {
          pm <- -1
        }
        o[[g]]$strands[[mc]] <- append(o[[g]]$strands[[mc]],pm)
        if (motif.cluster.info[[mc]][[m]][[1]]$posns[g.m,"strand"]=="-") {
          s.match <- e$rev.comp(s.match)
        }
        loc <- str_locate(gene.sequences[g],s.match)
        if (sum(is.na(loc))==0) {
          seq.loc <- seq(loc[1],loc[2],1)
          o[[g]]$counts[mc,seq.loc] = o[[g]]$counts[mc,seq.loc]+1
        }
      }
      if (sum(o[[g]]$counts[mc,])==0) {
        o[[g]]$consensus[[mc]] <- NA
      } else {
      # Find Peak
      peak <- names(which(o[[g]]$counts[mc,]==max(o[[g]]$counts[mc,]))[1])
      # Find left boundary
      left <- peak
      if (which(names(o[[g]]$counts[mc,])==left) != 1) {
        # This is the left boundary
        left <- peak
        while ( (which(names(o[[g]]$counts[mc,])==left) != 1) && (o[[g]]$counts[mc,left]>0) ) {
          if (strand == "R") {
            left <- as.character(as.integer(left)+1)
          } else {
            left <- as.character(as.integer(left)-1)
          }
        }
        if (strand == "R") {
            left <- as.character(as.integer(left)-1)
          } else {
            left <- as.character(as.integer(left)+1)
          }
      }
      right <- peak
      if (which(names(o[[g]]$counts[mc,])==right) != length(o[[g]]$counts[mc,])) {
        # This is the right boundary
        right <- peak
        while ( (which(names(o[[g]]$counts[mc,])==right) != length(o[[g]]$counts[mc,])) && (o[[g]]$counts[mc,right]>0) ) {
          if (strand == "R") {
            right <- as.character(as.integer(right)-1)
          } else {
            right <- as.character(as.integer(right)+1)
          }
        }
        if (strand == "R") {
            right <- as.character(as.integer(right)+1)
          } else {
            right <- as.character(as.integer(right)-1)
          }
      }
      
      if (strand == "R") {
        l.b <- left
        r.b <- right
        left <- r.b
        right <- l.b
      }
      c.seq <- as.character(substr(e$genome.info$genome.seqs[as.character(chr)],left,right))
      if (strand == "R") {
        c.seq <- e$rev.comp(c.seq)
      }
#       # Vote on the strand. Pull out the consensus sequence.
#       consensus <- mean(o[[g]]$strands[[mc]])
#       if (consensus == 0) {
#         consensus <- sample(x=c(-1,1),size=1)
#       }
#       if (consensus > 0) {
#         o[[g]]$consensus[[mc]] <- c.seq
#       } else if (consensus < 0) {
#         o[[g]]$consensus[[mc]] <- e$rev.comp(c.seq)
#       } 
      o[[g]]$consensus[[mc]] <- c.seq
    }               
  }
}
setwd(wd)
return(o)
}

plot.promoter.architecture.seq <- function( gene, n.motifs=10, return.em=T, window=c(200,20), plot=T,plot.all=T, include.bad=F, meme=F,... ) {
  require(matlab)
  require(gplots)
  coo <- e$get.gene.coords( gene, op.shift=T )
  if ( is.null( coo ) ) return( NULL )
  st.st <- c( coo$start_pos, coo$end_pos )
  if ( coo$strand == "R" ) {
    st.st <- rev( st.st )
    st.st <- st.st[ 1 ] + c( -(window[2]-1), +window[1] )
  } else {
    st.st <- st.st[ 1 ] + c( -window[1], +(window[2]-1) )
  }
  chr=as.character( coo$contig )
  names( st.st )[ 1 ] <- chr
  ##tmp1=get.biclusters(gene=gene)[[1]]
  ##tmp1a=unlist(get.motifs(bic=tmp1))
  tmp1a <- unique( unlist( get.motifs( pos=st.st ) ) )
  if ( length( tmp1a ) <= 0 ) return( NULL )
  if ( ! include.bad && exists( "bad.clusts" ) )
    tmp1a <- tmp1a[ ! tmp1a %in% unlist( get.motifs( motif.clust=bad.clusts, expand=F ) ) ]
  if ( length( tmp1a ) <= 0 ) return( NULL )
  tmpa1 <- get.motif.positions( tmp1a, start=st.st[1], stop=st.st[2], chr=as.character( coo$contig ),
                              return=T, plot=F )
  ##if ( is.null( tmpa ) ) return( NULL )
  out <- list()
  if ( return.em ) out[[ 'ALL' ]] <- tmpa1
  tmp2 <- lapply( motif.clusts[ 1:mc.length ], function( i ) i[ i %in% tmp1a ] )
  names( tmp2 ) <- paste( "MOTC", 1:length( tmp2 ), sep="_" )
  ##tmp2=unlist(get.motif.clusters(mot=tmp1a))
  tmp2a <- sort( sapply( tmp2, length ), decreasing=T ) ##table(tmp2),decreasing=T)
  if ( ! include.bad && exists( "bad.clusts" ) )
    tmp2a <- tmp2a[ ! names( tmp2a ) %in% bad.clusts ] ##paste( "MOTC", bad.clusts, sep="_" ) ]
  ##tmp2a <- sapply( motif.clusts, function( i ) mean( i %in% tmp1a ) ); names( tmp2a ) <- names( tmp2 )
  ##tmp2a <- sort( tmp2a, decreasing=T )
  ##names( tmp2a ) <- names( tmp2 )
  print( tmp2a[ 1:n.motifs ] )
  if ( length( tmp2a ) <= 0 || all( tmp2a <= 0 ) ) return( NULL )
  if (meme) n.motifs <- 1000
  if ( length( tmp2a ) < n.motifs ) n.motifs <- length( tmp2a )
  for ( i in 1:n.motifs ) {
    mcs <- names( tmp2a )[ i ]
    ## if ( exists( "motif.cluster.clusters" ) ) {
    ##   mccs <- unique( motif.cluster.clusters[ as.integer( gsub( "MOTC_", "", mcs ) ) ] )
    ##   mccs <- mccs[ mccs != 0 ]
    ##   if ( length( mccs ) > 0 ) mcs <- paste( "MOTC", which( motif.cluster.clusters %in% mccs ), sep="_" )
    ##   if ( ! include.bad && exists( "bad.clusts" ) ) mcs <- mcs[ ! mcs %in% bad.clusts ]
    ## }
    cat(mcs,length(unlist( get.motifs( motif.clusters=mcs ) )),"\n")
    tmp <- get.motif.positions( unlist( get.motifs( motif.clusters=mcs ) ), start=st.st[1], stop=st.st[2],
                               chr=chr ) ##, plot=i==1, ylim=range( tmpa1, na.rm=T ) )
    if ( return.em ) out[[ names(tmp2a)[i] ]] <- tmp
  }
  if (meme) {
    # restrict motif clusters to those found by meme (mast hits is default)
    m <- getMotifs(gene)
    out <- out[c("ALL",names(m$all.motifs))]
    out[["ALL"]][] = 0 
    for (i in 2:length(out)) {
      out[["ALL"]] = out[["ALL"]]+out[[i]]
    }
  }
  o <- c()
  for (i in 1:length(out)){
    print(i)
    o<-rbind(o,out[[i]])
  }
  rownames(o) <- names(out)
  if ( coo$strand == "R" ) {
    o <- fliplr(o)
  }
  out$counts <- o
  out$sequence <- e$get.sequences(gene,distance=c(-window[2],window[1]))
  if (plot) {
    if (plot.all) {
      plotGeneMotifNucleotideSequence(out$sequence,out$counts,seqName=names(out$sequence))
    } else {
      plotGeneMotifNucleotideSequence(out$sequence,out$counts[2:dim(out$counts)[1],],seqName=names(out$sequence))
    }
  }
  invisible( out )
}


getSig <- function(regulon=NULL,genes=NULL,corem.table=NULL,assoc=c("condition","motif")[1],n=1000,alt=T,pval=.05) {
  # If regulon is input get genes in regulon. Requires regulon lookup table
  require(multicore)
  if (!is.null(regulon)) {
    if (is.null(corem.table)) {
      print("Must specific regulon table for lookup")
      return(NULL)
    } else {
      genes <- getGenes(regulon,regulon.table)
    }
  }
    # Get hits to conditions or motifs
    if (assoc=="condition") {
      o <- lapply(genes,function(i){
        o.tmp<-egrin2.agglom(src=i,srcType="gene",targetType="condition",path="bicluster")
        o <- o.tmp[,3][o.tmp[,3]<=pval]
        names(o) <- rownames(o.tmp)[o.tmp[,3]<=pval]
        return(o)
        })
    } else if (assoc=="motif") {
      if (alt) {
        o <- lapply(genes,function(i){
          out.tmp <- try(egrin2.agglom(src=i,srcType="gene",targetType="motif.cluster",path="motif"))
          if (class(out.tmp) == "try-error") {
            return(NULL)
          } else {
          out <- out.tmp[,1]
          names(out) <- rownames(out.tmp)
          try(load("~/R/egrin2/gBg_antoine/motif_map.RData"))
          if (exists("motif.map")) {
            out<-out[names(out)%in%names(motif.map)]
            names(out) <- motif.map[names(out)]
            out <- out[out>5]
            # Should change later to return pval too. That needs some work.
            out <- unique(names(out))
          }
          return(out)
          }
        })
        names(o) <- genes
      } else {
        o <- egrin2.agglom(src=genes,srcType="gene",targetType="motif.cluster",path="motif")
        # Check to see if random sample file for this number of genes has already
        # been computed. If so, load it.
        if (file.exists(paste("~/R/egrin2/gBg_antoine/randomSampleMotifs/",length(genes),".RData",sep=""))) {
          load(paste("~/R/egrin2/gBg_antoine/randomSampleMotifs/",length(genes),".RData",sep=""))
        } else {
          #Compute random samples
          r <- mclapply(seq(1:n),function(i){
            if (i%%100) {
              print(i)
            }
            i <- egrin2.agglom(src=sample(rownames(e$row.membership),size=length(genes)),
                               srcType="gene",targetType="motif.cluster",path="motif")
          })
          save(r,file=paste("~/R/egrin2/gBg_antoine/randomSampleMotifs/",length(genes),".RData",sep=""))
        }
        r.vals <- mclapply(seq(1:dim(o)[1]),function(i){
          out <- sapply(seq(1:length(r)),function(j){j<-r[[j]][rownames(o)[i],1]})
          out[is.na(out)] <- 0
          i <- sum(out>=o[i,1])/length(out)
        })
        r.vals <- unlist(r.vals)
        r.vals <- p.adjust(r.vals,method="BH")
        r.vals <- sort(r.vals,decreasing=F)
        names(r.vals) <- rownames(o)
        o <- r.vals
        try(load("~/R/egrin2/gBg_antoine/motif_map.RData"))
        if (exists("motif.map")) {
          o<-o[names(o)%in%names(motif.map)]
          names(o) <- motif.map[names(o)]
          # Should change later to return pval too. That needs some work.
          o <- unique(names(o)[o<0.05])
        }
      }
    }
  toRemove <- c()
  for (i in 1:length(o)) {
    if (is.null(o[[i]]) || length(o[[i]]) == 0 ) {
      toRemove <- append(toRemove,i)
    }
  }
  if (length(toRemove)>0) {
    o <- o[-toRemove]
  }
  return(o)  
}

randomConditions <- function(geneSetSize=seq(3,176,1),ratios,resamples=30000,method=c("sd","c.var")[1]) {
  require(multicore)
  genePool <- rownames(ratios)
  geneSetSize <- as.character(geneSetSize)
  if (method == "sd") {
    o<-mclapply(seq(1,length(geneSetSize)),function(i) {
      len = as.integer(geneSetSize[i])
      print(len)
      if (file.exists(paste("./sd/",len,method,resamples,".txt.gz",sep=""))) {
        return(NULL)
      } else {
        i <- matrix(0,nrow=resamples,ncol=dim(ratios)[2],dimnames = list(seq(1,resamples,1),colnames(ratios)))
        for (j in seq(1,resamples,1)) {
          g <- sample(genePool,len)
          i[j,] <- sd(ratios[g,])
        }
        write.table(i,file=paste("./sd/",len,method,resamples,".txt",sep=""))
        system(paste("gzip",paste("./cvar/",len,method,resamples,".txt",sep=""),sep=" "))
        return(NULL)
      }
    })
  } else if (method == "c.var") {
    o<-mclapply(seq(1,length(geneSetSize)),function(i) {
      len = as.integer(geneSetSize[i])
      print(len)
      if (file.exists(paste("./cvar/",len,method,resamples,".txt",sep=""))) {
        return(NULL)
      } else {
        i <- matrix(0,nrow=resamples,ncol=dim(ratios)[2],dimnames = list(seq(1,resamples,1),colnames(ratios)))
        for (j in seq(1,resamples,1)) {
          g <- sample(genePool,len)
          m <- abs(colMeans(ratios[g,]))
          m[which(m==0)] = 1e-6
          i[j,] <- sd(ratios[g,])/m
        }
        write.table(i,file=paste("./cvar/",len,method,resamples,".txt",sep=""))
        system(paste("gzip",paste("./cvar/",len,method,resamples,".txt",sep=""),sep=" "))
        return(NULL)
      }
    })
  }
  invisible(o)
}

getSignificantConditions <- function(genes,ratios,resamples=50000,method=c("sd","c.var")[1],ref=NULL) {
  if (!file.exists(paste("~/R/egrin2/gBg_antoine/conditionResamples/",length(genes),method,resamples,".RData",sep="")) && is.null(ref)) {
    print("Computing values")
    ref <- randomConditions(length(genes),rownames(ratios),ratios,method=method)
  } else if (file.exists(paste("~/R/egrin2/gBg_antoine/conditionResamples/",length(genes),method,resamples,".RData",sep=""))) {
    print("Loading values from file")
    load(paste("~/R/egrin2/gBg_antoine/conditionResamples/",length(genes),method,resamples,".RData",sep=""))
    ref <- o
  } else {
    print("Using provided values")
  }
  if (method == "sd") {
    var <- sd(ratios[genes,])
  } else if (method == "c.var") {
    m <- abs(colMeans(ratios[genes,]))
    m[which(m==0)] = 1e-6
    var <- sd(ratios[genes,])/m
  }
  o<-mclapply(var,function(i){i<-sum(ref<=i)/length(ref)})
  o <- unlist(o)
  names(o) <- names(var)
  o <- o[o<=.05]
  return(o)
}

getSignificantConditions.2 <- function(genes,ratios,ratios.normalized=F,resamples=30000,method=c("sd","c.var")[1],
                                       all=F,pval=0.05,enforce.diff=F,diff.cutoff=2,...) {
  require("bigmemory")
  lookup.table = read.big.matrix(paste("/media/OS_Install/Ubuntu_docs/conditionResamples/cvar/",length(genes),method,resamples,".txt",sep=""),
                                 sep=" ",has.row.names=T,ignore.row.names=T,
                                 skip=1,type="double")
  if (method == "sd") {
    var.m <- matrix(apply(ratios[genes,],2,sd),nrow=resamples,ncol=dim(ratios)[2],dimnames = list(seq(1,resamples,1),colnames(ratios)),byrow=T)
    comp <- var.m>lookup.table
    o <- apply(comp,2,sum)/dim(comp)[1]  
  } else if (method == "c.var") {
    m <- abs(colMeans(ratios[genes,]))
    m[which(m==0)] = 1e-6
    var.m <- matrix(apply(ratios[genes,],2,sd)/m,nrow=resamples,ncol=dim(ratios)[2],dimnames = list(seq(1,resamples,1),colnames(ratios)),byrow=T)
    comp <- var.m>as.matrix(lookup.table)
    o <- apply(comp,2,sum)/dim(comp)[1]
  }
  if (!all) {
    # Only report conds with p<=pval
    o <- o[o<=pval]
  }
  if (enforce.diff) {
    normRatios <- function(ratios) {
      o <- t(scale(t(ratios), 
                center = apply(ratios, 1, median, 
                  na.rm = T), scale = apply(ratios, 
                  1, sd, na.rm = T)))
      return(o)
    }
    # normalize ratios (just in case they're not)
    if (ratios.normalized == F) {
      ratios <- normRatios(ratios)
    }
    meanExp <- abs(colMeans(ratios[genes,names(o)]))
    meanExp <- which(meanExp>=diff.cutoff)
    o <- o[meanExp]
  }
  return(o)
}

# Cluster Motifs

scoringFunction <- function(cut,clustObj,motifObj) {
  # Take pariwise motif simialrities from hclust cut
  # Determines the mean/clust size
  # Get motifs
  mots <- lapply(seq(1:length(unique(cut))),function(i){motifObj[,clustObj$labels[which(cut==i)]]})
  mot.score <- function(m,method=c("cosine")[1]) {
    require(lsa)
    # m is motif occurence matrix
    # rows = motif
    # columns = genes
    # value = # motifs X/ sum(# motifs)
    if (method == "cosine") {
      mot.score <- cosine(m)
    }
    m <- mean(mot.score[upper.tri(mot.score)])*dim(mot.score)[1]
    return(m)
  }
  mot.score <- mean()
}

# Bicluster quality scores
# From Zhao et al 2011
sd.within <- function(matrixIn) {
  # According to Zhao et al 2011
  o<-sqrt(sum((matrixIn-matrix(colMeans(matrixIn),nrow=dim(matrixIn)[1],ncol=dim(matrixIn)[2],
                            byrow=T))^2)/(dim(matrixIn)[1]*dim(matrixIn)[2]-1))
  return(o)
}

getActiveCorems <- function(condition,regulonConditionList,egrin=F) {
  if (egrin) {
    o <- names(regulonConditionList)[unlist(lapply(seq(1:length(regulonConditionList)),
                                                 function(i){condition%in%rownames(regulonConditionList[[i]])}))]
  } else {
    o <- names(regulonConditionList)[unlist(lapply(seq(1:length(regulonConditionList)),
                                                 function(i){condition%in%names(regulonConditionList[[i]])}))]
  }
  return(o)
}

makeConditionOntology <- function(oboFile) {
  require(ontoCAT)
  o <- getOntology(oboFile)
  return(o)
}

conditionEnrichment <- function(conditions,annotations,ontology,withParents=F,
                                pval.correct=T,method=c("BH","bonferroni")[1],return.all=F,
                                c.tot = NULL) { 
  #library(rJava)
  #options(java.parameters="-Xmx512")
  #.jinit()
  require(ontoCAT)
  if (withParents) {
    c.set<-lapply(annotations[conditions],function(i){
      i <- unlist(lapply(i,function(j){
        org <- j
        j <- getAllTermParentsById(ontology,gsub(":","_",j))
        j <- unlist(lapply(j,function(m){getAccession(m)}))
        j <- c(j,gsub(":","_",org))
        }))
    })
    c.set <- table(unlist(c.set))
    if (is.null(c.tot)) {
      c.tot<-lapply(annotations,function(i){
       i <- unlist(lapply(i,function(j){
         org <- j
          j <- getAllTermParentsById(ontology,gsub(":","_",j))
          j <- unlist(lapply(j,function(m){getAccession(m)}))
          j <- c(j,gsub(":","_",org))
         }))
      })
      c.tot <- table(unlist(c.tot))
    } else {
      c.tot = c.tot
    }
  } else {
    # Determine how many times each term occurs
    c.set <- table(unlist(annotations[conditions]))
    names(c.set) <- gsub(":","_",names(c.set))
    c.tot <- table(unlist(annotations))
    names(c.tot) <- gsub(":","_",names(c.tot))
  }
  o <- unlist(lapply(names(c.set),function(i){
    i <- phyper(c.set[i],c.tot[i],sum(c.tot)-c.tot[i],sum(c.set),lower.tail=F)
  }))
  # translate names
  n <- unlist(lapply(names(c.set),function(i){getTermNameById(ontology,i)}))
  names(o) <- n
  if (pval.correct) {
    p.adjust(o,method)
  }
  o <- sort(o)
  if (!return.all) {
    o <- o[o<0.05]
  }
  return(o)
}

getActiveCorems <- function(conditions,referenceConditionList=gBg_backbone_0.59_clean_list$conditions.cvar) {
  lapply(conditions,function(i){
    i <- names(referenceConditionList)[which(unlist(lapply(seq(1:length(referenceConditionList)),function(j){
      j<-i%in%names(referenceConditionList[[j]])
    })))]
  })
}

plotMotif <- function(motc,save=F,file="tmp.png") {
  require(Cairo)
  m=get.motif.cluster.info(motif.clusters=motc)[1]
  pssm=attr(m[[1]]$tttt.out[[1]],"combined.pssm")
  if (save) {
    png(file=file,bg="transparent")
    e$viewPssm(pssm,main=names(m))
    dev.off()
  }
  e$viewPssm(pssm,main=names(m))
}

findJumpingGenes <- function(r=o$corem_list,ratios=ratios.norm) {
  require(multicore)
  corems = r$corems
  # Compute overlap matrix
  m = t(combn(corems,2))
  # compute num genes in overlap
  overlap = unlist(lapply(seq(1:dim(m)[1]),function(i){i <- length(intersect(r$genes[[m[i,1]]],r$genes[[m[i,2]]]))}))
  m <- cbind(m,overlap)
  # recompute index
  m = m[which(m[,3]>0),]
  # Calculate score, mean cvar
  cvar <- function(genes,conditions) {
    genes <- intersect(rownames(ratios),genes)
    m <- abs(colMeans(ratios[genes,conditions,drop=F]))
    m[which(m==0)] = 1e-6
    var.m <- mean(apply(ratios[genes,conditions,drop=F],2,sd)/m)
    return(var.m)
  }
  jScore <- function(genes1,conditions1,genes2,conditions2) {
    if ((length(setdiff(conditions1,conditions2))==0) || (length(setdiff(conditions2,conditions1))==0)) {
      return(0)
    } else {
      e1.1 <- cvar(genes1,setdiff(conditions1,conditions2))
      e1.2 <- cvar(genes1,setdiff(conditions2,conditions1))
      e2.1 <- cvar(genes2,setdiff(conditions1,conditions2))
      e2.2 <- cvar(genes2,setdiff(conditions2,conditions1))
      if (length(intersect(conditions1,conditions2))==0) {
        e1.12=1; e2.12=1
      } else {
        e1.12 <- cvar(genes1,intersect(conditions1,conditions2))
        e2.12 <- cvar(genes2,intersect(conditions1,conditions2))
      }
      score <- (e2.1/e1.1) * (e1.2/e2.2) * (e1.12/e2.12)
      return(score)
    }
  }
  score = unlist(lapply(seq(1:dim(m)[1]),function(i){i <- jScore(r$genes[[m[i,1]]],intersect(names(r$conditions[[m[i,1]]]),colnames(ratios)),r$genes[[m[i,2]]],intersect(names(r$conditions[[m[i,2]]]),colnames(ratios)))}))
  l.r1 <- unlist(lapply(seq(1:dim(m)[1]),function(i){length(r$genes[[m[i,1]]])}))
  l.r2 <- unlist(lapply(seq(1:dim(m)[1]),function(i){length(r$genes[[m[i,2]]])}))
  m <- cbind(m,score,l.r1,l.r2)
  m <- m[order(as.numeric(m[,4]),decreasing=T),]
  colnames(m) <- c("Corem1","Corem2","Overlap","Score","#Genes1","#Genes2")
  return(m)
}

plotJumpingGenes <- function(r1,r2,r=gBg_backbone_0.59_clean_list,ratios=ratios.norm,...) {
  require(gplots)
  require(colorRamps)
  g1 <- r$genes[[r1]] 
  g2 <- r$genes[[r2]]
  c1 <- names(r$conditions.cvar[[r1]])
  c2 <- names(r$conditions.cvar[[r2]])
  rowmemb = unique(g1,g2); names(rowmemb) = rowmemb
  rowmemb[g1] = "red"
  rowmemb[g2] = "black"
  rowmemb[intersect(g1,g2)] = "yellow"
  #par(mfrow=c(1,3))
  rowmemb = c(rowmemb[which(rowmemb=="red")],rowmemb[which(rowmemb=="yellow")],rowmemb[which(rowmemb=="black")])
  heatmap.2(ratios[names(rowmemb),setdiff(c1,c2),drop=F],trace="none",dendrogram="none",col=blue2yellow(200),breaks=seq(-2,2,.02),RowSideColors=rowmemb,cexRow=.8,Rowv=F,main=paste("Red Conditions\n",r1,"\nYellow are shared",sep=""),...)
  readline()
  heatmap.2(ratios[names(rowmemb),intersect(c1,c2),drop=F],trace="none",dendrogram="none",col=blue2yellow(200),breaks=seq(-2,2,.02),RowSideColors=rowmemb,cexRow=.8,Rowv=F,main=paste("Red & Black Conditions\n",r1," & ",r2,"\nYellow are shared",sep=""))
  readline()
  heatmap.2(ratios[names(rowmemb),setdiff(c2,c1),drop=F],trace="none",dendrogram="none",col=blue2yellow(200),breaks=seq(-2,2,.02),RowSideColors=rowmemb,cexRow=.8,Rowv=F,main=paste("Black Conditions\n",r2,"\nYellow are shared",sep=""))
}

compareRegulonConditionality <- function(r=gBg_backbone_0.59_clean_list,ratios=ratios.norm) {
  m <- matrix(0,nrow=length(r$corems),ncol=length(colnames(ratios)),dimnames=list(r$corems,colnames(ratios)))
  for (regulon in r$corems) {
    g <- r$genes[[regulon]]
    c <- names(r$conditions.cvar[[regulon]])
    exp <- colMeans(ratios[g,c,drop=F])
    exp[exp>0] = 1
    exp[exp<0] = -1
    m[regulon,names(exp)] = exp  
  }
  return(m)
}

cosineRegulonConditionality <- function(m) {
  # m is matrix from compareRegulonCondiitonality. Computes the cosine similarity between all conditions in m.
  require(lsa)
  out <- matrix(0,ncol=dim(m)[2],nrow=dim(m)[2],dimnames=list(colnames(m),colnames(m)))
  for (i in 1:dim(m)[2]) {
    for (j in 1:dim(m)[2]) {
      out[i,j] = cosine(m[,i],m[,j])
    }
  }
  return(out)
}

reclusterMotifs <- function(motc,plot=T) {
  motifs <- unlist(get.motifs(motif.clusters=motc))
  ks=as.integer(sapply(strsplit(motifs,"_"),"[",2)) ## get the bicluster indexes
  mot.inds=as.integer(sapply(strsplit(motifs,"_"),"[",3)) ## get the motif indexes for each bicluster
  tt.out <- e$motif.similarities.tomtom(ks,ks,mot.inds,mot.inds,desymm=F,e.value.cutoff=Inf)
  tt.outa=subset(tt.out,p.value<=0.001) ## this cutoff seems to work well
  tt.out2 <- e$cluster.tomtom.results( tt.outa, n.cutoff=3,min.size=3,k.cut=0.99) ## uses hclust by default
  out=new.env()
  out$tt.out=tt.out
  out$tt.out2=tt.out2
  out$e=e
  sys.source("cmonkey-ensemble-funcs.R",envir=out)
  par( mfrow=c( 4, 4 ) )
  for ( i in 1:length( out$tt.out2 ) ) {
    print( i )
    tmp <- out$get.motif.cluster.info( paste( "MOTC", i, sep="_" ) )[[ 1 ]]
    out$e$viewPssm( attr( tmp, 'combined.pssm' ), main=paste( i, length( attr( tmp, 'mot.names' ) ) ) )
  }
  return(out)
}

getIntergenicSequences <- function(file="intergenic_bg.fa") {
  # [[1]] is pNRC100
  # [[2]] is chr
  # [[3]] is pNRC200
  # Remove entries without annotated start/stop
  genome.info <- e$genome.info$feature.tab[!is.na(e$genome.info$feature.tab[,"Start.new"]),]
  pNRC100 <- genome.info[which(genome.info[,"where"]=="pNRC100"),]
  pNRC100 <- pNRC100[order(pNRC100[,"Start.new"]),]
  chr <- genome.info[which(genome.info[,"where"]=="Chr"),]
  chr <- chr[order(chr[,"Start.new"]),]
  pNRC200 <- genome.info[which(genome.info[,"where"]=="pNRC200"),]
  pNRC200 <- pNRC100[order(pNRC200[,"Start.new"]),]
  l = list()
  findIntergenic <- function(table,l,index,file) {
    for (i in seq(1:dim(table)[1])) {
      if (i == dim(table)[1]) {
        # Last gene
        # Check orientations
        orient <- paste(c(as.character(table[1,"strand"]),as.character(table[i,"strand"])),collapse="")
        if (orient=="DD") {
          # point same direction, forward
          end.new <- table[1,"Start.new"]
          start.new <- table[i,"Stop.new"]
          nt <- paste(substr(e$genome.info$genome.seqs[[index]],start.new,nchar(e$genome.info$genome.seqs)),
                      substr(e$genome.info$genome.seqs[[index]],1,end.new),sep="")
          if (nchar(nt)>50) {
            l = c(l,nt)
            write(paste(">",nt,sep="\n"),file=file,append=T)
          }
        } else if (orient=="RR") {
          end.new <- table[1,"Stop.new"]
          start.new <- table[i,"Start.new"]
          nt <- paste(substr(e$genome.info$genome.seqs[[index]],start.new,nchar(e$genome.info$genome.seqs)),
                      substr(e$genome.info$genome.seqs[[index]],1,end.new),sep="")
          if (nchar(nt)>50) {
            l = c(l,nt)
            write(paste(">",nt,sep="\n"),file=file,append=T)
          }
        } else if (orient=="DR") {
          end.new <- table[1,"Start.new"]
          start.new <- table[i,"Start.new"]
          nt <- paste(substr(e$genome.info$genome.seqs[[index]],start.new,nchar(e$genome.info$genome.seqs)),
                      substr(e$genome.info$genome.seqs[[index]],1,end.new),sep="")
          if (nchar(nt)>50) {
            l = c(l,nt)
            write(paste(">",nt,sep="\n"),file=file,append=T)
          }
        }
      } else {
        orient <- paste(c(as.character(table[i,"strand"]),as.character(table[i+1,"strand"])),collapse="")
        if (orient=="DD") {
          # point same direction, forward
          end.new <- table[i+1,"Start.new"]
          start.new <- table[i,"Stop.new"]
          nt <- paste(substr(e$genome.info$genome.seqs[[index]],start.new,end.new),sep="")
          if (nchar(nt)>50) {
            l = c(l,nt)
            write(paste(">",nt,sep="\n"),file=file,append=T)
          }
        } else if (orient=="RR") {
          end.new <- table[i+1,"Stop.new"]
          start.new <- table[i,"Start.new"]
          nt <- paste(substr(e$genome.info$genome.seqs[[index]],start.new,end.new),sep="")
          if (nchar(nt)>50) {
            l = c(l,nt)
            write(paste(">",nt,sep="\n"),file=file,append=T)
          }
        } else if (orient=="RD") {
          end.new <- table[i+1,"Start.new"]
          start.new <- table[i,"Start.new"]
          nt <- paste(substr(e$genome.info$genome.seqs[[index]],start.new,end.new),sep="")
          if (nchar(nt)>50) {
            l = c(l,nt)
            write(paste(">",nt,sep="\n"),file=file,append=T)
          }
        }
      }
    }
    return(l)
  }
  l = findIntergenic(pNRC100,l,1,file)
  l = findIntergenic(chr,l,2,file)
  l = findIntergenic(pNRC200,l,3,file)
  invisible(l)
}

pCoReg <- function(r = gBg_backbone_0.59_clean_list,ref=gBg_backbone_0.59_clean) {
  require(multicore)
  require(bigmemory)
  g <- sort(unique(unlist(r$genes)))
  m <- matrix(0,nrow=length(g),ncol=length(g),dimnames=list(g,g))
  c <- unique(unlist(sapply(r$conditions.cvar,names)))
  g.c<-t(combn(g,2))
  s <- mclapply(seq(1:dim(g.c)[1]),
              function(i){
                r1 <- getCorems(g.c[i,1],ref)
                r2 <- getCorems(g.c[i,2],ref)
                r.c <- intersect(r1,r2)
                if (length(r.c)>0) {
                  conds <- unique(unlist(sapply(r$conditions.cvar[r.c],names)))
                  i <- length(conds)/length(c)
                } else {
                  i <- 0
                }
                return(i)
                })
  m[g.c] <- unlist(s)
  m[cbind(g.c[,2],g.c[,1])] <- unlist(s)
  return(m)
}

plotGene.p.reg <- function(gene,p=p.coreg,ref=gBg_backbone_0.59_clean) {
  require(gplots)
  # Find nonzero xpression values
  g.p <- sort(p[gene,][p[gene,]>0])
  r <- getCorems(gene,ref)
  r.g <- lapply(r,function(i){getGenes(i,ref)})
  p.regs <- lapply(names(g.p),function(i){
    i<-which(unlist(lapply(r.g,function(j){i%in%j})))})
  id <- sapply(p.regs,sum)
  names(id) <- names(g.p)
  colTable <- sort(table(id))
  colPal <- c(brewer.pal(8,"Accent"),brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"),
              brewer.pal(8,"Pastel2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))
  colRef <- sample(colPal,length(colTable)); names(colRef) <- names(colTable)
  id.col <- unlist(lapply(id,function(i){colRef[which(names(colRef)==i)]}))
  plot(g.p,col=id.col,pch=19,xlab="sorted gene index",ylab="conditional probability",main=gene)
}

calculateEigenGene <- function(corem,ratios=ratios.norm,ref=o$corem_list,alt.c=F,remove=c("genes","conditions")[1]) {
  #print(corem)
  g <- intersect(ref$genes[[corem]],rownames(ratios))
  if (!is.logical(alt.c)) {
    conds <- intersect(alt.c,colnames(ratios))
  } else {
    conds <- names(ref$conditions.cvar[[regulon]])
  }
  if (remove=="genes") {
    # remove genes with na values 
    ratios <- ratios[g,conds,drop=F]
    to.remove <- which(apply(ratios[g,conds],1,function(i)sum(is.na(i)))>0)
    if (length(to.remove)>0) {
      ratios <- ratios[-to.remove,,drop=F]
      g <- g[-to.remove]
    } 
  } else if (remove=="conditions") {
    # remove conditions with na values 
    ratios <- ratios[g,conds,drop=F]
    to.remove <- which(apply(ratios[g,conds],1,function(i)sum(is.na(i)))>0)
    if (length(to.remove)>0) {
      ratios <- ratios[,-to.remove,drop=F]
      conds <- conds[-to.remove]
    }
  }
  if ((dim(ratios)[1]<=1||dim(ratios)[2]<=1)) {
    return(NULL)
  } else {
    exp.svd <- svd(ratios)
    # Take only first eigengene
    o <- list()
    o$profile <- -exp.svd$v[,1]; names(o$profile) <- conds
    o$profile.component <- o$profile*exp.svd$d[1]; names(o$profile.component) <- conds
    o$var.explained <- exp.svd$d[1]^2/sum(exp.svd$d^2)
    # Best match by correlation
    exp.cor <- unlist(lapply(g,function(i){cor(ratios[i,conds],o$profile)})); names(exp.cor) <- g
    o$representative <- exp.cor[which(exp.cor==max(exp.cor))]
    names(o$representative) <- names(exp.cor)[which(exp.cor==max(exp.cor))]
    return(o)
  }
}

conditionCoherenceTest <- function(coremStruct=o$corem_list,cutoff=0.05,
                                   ratios=ratios.norm,method=c("sd","cvar")[2],type=c("egrin2","cvar")[1]) {
  # Find
  require(multicore)
  cvar <- function(genes,conditions) {
    m <- abs(colMeans(ratios[genes,conditions,drop=F]))
    m[which(m==0)] = 1e-6
    var.m <- median(apply(ratios[genes,conditions,drop=F],2,sd)/m)
    return(var.m)
  }
  r <- coremStruct$corems
  o <- unlist(mclapply(r,function(i){
    g <- coremStruct$genes[[i]]
    if (type == "egrin2") {
      c <- coremStruct$conditions.egrin2[[i]]
      c <- rownames(c)[c[,2]<cutoff]
    } else if (type == "cvar") {
      c <- coremStruct$conditions.cvar[[i]]
      c <- names(c)[c<cutoff]
    }
    if (length(c)==0) {
      return(NULL)
    } else {
      if (method=="sd") {
       i <- median(apply(ratios[g,c,drop=F],2,sd))  
     } else if (method=="cvar") {
        i <- cvar(g,c)
      }
      return(i)
    }
  }))
  return(o)
}

writeMapEquation <- function(condition1,condition2=NULL,file="tmp.map",regCond=m,
                             ref=gBg_backbone_0.59_clean_list,ratios=ratios.norm) {
  # load("/isb-1/R/egrin2/gBg_antoine/gBg_backbone_corems.RData")
  # load("/isb-1/R/egrin2/gBg_antoine/regulonConditionality/regCond.RData")
  # load("/isb-1/R/egrin2/egrin2_ratios.RData")
  # condition1 = "light2__1440m_vs_light1_3960m"
  # condition2 = "dark1__0240m_vs_NRC-1c"
  # writeMapEquation()
  # get active corems for this condition
  r <- names(regCond[,condition1])[regCond[,condition1]!=0]
  if (!is.null(condition2)) {
    r2 = names(regCond[,condition2])[regCond[,condition2]!=0]
    r = setdiff(r,r2)
  }
  g <- ref$genes[r]
  eigen <- ref$eigengene[r]
  g.unique <- unique(unlist(g))
  g.unique.r <- lapply(g.unique,function(i){
    r.tmp <- intersect(r,getCorems(i,gBg_backbone_0.59_clean))
    s <- sapply(r.tmp,function(j){
      cor(ratios[i,names(eigen[[j]]$profile)],eigen[[j]]$profile)
    })
    s <- s[which(s==max(s))][1]
    return(names(s))
    })
  names(g.unique.r) <- g.unique
  r <- unique(unlist(g.unique.r))
  g <- lapply(r,function(i){
    i <- names(unlist(g.unique.r))[grep(i,unlist(g.unique.r))]
  })
  names(g) <- r
  r.links <- ref$regulon.matrix.gene[r,r]
  # Write header
  write(paste("# modules: ",length(r),sep=""),file)
  write(paste("# modulelinks: ",sum(r.links>0)/2,sep=""),file,append=T)
  write(paste("# nodes: ",length(unique(unlist(g))),sep=""),file,append=T)
  write(paste("# links: ",length(unique(unlist(g))),sep=""),file,append=T)
  write(paste("# codelength: ","2.51912",sep=""),file,append=T)
  write("*Undirected",file,append=T)
  write(paste("*Modules ",length(r),sep=""),file,append=T)
  for (i in seq(1:length(r))) {
    module = r[i]
    # Var explained by eigengene
    var = ref$eigengene[[module]]$var.explained
    # "Flow" -- mean overlap with other corems
    flow = sum(r.links[module,])/length(r.links[module,]>0)
    write(paste(i,paste('"',module,'"',sep=""),var,flow,sep=" "),file,append=T)
  }
  write(paste("*Nodes ",length(unique(unlist(g))),sep=""),file,append=T)
  for (i in seq(1:length(r))) {
    module = r[i]
    eigengene = ref$eigengene[[module]]$profile
    g.exp = ratios.norm[g[[module]],names(eigengene),drop=F]
    g.cor = sort(apply(g.exp,1,function(i){cor(i,eigengene)}),decreasing=T)
    for (j in seq(1:length(g.cor))) {
      write(paste(paste(i,j,sep=":"),paste('"',names(g.cor)[j],'"',sep=""),g.cor[j]),sep=" ",file,append=T)
    }
  }
  write(paste("*Links ",sum(r.links>0)/2,sep=" "),file,append=T)
  for (i in seq(1:dim(r.links)[1])) {
    cm <- which(r.links[i,]>0)
    if (length(cm)>0) {
      for (j in cm) {
        write(paste(i,j,r.links[i,j]),file,append=T)
      }
    }
  }
}

corem_coherence <- function(ref=gBg_backbone_0.59_clean_list,ratios=ratios.norm,groupID=groupID,
                            resamples=15000) {
  # Are corems coherent across experimental sets? or are they
  # an artifact of coincidental coexpression
  # Use resampling to find out
  require(multicore)
  cvar <- function(genes,conditions) {
    m <- abs(colMeans(ratios[genes,conditions,drop=F]))
    m[which(m==0)] = 1e-6
    var.m <- mean(apply(ratios[genes,conditions,drop=F],2,sd)/m)
    return(var.m)
  }
  # Find full complement of conditions. All experiments in a particular group
  full.conds <- mclapply(seq(1:length(ref$conditions.cvar)),function(i){
    ids <- unique(groupID[names(ref$conditions.cvar[[i]])])
    i <- unlist(lapply(ids,function(j){names(groupID)[which(groupID==j)]}))
    return(i)
  })
  names(full.conds) <- ref$corems
  cvar.test <- mclapply(ref$corems,function(i){
    i <- cvar(ref$genes[[i]],full.conds[[i]]) 
  })
  names(cvar.test) <- ref$corems
  # resample random gene sets
  cvar.resample <- mclapply(ref$corems,function(i) {
    out <- lapply(seq(1:resamples),function(j) {
      j <- cvar(sample(rownames(ratios),length(ref$genes[[i]])),full.conds[[i]])
    })
  })
  names(cvar.resample) <- ref$corems
  out <- mclapply(names(cvar.test),function(i){
    i <- sum(cvar.test[[i]]>=cvar.resample[[i]])/length(cvar.resample[[i]])
  })
  names(out) <- ref$corems
  return(out)
}

assessBrokenOperons <- function() {
  tmp.operon <- lapply(as.character(operons[,3]),function(a){
    o <- strsplit(a,"[\\*\\| ]")[[1]]
    o.n <- sapply(o,nchar)
    o <- o[o.n>0]
    return(o)
    })
  tmp.operon.broken <- lapply(seq(1:length(tmp.operon)),function(i){
    r <- table(unlist(lapply(tmp.operon[[i]],function(j){
      o<-getCorems(j,gBg_backbone_0.59_clean)
      return(o)
    })))
    if (length(which(r<length(tmp.operon[[i]])))>0) {
      if (length(r)>0) {
        # broken
        return(1)
      } else {
        return(NA)
      }
    } else {
      if (length(r)>0) {
        # not broken
        return(0)
      } else {
        return(NA)
      }
    }
  })
}

motif2 <- function(regulon=NULL,genes=NULL,quant=.1,cutoff=NULL,window=500,shift=250,motif.cutoff = 1e-6,corem.table=o$corems,
                   freq.cutoff=.05) {
  # Tries to find the best motif.cluster explaining coregulation of
  # genes in a corem.
  # Uses intersection of user defined upper quantile of biclusters and mast hits 
  # to find potentially relevant motif clusters
  # Requires egrin2 to be loaded
  # make sure envs are right
  environment(get.mast) <- out
  if (!is.null(regulon)) {
    g<-getGenes(regulon,corem.table)
  } else if (!is.null(genes)) {
    g<-genes
  } else {
    cat("you need to provide either a vector of genes or a regulon number and regulon table\n")
    return(invisible(NULL))
  }
  bcs <- sort(table(unlist(get.biclusters(g))),decreasing=T)
  q <- rev(quantile(bcs,probs=seq(0,1,quant)))[2]
  # filter bcs below user supplied quantile
  bcs <- names(bcs)[bcs>=q]
  out<-list()
  out$bc.motifs <- agglom(src=bcs,srcType="bicluster",targetType="motif.cluster",path="motif")
  if (!is.null(cutoff)) {
    out$bc.motifs<- out$bc.motifs[out$bc.motifs[,2]<=cutoff,]
  } else {
    out$bc.motifs <- out$bc.motifs[out$bc.motifs[,3]<=0.05,]
  }
  out$mast.motifs <- lapply(g,function(i){
    #print(i)
    n = i
    genome.loc <- e$get.gene.coords(n,op.shift=F)
    if (is.null(genome.loc)) {
      return(0)
    }
    if (genome.loc$strand == "D") {
      ss <- genome.loc$start_pos
      i <- get.mast(ss,window=window,shift=shift,p.value.cutoff=motif.cutoff)
    } else if (genome.loc$strand == "R") {
      ss <- genome.loc$end_pos
      i <- get.mast(ss,window=window,shift=shift,p.value.cutoff=motif.cutoff)
    }
    i <- as.data.frame(i)
    #m <- agglom(i,"motif.cluster","bicluster,motif","motif")
    if (dim(i)[1]>0) {
      motclust <- sapply(seq(1,dim(i)[1]),function(row) {
        return(try(unlist(get.motif.clusters(paste("MOT",i[row,1],abs(i[row,2]),sep="_")))))
      })
      empty.ind <- sapply(motclust,is.null)
      i <- i[!empty.ind,]
      i$motif.clusters <- unlist(motclust)
      genome.loc <- e$get.gene.coords(n,op.shift=F)
      if (genome.loc$strand == "D") {
        ss <- genome.loc$start_pos
        i$distance <- i$posns-ss
        j<-i[intersect(which(i$distance>=(-window-shift)),which(i$distance<=(window-shift))),]
      } else if (genome.loc$strand == "R") {
        ss <- genome.loc$end_pos
        i$distance <- ss-i$posns
        j<-i[intersect(which(i$distance>=(-window-shift)),which(i$distance<=(window-shift))),]
      }
      if (dim(j)[1]>0) {
        return(j)
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  })
  names(out$mast.motifs) <- g
  out$common.motifs.genes <- lapply(g,function(i){
    #print(i)
    if (out$mast.motifs[[i]]==0) {
      return(0)
    } 
    common.m <- intersect(unique(out$mast.motifs[[i]]$motif.clusters),rownames(out$bc.motifs))
    # score is fract m1 * fract m2
    if (length(common.m)>0) {
      m.score <- sapply(common.m,function(j){
        #mast.sum <- sum(out$mast.motifs[[i]])
        bc.sum <- sum(out$bc.motifs[,1])
        #return(((out$mast.motifs[[i]][j]/mast.sum)*(out$bc.motifs[j,1]/bc.sum)))
        s <- out$bc.motifs[j,1]/bc.sum
        if (s>=freq.cutoff) {
          return(s)
        } else {
          return(0)
        }
      })
      names(m.score) <- common.m
    } else {
      m.score = 0
    }
    return(m.score)
  })
  names(out$common.motifs.genes) <- g
  # clean up -- remove MOTC with 0 weight
  out$common.motifs.genes <- lapply(out$common.motifs.genes,function(i)i<-i[which(i>0)])
  # try to save genes with no MOTC annotation. if in operon, assign leader gene
  empty.ind <- names(sapply(out$common.motifs.genes,length))[which(sapply(out$common.motifs.genes,length)==0)]
  # which of these is in an operon
  translation <- sapply(empty.ind,function(i){
    ind <- grep(i,e$operon.list())
    if (length(ind)>0) {
      # get header gene
      return(names(e$operon.list())[ind])
    }
  })
  # remove entries with no matches
  translation <- translation[which(sapply(translation,length)>0)]
  # keep only entries that are also in corem
  translation <- translation[translation%in%names(out$common.motifs.genes)]
  # replace translated entries
  out$common.motifs.genes[names(translation)] <-  out$common.motifs.genes[unlist(translation)]
  out$common.motifs.all <- sort(table(unlist(sapply(out$common.motifs.genes,names))),decreasing=T)/sum(table(unlist(sapply(out$common.motifs.genes,names))))
  # clean up mast hits by collapsing hits that occur within 25 bp
  out$mast.motifs <- lapply(seq(1,length(out$mast.motifs)),function(i){
    if(!is.null(dim(out$mast.motifs[[i]])[1])) {
      unique.m <- unique(out$mast.motifs[[i]][,"motif.clusters"])
      # sort
      unique.m <- unique.m[order(as.numeric(sapply(strsplit(unique.m,"_"),"[[",2)))]
      unique.m.i <- lapply(unique.m,function(j){which(sapply(out$mast.motifs[[i]][,"motif.clusters"],function(x)(j==x)))})
      names(unique.m.i) <- unique.m
      toReturn <- data.frame()
      for (m in unique.m) {
        sub.t <- out$mast.motifs[[i]][unique.m.i[[m]],]
        while(min(as.matrix(dist(sub.t[,"distance"]))[lower.tri(as.matrix(dist(sub.t[,"distance"])))])<25) {
          # collpase row with greatest similarities
          row.i <- which(apply(as.matrix(dist(sub.t[,"distance"])),2,function(x)sum(x<25))==max(apply(as.matrix(dist(sub.t[,"distance"])),2,function(x)sum(x<25))))[1]
          row.i.toc <- which(as.matrix(dist(sub.t[,"distance"]))[,row.i]<25)
          toCollapse <- sub.t[row.i.toc,]
          sub.t <- sub.t[-row.i.toc,]
          toInsert <- toCollapse[1,,drop=F]; toInsert[1,"mots"] = median(sign(toCollapse[,"mots"])); toInsert[1,"pvals"] = mean(toCollapse[,"pvals"])
          toInsert[1,"posns"] <- median(toCollapse[,"posns"]); toInsert[1,"distance"] <- median(toCollapse[,"distance"])
          sub.t <- rbind(sub.t,toInsert)
        }
        toReturn <- rbind(toReturn,sub.t)
      }
      if (dim(toReturn)[1]>0) {
        return(toReturn)
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  })
  names(out$mast.motifs) <- g
  out$mast.motifs.common <- lapply(g,function(i){
    if(!is.null(dim(out$mast.motifs[[i]])[1])) {
      mots2keep <- names(out$common.motifs.genes[[i]])
      unique.m.i <- sort(unlist(lapply(mots2keep,function(j){which(sapply(out$mast.motifs[[i]][,"motif.clusters"],function(x)(j==x)))})))
      if (dim(out$mast.motifs[[i]][unique.m.i,])[1]>0) {
        return(out$mast.motifs[[i]][unique.m.i,])
      } else {
        return(0) 
      }
    } else {
      return(0)
    }
  })
  names(out$mast.motifs.common) <- g
  # only keep motifs that are common -- i.e. also in biclusters
  # make final motif table
  out$table <- data.frame()
  for (i in seq(1,length(out$mast.motifs.common))) {
    if (!is.null(dim(out$mast.motifs.common[[i]])[1])) {
      if (names(out$mast.motifs.common)[i]%in%names(translation)) {
        ref.gene = translation[names(out$mast.motifs.common)[i]]
      } else {
        ref.gene = names(out$mast.motifs.common)[i]
      }
      toAdd <- data.frame("Gene"=names(out$mast.motifs.common)[i],"RefGene"=ref.gene,"GRE"=out$mast.motifs.common[[i]][,"motif.clusters"],"DistanceToTSS"=out$mast.motifs.common[[i]][,"distance"],"Orientation"=sign(out$mast.motifs.common[[i]][,"mots"]))
      out$table <- rbind(out$table,toAdd)
    }
  }
  return(out)
}


#######################
#
# Motif regulatory network
#
#
########################
# This needs to be attached into EGRIN2 env. Usually called "out"
# function motif2 also needs to be in "out"

make.motif.reg.network <- function(corems=NULL,genes=NULL,gBg_backbone = NULL,corem.table = NULL,outdir = NULL, 
                                   genesTohighlight=NULL,fitness=F,multicore=T,cores=8,altCoremPie=F,
                                   quant=.1,cutoff=NULL,window=1000,shift=375,motif.cutoff = 1e-6,
                                   freq.cutoff=.05,pie.pdf=F) {
  require(data.table)
  require(igraph)
  require(multicore)
  options(cores=cores)
  require(RColorBrewer)
  if (fitness) {
    load("/isb-1/R/ecoli/chemgen/chemgen.RData")
  }
  if (!is.null(corems)) {
    gene.list <- mclapply(corems,function(i){getGenes(i,corem.table)})
    names(gene.list) <- corems
    all.genes <- unique(unlist(gene.list))
    # get edges to add
    edge.ind <- c()
    for (i in corems) {
      edge.ind<-rbind(edge.ind,cbind(as.matrix(corem.table[i,])[,2],"(pp)",as.matrix(corem.table[i,])[,3],"=",i))
    }
    edge.ind <- rbind(edge.ind,cbind(edge.ind[,3],edge.ind[,2],edge.ind[,1],edge.ind[,4:5]))
    
    #   motif.comp <- list()
    #   for (i in corems) {
    #     print(i)
    #     motif.comp[[i]] <- motif2(i,corem.table=corem.table)
    #   }
    if (multicore) {
      motif.comp <- mclapply(seq(1,length(corems)),function(i){
        i <- motif2(corems[i],corem.table=corem.table,quant=quant,cutoff=cutoff,window=window,shift=shift,motif.cutoff = motif.cutoff,
                    freq.cutoff=freq.cutoff)
        return(i)
      })
    } else {
      motif.comp <- lapply(seq(1,length(corems)),function(i){
        i <- motif2(corems[i],corem.table=corem.table,quant=quant,cutoff=cutoff,window=window,shift=shift,motif.cutoff = motif.cutoff,
                    freq.cutoff=freq.cutoff)
        return(i)
      })
    }
    names(motif.comp) <- corems
    # make graph
    sub.m <- gBg_backbone[all.genes,all.genes]; sub.m[] = 0
    sub.m[cbind(edge.ind[,1],edge.ind[,3])] <- gBg_backbone[cbind(edge.ind[,1],edge.ind[,3])]
  } else if(!is.null(genes)) {
    all.genes <- genes
    motif.comp <- list(motif2(genes=genes,quant=quant,cutoff=cutoff,window=window,shift=shift,motif.cutoff = motif.cutoff,
                              freq.cutoff=freq.cutoff))
    names(motif.comp) = c("1")
    sub.m <- gBg_backbone[all.genes,all.genes]
  }
  dir.create(outdir)
  
  graph <- graph.adjacency(sub.m,mode="undirected",weighted=T)
  write.graph(graph,file=paste(outdir,"/graph.graphml",sep=""),format="graphml")
  
  # node mapping
  g2n <- paste("n",seq(0,(length(all.genes)-1)),sep=""); names(g2n) <- V(graph)$name
  
  # write edge attribute corems
  # convert table names
  if (!is.null(corems)) {
    table.toWrite <- cbind(g2n[edge.ind[,1]],edge.ind[,2],g2n[edge.ind[,3]],edge.ind[,4:5])
    
    write("Corem",file=paste(outdir,"/corem.txt",sep=""))
    write.table(table.toWrite,
                file=paste(outdir,"/corem.txt",sep=""),
                append=T,row.names=F,col.names=F,sep=" ",quote=F)
  }
  
  if (fitness) {
    inChemgen<- which(apply(cbind(edge.ind[,1]%in%rownames(chemgen.cc),edge.ind[,3]%in%rownames(chemgen.cc)),1,sum)==2)
    edge.ind.sub <- edge.ind[inChemgen,]
    table.toWrite <- cbind(g2n[edge.ind.sub[,1]],edge.ind.sub[,2],g2n[edge.ind.sub[,3]],edge.ind.sub[,4],chemgen.cc[cbind(edge.ind.sub[,1],edge.ind.sub[,3])])
    write("Fitness",file=paste(outdir,"/fitness.txt",sep=""))
    write.table(table.toWrite,
                file=paste(outdir,"/fitness.txt",sep=""),
                append=T,row.names=F,col.names=F,sep=" ",quote=F)
  }
  
  # if present make node attributes for genesTohighlight
  if (!is.null(genesTohighlight)) {
    g <- intersect(genesTohighlight,names(g2n))
    table.toWrite <- cbind(g2n[g],"=","TRUE")
    write("Genes to highlight",file=paste(outdir,"/genesTohighlight.txt",sep=""))
    write.table(table.toWrite,
                file=paste(outdir,"/genesTohighlight.txt",sep=""),
                append=T,row.names=F,col.names=F,sep=" ",quote=F)
  }
  
  # Construct motif pie charts
  motifs <- unique(unlist(mclapply(motif.comp,function(i)names(i$common.motifs.all))))
  # order motif clusters
  m.order <- order(as.numeric(sapply(motifs,function(i)strsplit(i,split="_")[[1]][2])))
  motifs <- motifs[m.order]
  names(motifs) <- motifs; motifs[] <- 1/length(motifs)
  
  b.colors <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Pastel1"))
  motif.colors <- motifs; names(motif.colors) <- names(motifs); motif.colors[] <- b.colors
  
  # Key
  pdf(paste(outdir,"/key.pdf",sep=""))
  pie(as.numeric(motifs),col=motif.colors,labels=names(motifs))
  dev.off()
  
  # For genes
  gene.motif.prevalence <- matrix(0,nrow=length(all.genes),ncol=length(motifs),dimnames=list(all.genes,names(motifs)))
  for (i in names(motif.comp)) {
    for (j in names(motif.comp[[i]]$common.motifs.genes)) {
      print(j)
      for (m in names(motif.comp[[i]]$common.motifs.genes[[j]])) {
        print(m)
        gene.motif.prevalence[j,m] <- motif.comp[[i]]$common.motifs.genes[[j]][m]
      }
    }
  }
  
  # normalize rows
  gene.motif.prevalence <- t(apply(gene.motif.prevalence,1,function(i){return(i/sum(i))}))
  # check to see if need to retranspose
  if (dim(gene.motif.prevalence)[1]==1) {
    gene.motif.prevalence <- t(gene.motif.prevalence)
    colnames(gene.motif.prevalence) <- names(motifs)
  }
  # remove NaNs -- make 0
  gene.motif.prevalence[is.nan(gene.motif.prevalence)] <- 0
  
  if (!is.null(corems)){
    # For corems
    if (altCoremPie == F) {
      gene.motif.prevalence.corem <- matrix(0,nrow=length(motif.comp),ncol=length(motifs),dimnames=list(names(motif.comp),names(motifs)))
      for (i in names(motif.comp)) {
        for (j in names(motif.comp[[i]]$common.motifs.all)) {
          # get all gene scores
          if (dim(gene.motif.prevalence)[2]==1) {
            score = 1
            names(score) = j
          } else {
            score <- apply(gene.motif.prevalence[gene.list[[i]],],2,sum)
          }
          gene.motif.prevalence.corem[i,j] <- score[j]
        }
      }
      # normalize rows
      gene.motif.prevalence.corem <- t(apply(gene.motif.prevalence.corem,1,function(i){return(i/sum(i))}))
      # remove NaNs -- make 0
      gene.motif.prevalence.corem[is.nan(gene.motif.prevalence.corem)] <- 0
    } else if (altCoremPie == T) {
      # This is a little hacked together. Would like to rework and make more robust. 
      gene.motif.prevalence.corem <- matrix(0,nrow=length(motif.comp),ncol=length(motifs),dimnames=list(names(motif.comp),names(motifs)))
      for (i in names(motif.comp)) {
        for (j in names(motif.comp[[i]]$common.motifs.all)) {
          # get all gene scores
          score <- motif.comp[[i]]$bc.motifs[j,1]/sum(motif.comp[[i]]$bc.motifs[names(motif.comp[[i]]$common.motifs.all),1])
          gene.motif.prevalence.corem[i,j] <- score
        }
      }
    }
  }
  
  
  dir.create(paste(outdir,"/pie/",sep=""))
  for (i in 1:dim(gene.motif.prevalence)[1]) {
    if (pie.pdf) {
      pdf(paste(outdir,"/pie/",rownames(gene.motif.prevalence)[i],".pdf",sep=""),
          width=2,height=2)
    } else {
      png(paste(outdir,"/pie/",rownames(gene.motif.prevalence)[i],".png",sep=""),
          type="quartz",bg="transparent",width=960,height=960,units="px",res=300)
    }
    
    if (sum(gene.motif.prevalence[i,])>0) {
      pie(gene.motif.prevalence[i,names(motifs)],labels=NA,col=motif.colors)
    } else {
      pie(1,labels=NA,col="white")
    }
    dev.off()
  }
  
  if (!is.null(corems)){
    for (i in corems) {
      if (pie.pdf) {
        pdf(paste(outdir,"/pie/",i,".pdf",sep=""),
            width=2,height=2)
      } else {
        png(paste(outdir,"/pie/",i,".png",sep=""),
            type="quartz",bg="transparent",width=960,height=960,units="px",res=300)
      } 
      if (sum(gene.motif.prevalence.corem[i,])>0) {
        pie(gene.motif.prevalence.corem[i,names(motifs)],labels=NA,col=motif.colors)
      } else {
        pie(1,labels=NA,col="white")
      }
      dev.off()
    }
  }
  # write node attribute motif pie chart
  write("Motif pie chart",file=paste(outdir,"/motifs.txt",sep=""))
  write.table(cbind(g2n[all.genes],"=",paste("file://",outdir,"/pie/",all.genes,".png",sep="")),
              file = paste(outdir,"/motifs.txt",sep=""),
              append=T,row.names=F,col.names=F,sep=" ",quote=F)
  
  # write table of motif posns
  write("",file=paste(outdir,"/motifs_posns.txt",sep=""))
  for (i in seq(1,length(motif.comp))) {
    write(names(motif.comp)[i],file=paste(outdir,"/motifs_posns.txt",sep=""),append=T)
    write.table(motif.comp[[i]]$table,file = paste(outdir,"/motifs_posns.txt",sep=""),
                append=T,row.names=F,col.names=T,sep="\t",quote=F)
  }
}

makeConditionOntology <- function(oboFile) {
  require(ontoCAT)
  o <- getOntology(oboFile)
  return(o)
}

conditionEnrichment <- function(conditions,annotations,ontology=NULL,withParents=F,
                                pval.correct=T,method=c("BH","bonferroni")[1],return.all=F,
                                c.tot = NULL) { 
  #library(rJava)
  #options(java.parameters="-Xmx512")
  #.jinit()
  if (!is.null(ontology)) {
    require(ontoCAT)
    if (withParents) {
      c.set<-lapply(annotations[conditions],function(i){
        i <- unlist(lapply(i,function(j){
          org <- j
          j <- getAllTermParentsById(ontology,gsub(":","_",j))
          j <- unlist(lapply(j,function(m){getAccession(m)}))
          j <- c(j,gsub(":","_",org))
        }))
      })
      c.set <- table(unlist(c.set))
      if (is.null(c.tot)) {
        c.tot<-lapply(annotations,function(i){
          i <- unlist(lapply(i,function(j){
            org <- j
            j <- getAllTermParentsById(ontology,gsub(":","_",j))
            j <- unlist(lapply(j,function(m){getAccession(m)}))
            j <- c(j,gsub(":","_",org))
          }))
        })
        c.tot <- table(unlist(c.tot))
      } else {
        c.tot = c.tot
      }
    } else {
      # Determine how many times each term occurs
      c.set <- table(unlist(annotations[conditions]))
      names(c.set) <- gsub(":","_",names(c.set))
      c.tot <- table(unlist(annotations))
      names(c.tot) <- gsub(":","_",names(c.tot))
    }
    o <- unlist(lapply(names(c.set),function(i){
      i <- phyper(c.set[i],c.tot[i],sum(c.tot)-c.tot[i],sum(c.set),lower.tail=F)
    }))
    # translate names
    n <- unlist(lapply(names(c.set),function(i){getTermNameById(ontology,i)}))
    names(o) <- n
  } else {
    # Determine how many times each term occurs
    c.set <- table(unlist(annotations[conditions]))
    c.tot <- table(unlist(annotations))
    o <- unlist(lapply(names(c.set),function(i){
      i <- phyper(c.set[i],c.tot[i],sum(c.tot)-c.tot[i],sum(c.set),lower.tail=F)
    }))
  }
  if (pval.correct) {
    p.adjust(o,method)
  }
  o <- sort(o)
  if (!return.all) {
    o <- o[o<0.05]
  }
  return(o)
}

getGO <- function(genes,gene2entrez,class=c("BP","MF","CC")[1],return.all=F,pval=5) {
  require(DAVIDQuery)
  g <- gene2entrez[intersect(genes,names(gene2entrez))]
  to.r <- try(DAVIDQuery(g,type="ENTREZ_GENE_ID",tool="chartReport",annot=paste("GOTERM_",class,"_FAT",sep=""),verbose=F,URLlengthLimit = 20004800,formatIt=T)$DAVIDQueryResult)
  if (class(to.r)=="try-error") {
    return(NULL)
  } else {
    if (!return.all) {
      if (dim(to.r)[1]>0) {
        if (to.r[1,"FDR"]=="FDR") {
          to.r<-to.r[2:dim(to.r)[1],]
        }
        to.r <- to.r[]
        to.r <- to.r[as.numeric(to.r[,"FDR"])<=pval,]
        if (dim(to.r)[1]==0) {
          to.r <- NULL
        }
      } else {
        to.r <- NULL
      }
    }
    return(to.r)
  }
}

