# To address reviewer concerns for Brooks, Reiss et al 2014 MSB regarding model overfitting.

# Basic idea: resample backbone and corems when some of the runs (with specific conditions) are left out

# load condition annotation file from Dave 

run <- function(file="./validation/gBg_stability.RData",type=c("runs","conditions")[1],n = 100, by = 10, rep = 10){
  doBackbone <- function(i,j,cmonkey_annotations,allBCs,to.r) {
    cat(paste(i,"-",j,"\n"))
    # select conditions
    runs.out = sample(seq(1,length(cmonkey_annotations$clusts)),replace=F,size=i)
    bc.remain = allBCs[-as.numeric(unlist(sapply(cmonkey_annotations$clusts[runs.out],function(x){strsplit(as.character(x),split=";")[[1]]})))]
    to.r$bcs[[j]] = bc.remain
    to.r$frac[[j]] = round(length(bc.remain)/length(allBCs),digits=2)
    if (length(to.r$bcs)>0) {
      gBg <- make.r.gBg(clusterStack=e$clusterStack[as.numeric(sapply(strsplit(to.r$bcs,split="_"),"[",2))])
      to.r$gBg[[j]] <- cor(as.vector(gBg),as.vector(o$gg$gBg))
      gBg.backbone <- multiscaleBackbone(gBg,multicore=MULTICORE)
      to.r$gBg.backbone[[j]] <- cor(as.vector(gBg.backbone),as.vector(o$gg$gBg.backbone))
      # write backbones to file
      write.table(gBg,file=paste("./validation/runStability_",i,"_gBg.txt",sep=""),col.names=T,row.names=T,sep="\t")
      write.table(gBg.backbone,file=paste("./validation/runStability_",i,"_gBg_backbone.txt",sep=""),col.names=T,row.names=T,sep="\t")
      cat(paste(to.r$gBg[[j]],":",to.r$gBg.backbone[[j]],"\n"))
    }
  }
	root = "/home/abrooks/Documents/EGRIN2/Eco_ensemble_2"
	BACKBONE.PVAL=0.05
	FULLY.CONNECTED=F 
	MULTICORE=F
	CORES=1
	RESID.FILTER=1
	MINMOTIFSCORE=Inf
	HAL=F
	# load EGRIN2 environment
	source(".Rprofile_backup")
	cmonkey_annotations <- read.table("./validation/info_for_aaron_cross_validation.tsv",header=T)
	# assess distribution of conditions
	# dist of conditions 
	conditions.per.run <- lapply(as.character(cmonkey_annotations[,3]),function(i){strsplit(i,split=";")[[1]][strsplit(i,split=";")[[1]]%in%colnames(e$ratios[[1]])]})
	all.conditions <- unique(unlist(conditions.per.run))
	dist.conditions <- sort(table(unlist(conditions.per.run)),decreasing=T)
	cond2bc <- lapply(all.conditions,function(i){
		is.it=sapply(conditions.per.run,function(j)i%in%j)
		return(which(is.it))
		})
	names(cond2bc) <- all.conditions
	allBCs = paste("BIC",seq(1,length(e$clusterStack),1),sep="_")
  require(multicore)
  if (type == "conditions") {
    after.condition.removal <- lapply(seq(1,n,by),function(i){
      cat(paste(i,"\n"))
      if (i>1) {
        i = i-1
      }
      to.r = sapply(seq(1,rep,1),function(j){
        # select conditions
        conds.out = sample(colnames(e$ratios[[1]]),replace=F,size=i)
        bc.remain = setdiff(allBCs,unique(unlist(out$get.biclusters(conditions=conds.out))))
        r = list()
        r$bcs = bc.remain
        r$frac = length(bc.remain)/length(allBCs)
        if (length(r$bcs)>0) {
          r$gBg <- make.r.gBg(clusterStack=e$clusterStack[as.numeric(sapply(strsplit(r$bcs,split="_"),"[",2))])
          r$gBg.cor <- cor(as.vector(r$gBg),as.vector(o$gg$gBg))
          #r$gBg.w <- wilcox.test(as.vector(r$gBg),as.vector(o$gg$gBg))
          r$gBg.backbone <- multiscaleBackbone(r$gBg)
          r$gBg.backbone.cor <- cor(as.vector(r$gBg.backbone),as.vector(o$gg$gBg))
          #r$gBg.backbone.w <- wilcox.test(as.vector(r$gBg.backbone),as.vector(o$gg$gBg.backbone))
        }
        })
      return(to.r)
    })
    names(after.condition.removal) = sapply(seq(1,n,by),function(i)if(i>1){return(i-1)}else{return(i)})
    save(after.condition.removal,file=file)
  } else if (type == "runs") {
    # remove some number of runs, see how this affects backbone
    
    after.run.removal <- lapply(seq(1,n,by),function(i){
      if (i>1) {
        i = i-1
      }
      to.r <- list()
      to.r$gBg <- list()
      to.r$gBg.backbone <- list()
      to.r$bcs = list()
      to.r$frac = list()
      j = 1
      while (j<=rep) {
        tmp = try(doBackbone(i,j,cmonkey_annotations,allBCs,to.r))
        if (class(tmp)!="try-error") {
          j = j+1
        } 
      }
      save(to.r,file=paste("./validation/runStability_",i,".RData",sep=""))
      return(to.r)
    })
    names(after.run.removal) = sapply(seq(1,n,by),function(i)if(i>1){return(i-1)}else{return(i)})
    save(after.run.removal,file=file)
  }
}



