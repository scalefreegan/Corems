# To address reviewer concerns for Brooks, Reiss et al 2014 MSB regarding model overfitting.

# Basic idea: resample backbone and corems when some of the runs (with specific conditions) are left out

# load condition annotation file from Dave 

run <- function(file="./validation/gBg_stability.RData",type=c("runs","conditions")[1],n = 104, by = 10, rep = 5,multicore=F,s=1){
  doBackbone <- function(i,j,cmonkey_annotations,allBCs) {
    cat(paste(i,"-",j,"\n"))
    # select conditions
    runs.out = sample(seq(1,length(cmonkey_annotations$clusts)),replace=F,size=i)
    bc.remain = allBCs[-as.numeric(unlist(sapply(cmonkey_annotations$clusts[runs.out],function(x){strsplit(as.character(x),split=";")[[1]]})))]
    to.r <- list()
    to.r$bcs = bc.remain
    to.r$frac = round(length(bc.remain)/length(allBCs),digits=2)
    cat(paste(to.r$frac,"\n"))
    if (length(to.r$bcs)>0) {
      gBg <- make.r.gBg(clusterStack=e$clusterStack[as.numeric(sapply(strsplit(to.r$bcs,split="_"),"[",2))])
      to.r$gBg <- cor(as.vector(gBg),as.vector(o$gg$gBg))
      gBg.backbone <- multiscaleBackbone(gBg,multicore=F)
      to.r$gBg.backbone <- cor(as.vector(gBg.backbone),as.vector(o$gg$gBg.backbone))
      # write backbones to file
      write.table(gBg,file=paste("./validation/runStability/runStability_",i,"_gBg.txt",sep=""),col.names=T,row.names=T,sep="\t")
      write.table(gBg.backbone,file=paste("./validation/runStability/runStability_",i,"_gBg_backbone.txt",sep=""),col.names=T,row.names=T,sep="\t")
      cat(paste(to.r$gBg,":",to.r$gBg.backbone,"\n"))
      return(to.r)
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
	#allBCs = paste("BIC",seq(1,length(e$clusterStack),1),sep="_")
  load("./validation/allBCs.RData")
  require(multicore)
  source("./scripts/main.R")
  source("./scripts/processEGRIN.R")
  source("./scripts/analyze_corems.R")
  o <- loadEnv()
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
          r$gBg.backbone.cor <- cor(as.vector(r$gBg.backbone),as.vector(o$gg$gBg.backbone))
          #r$gBg.backbone.w <- wilcox.test(as.vector(r$gBg.backbone),as.vector(o$gg$gBg.backbone))
        }
        })
      return(to.r)
    })
    names(after.condition.removal) = sapply(seq(1,n,by),function(i)if(i>1){return(i-1)}else{return(i)})
    save(after.condition.removal,file=file)
  } else if (type == "runs") {
    # remove some number of runs, see how this affects backbone
    
    after.run.removal <- lapply(c(n+1,rev(seq(s,n,by))),function(i){
      if (i>1) {
        i = i-1
      }
      j = 1
      to.r <- list()
      while (j<=rep) {
        tmp = try(doBackbone(i,j,cmonkey_annotations,allBCs))
        if (class(tmp)!="try-error") {
          to.r[[as.character(i)]][[as.character(j)]] <- tmp
          j = j+1
        }  
      }
      save(to.r,file=paste("./validation/runStability/runStability_",i,".RData",sep=""))
      return(to.r)
    })
    names(after.run.removal) = sapply(seq(1,n,by),function(i)if(i>1){return(i-1)}else{return(i)})
    save(after.run.removal,file=file)
  }
}

compileResults <- function(dir="./validation/runStability/") {
  f = list.files(dir)
  n = sort(unique(as.numeric(unlist(sapply(f,function(i){strsplit(strsplit(i,"_")[[1]][2],"\\.")[[1]][1]})))))
  o <- lapply(n,function(i){
    load(paste(dir,"runStability_",i,".RData",sep=""))
    return(to.r[[1]])
    })
  names(o) <- n
  return(o)
}

plotResults <- function(data) {
  require(ggplot2)
  require(reshape2)
  x = as.numeric(unlist(lapply(names(data),function(i){sapply(data[[i]],function(j){return(i)})})))
  x_frac = as.numeric(unlist(lapply(names(data),function(i){sapply(data[[i]],function(j){return(j$frac)})})))
  y_gBg_cor = as.numeric(unlist(lapply(names(data),function(i){sapply(data[[i]],function(j){return(j$gBg)})})))
  y_gBg_backbone_cor = as.numeric(unlist(lapply(names(data),function(i){sapply(data[[i]],function(j){return(j$gBg.backbone)})})))
  # convert to data.frame
  df <- data.frame(x=x,x_frac=x_frac,y_gBg_cor=y_gBg_cor,y_gBg_backbone_cor=y_gBg_backbone_cor)
  df.2=melt(df,measure.vars=(c("y_gBg_cor","y_gBg_backbone_cor")))
    c <- ggplot(df.2,aes(factor(round(1-x/105,2)),value,fill=variable,color=variable))  
    p1<-c+geom_boxplot()+labs(x="Fraction of runs included",y="Network correlation")+scale_fill_discrete(name="Network",
                         breaks=c("y_gBg_cor", "y_gBg_backbone_cor"),
                         labels=c("gBg", "gBg backbone"))+expand_limits(y=c(0,1))
    ggsave("./validation/runStability/stability.pdf", plot = p1)
}


