
make_roc_plot <- function(data, title) {
  data <- data.frame(data)
  final.data.melted <- melt(as.data.frame(data), variable.name = "class", value.name = "TPR", id.vars = colnames(data)[1])
  print(final.data.melted[1:10,])
  p <- ggplot(final.data.melted, aes(x = FPR, y = TPR, colour= class)) + geom_line() #+ geom_hline(yintercept = 1.0, color = "blue", linetype = "longdash") #+ geom_vline(xintercept=50, color = "blue", linetype = "longdash")
  p <- p + ylab("True Positive Rate") + xlab("False Positive Rate") + ggtitle(title)
  p <- p + geom_abline(intercept = 0.0, slope = 1)
  #p <- p +  geom_line(data=final.data.melted,thickness=3) #+ ylab(colnames(data)[2]) + xlab(colnames(data)[1])
  #p <- p + theme(axis.title = element_text(size = 12,face="bold"))
  #p <- p + theme(plot.title = element_text(size=12,face="bold"))
  #p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #               panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #p <- p + scale_x_continuous(expand=c(.002,0))
  #p <- p + theme(aspect.ratio=.5)
  print(p)
  #pdf(paste(title,".pdf"),width=12,height=6)
  #plot(p)
  #dev.off()
  return(p)
}

get_cluster_average <- function(cluster_name, df) {
  cluster_rows <- grep(paste(cluster_name, "_", sep=""), rownames(df), fixed=T)
  cluster_df <- as.data.frame(df[cluster_rows,])
  cluster_averages <- apply(cluster_df, 2, mean)
  if(cluster_name == "cluster0") {
    print(cluster_df)
    print(cluster_averages)
  }
  return(cluster_averages)
}


cluster_averages <- function(df) {
  nrows <- dim(df)[1]
  splitrows <- do.call(rbind, strsplit(rownames(df), "_"))
  cluster.names <- unique(splitrows[,1])
  df_averages <- as.data.frame(do.call(rbind,lapply(cluster.names, get_cluster_average, df)))
  rownames(df_averages) <- cluster.names
  return(df_averages)
}

get.cluster.percentile <- function(cluster_name, df, percentile) {
  cluster_rows <- grep(paste(cluster_name, "_", sep=""), rownames(df), fixed=T)
  cluster_df <- as.data.frame(df[cluster_rows,])
  cluster_averages <- apply(cluster_df, 2, quantile, probs=c(percentile))
  if(cluster_name == "cluster0") {
    print(cluster_df)
    print(cluster_averages)
  } 
  return(cluster_averages)
}

cluster.percentiles <- function(df, percentile) {
  nrows <- dim(df)[1]
  splitrows <- do.call(rbind, strsplit(rownames(df), "_"))
  cluster.names <- unique(splitrows[,1])
  df_averages <- as.data.frame(do.call(rbind,lapply(cluster.names, get.cluster.percentile, df, percentile)))
  rownames(df_averages) <- cluster.names
  return(df_averages)
}

cluster.averages.new <- function(df) {
  nrows <- dim(df)[1]
  splitrows <- do.call(rbind, strsplit(rownames(df), "_"))
  cluster.names <- unique(splitrows[,1])
  n.clusters <- length(cluster.names)
  df.averages <- matrix(0.0,nrow=n.clusters,ncol=dim(df)[2])
  rownames(df.averages) <- cluster.names
  for (i in 1:n.clusters) {
    cluster <- cluster.names[i]
    rows <- grep(cluster, rownames(df),fixed=T)
    if(i == 1) {
      print(rows)
    }
    if(length(rows) < 2) {
      print(dim(df.averages))
      print("1")
      print(df[rows,])
      print(dim(df[rows,]))
      df.averages[cluster,] <- unlist(df[rows,])
      print(dim(df.averages))
    } else {
      print("else")
      print(rows)
      cluster.averages <- apply(df[rows,],2,mean)
      print(cluster.averages)
      df.averages[cluster,] <- cluster.averages
    }
    
  }
  colnames(df.averages) <- colnames(df)
  return(df.averages)
}

cluster.percentiles.new <- function(df, percentile) {
  nrows <- dim(df)[1]
  splitrows <- do.call(rbind, strsplit(rownames(df), "_"))
  cluster.names <- unique(splitrows[,1])
  n.clusters <- length(cluster.names)
  df.averages <- matrix(0.0,nrow=n.clusters,ncol=dim(df)[2])
  rownames(df.averages) <- cluster.names
  for (i in 1:n.clusters) {
    cluster <- cluster.names[i]
    rows <- grep(cluster, rownames(df),fixed=T)
    if(i == 1) {
      print(rows)
    }
    print(rows)
    if(length(rows) < 2) {
      cluster.averages <- df[rows,]
    } else {
      cluster.averages <- apply(df[rows,],2,quantile,c(percentile))
    }
    df.averages[cluster,] <- cluster.averages
  }
  colnames(df.averages) <- colnames(df)
  return(df.averages)
}

cluster.sums <- function(df, percentile) {
  nrows <- dim(df)[1]
  splitrows <- do.call(rbind, strsplit(rownames(df), "_"))
  cluster.names <- unique(splitrows[,1])
  n.clusters <- length(cluster.names)
  df.averages <- matrix(0.0,nrow=n.clusters,ncol=dim(df)[2])
  rownames(df.averages) <- cluster.names
  for (i in 1:n.clusters) {
    cluster <- cluster.names[i]
    rows <- grep(cluster, rownames(df),fixed=T)
    if(i == 1) {
      print(rows)
    }
    cluster.averages <- apply(df[rows,],2, sum)
    df.averages[cluster,] <- cluster.averages
  }
  colnames(df.averages) <- colnames(df)
  return(df.averages)
}

which.set <- function(value, mylist) {
  is.in <- sapply(mylist, function(x, value) is.element(value, x), value, simplify=T)
  return(names(is.in[is.in==T]))
}

combine.columns <- function(column.1, column.2) {
  column.1.2 <- as.matrix(cbind(column.1, column.2))
  column.1.2 <- t(apply(column.1.2, 1, sort))
  
  print(column.1.2[1:5,])
  column.1.2 <- apply(column.1.2, 1, function(x) paste(x[1], x[2], sep="::"))
  return(column.1.2)
}

interhelix.distance <- function(df, only.helices = F) {
  helix.residues <- list()
  helix.residues[["tm1"]] <- seq(1,60)
  helix.residues[["icl1"]] <- seq(61,66)
  helix.residues[["tm2"]] <- seq(67,96)
  helix.residues[["ecl1"]] <- seq(97,101)
  helix.residues[["tm3"]] <- seq(102,136)
  helix.residues[["icl2"]] <- seq(137,146)
  helix.residues[["tm4"]] <- seq(147,170)
  helix.residues[["ecl2"]] <- seq(171,196)
  helix.residues[["tm5"]] <- seq(197,229)
  helix.residues[["icl3"]] <- seq(230,266)
  helix.residues[["tm6"]] <- seq(267, 298)
  helix.residues[["ecl3"]] <- seq(299,304)
  helix.residues[["tm7"]] <- seq(305,329)
  helix.residues[["tm8"]] <- seq(330,400)
  
  residues.1 <- df[,dim(df)[2]-1]
  residues.2 <- df[,dim(df)[2]]
  
  print(residues.1[1:5])
  
  helices.1 <- residues.1
  helices.2 <- residues.2
  helices.1 <- sapply(residues.1,which.set, helix.residues)
  print(helices.1[1:5])
  helices.2 <- sapply(residues.2,which.set, helix.residues)
  helices.1.2 <- combine.columns(helices.1, helices.2)
  
  rownames(df) <- helices.1.2
  df <- df[,1:(dim(df)[2]-2)]
  
  if(only.helices == T) {
    loop.rows <- grep("tm", rownames(df))
    df <- df[loop.rows,]
  }
  
  return(df)
}

which.helix <- function(df, only.helices = F) {
  helix.residues <- list()
  helix.residues[["tm1"]] <- seq(1,60)
  helix.residues[["icl1"]] <- seq(61,66)
  helix.residues[["tm2"]] <- seq(67,96)
  helix.residues[["ecl1"]] <- seq(97,101)
  helix.residues[["tm3"]] <- seq(102,136)
  helix.residues[["icl2"]] <- seq(137,146)
  helix.residues[["tm4"]] <- seq(147,170)
  helix.residues[["ecl2"]] <- seq(171,196)
  helix.residues[["tm5"]] <- seq(197,229)
  helix.residues[["icl3"]] <- seq(230,266)
  helix.residues[["tm6"]] <- seq(267, 298)
  helix.residues[["ecl3"]] <- seq(299,304)
  helix.residues[["tm7"]] <- seq(305,329)
  helix.residues[["tm8"]] <- seq(330,400)
  
  residues <- rownames(df)
  print(residues[1:10])
  helices <- sapply(residues,which.set, helix.residues)
  print(helices[1:10])
  rownames(df) <- helices
  
  if(only.helices == T) {
    loop.rows <- grep("tm", rownames(df))
    df <- df[loop.rows,]
  }
  
  return(df)
}


sort.by.column <- function(df) {
  indices <- seq(1,dim(df)[2])
  sorted.list <- lapply(indices, function(x) return(df[order(-1.0*df[,x]),,drop=F]))
  #sorted.list <- apply(df, 2, function(x) return(df[order(-1.0*x),]))
  return(sorted.list)
}

compute.aggregate.docking <- function(data, inverse.agonists = c()) {
  means <- apply(data,2,mean)
  sds <- apply(data,2,sd)
  scores.minus.means <- t(t(data)-means)
  z.scores <- t(t(scores.minus.means)/sds)
  if (length(inverse.agonists) >= 1) {
    inverse.cols <- which(colnames(data) %in% inverse.agonists)
    if (length(inverse.cols) >= 1) {
      z.scores[,inverse.cols] <- -1.0 * z.scores[,inverse.cols]
    }
  }
  aggregate.z.scores <- apply(z.scores,1,mean)
  return(aggregate.z.scores)
}

is_active <- function(row) {
  if(row["tm6_tm3_dist"] > 12.0 &  row["connector_rmsd_active"] < 1.0 & row["npxxy_rmsd_active"] < .7 & row["connector_rmsd_active"] < row["connector_rmsd_inactive"] & row["npxxy_rmsd_inactive"] > 0.8) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

find_active <- function(df) {
  active_rows <- apply(df, 1, is_active)
  return(active_rows)
}

is.intermediate <- function(row) {
  if(row["tm6_tm3_dist"] > 12.0 &  row["connector_rmsd_active"] < 1.0 & row["npxxy_rmsd_active"] > .7) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

find.intermediate <- function(df) {
  active_rows <- apply(df, 1, is.intermediate)
  return(active_rows)
}

is.inactive <- function(row) {
  if(row["tm6_tm3_dist"] < 12.0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

find.inactive <- function(df) {
  active_rows <- apply(df, 1, is.inactive)
  return(active_rows)
}

total_active <- function(df, col) {
  num_active = length(df[df[,col] == T, col])
  return(num_active)
}

active_under_index <- function(index, df, col) {
  under_index = df[1:index,]
  num_active = length(under_index[under_index[,col] == T, col])
  return(num_active)
}

enrichment <- function(percent_active, percent_all) {
  return(percent_active/percent_all)
}

calc_tpr <- function(index, df, col) {
  TP <- active_under_index(index, df, col)
  TP_FN <- total_active(df, col)
  sensitivity <- TP / TP_FN
  return(sensitivity)
}

calc_fpr <- function(index, df, col) {
  under_index = df[1:index,]
  TP_FP <- dim(under_index)[1]
  TP <- active_under_index(index, df, col)
  FP <- TP_FP - TP
  over_index <- df[index:dim(df)[1],]
  over_index_over_cutoff <- over_index[over_index[col] == F,]
  TN <- dim(over_index_over_cutoff)[1]
  fpr <- FP / (TN + FP)
  return(fpr)
}

test_method <- function(docking, pnas.rows, title, save.dir) {
  common.names <- intersect(rownames(docking),names(pnas.rows))
  docking <- docking[common.names,]
  pnas.rows <- pnas.rows[common.names]
  print(docking[1:10])
  clusters.active.docking <- data.frame(docking,pnas.rows)
  print(clusters.active.docking[1:10,])
  #clusters.active.docking <- merge(docking, pnas_rows, by = "row.names", all = TRUE)
  #rownames(clusters.active.docking) <- clusters.active.docking[,1] 
  #clusters.active.docking <- clusters.active.docking[,2:dim(clusters.active.docking)[2]]
  colnames(clusters.active.docking) <-c("aggregate_docking_score", "is_active")
  
  clusters.active.docking.ordered <- clusters.active.docking[order(-1.0 * clusters.active.docking$aggregate_docking_score),]
  
  N <- dim(clusters.active.docking.ordered)[1]
  cutoff <- 0.20
  index.fpr.tpr <- cbind(seq(1, N, 1), rep(0, N), rep(0, N)) 
  colnames(index.fpr.tpr) <- c("Index", "FPR", "TPR")
  index.fpr.tpr[,2] <- sapply(index.fpr.tpr[,1], calc_fpr, df = clusters.active.docking.ordered, col=2)
  index.fpr.tpr[,3] <- sapply(index.fpr.tpr[,1], calc_tpr, df = clusters.active.docking.ordered, col=2)
  fpr.tpr <- index.fpr.tpr[,c(2,3)]
  colnames(fpr.tpr) <- colnames(index.fpr.tpr)[c(2,3)]
  plot <- make_roc_plot(fpr.tpr, title)
  pdf(paste(save.dir, "/", title, ".pdf", sep=""))
  print(plot)
  dev.off()
  
  pos.scores <- clusters.active.docking.ordered[clusters.active.docking.ordered[,2] == T, 1]
  neg.scores <- clusters.active.docking.ordered[clusters.active.docking.ordered[,2] == F, 1]
  auc <- mean(sample(pos.scores,1000,replace=T) > sample(neg.scores,1000,replace=T))
  print(auc)
}

plot.docking.vs.reference <- function(docking, reference.docking, save.dir) {
  for (i in 1:dim(docking)[2]) {
    ligand_name <- colnames(docking)[i]
    print(ligand_name)
    pdf(paste(save.dir, "/docking_", ligand_name, ".pdf", sep=""))
    hist(docking[,i], breaks=seq(0,15,.1), main = paste("Docking score ", ligand_name, sep=""), xlim = range(0,15))
    abline(v = -1.0*reference.docking[1,ligand_name], col="green")
    abline(v = -1.0*reference.docking[2, ligand_name], col = "blue")
    dev.off()
  }
}

combine.dfs <- function(df1, df2) {
  common.rows <- intersect(rownames(df1), rownames(df2))
  print(length(common.rows))
  print("rows in common")
  df1 <- df1[common.rows, colnames(df1), drop=F]
  df2 <- df2[common.rows, colnames(df2), drop=F]
  new.df <- data.frame(df1,df2)
  return(new.df)
}

plot.docking.vs.reaction.coord <- function(docking, pnas.coords, save.dir) {
  d <- combine.dfs(docking, pnas.coords)
  colnames(d)[dim(d)[2]] <- "aggregate.docking.score"
  d <- data.frame(log(1/pnas.coords.averages[,3]), docking.aggregated[rownames(pnas.coords.averages),1])
  colnames(d) <- c("log_inverse_npxxy_rmsd_active", "aggregate_docking_score")
  pdf(paste(save.dir, "/aggregate_docking_vs_npxxy_rmsd_active", sep=""))
  logEstimate <- lm(aggregate_docking_score ~ log_inverse_npxxy_rmsd_active,data=d)
  plot(d, xlab = "ln(1/npxxy_rmsd_active)", ylab = "aggregate docking score")
  abline(logEstimate)
  dev.off()
}



plot_coords_and_docking <- function(df, ori_df, refcoords, title, save.dir) {
  colors <- brewer.pal(9, "YlGnBu")
  brks<-classIntervals(ori_df[,3], n=9, style="fisher")
  brks <- brks$brks
  pdf(paste(save.dir, "/", title, ".pdf", sep=""))
  plot(df[,c(1,2)], col=colors[findInterval(df[,3], brks, all.inside=TRUE)], axes = T, pch=16, main = title)
  par(new = T)
  points(refcoords, col = "red", pch=18)
  par(new = F)
}

plot.colmap <- function(docking, pnas, refcoords, title, save.dir, top=0) {
  df <- combine.dfs(pnas, docking)
  print(df[1:10,])
  colors <- brewer.pal(9, "YlGnBu")
  brks <- classIntervals(df[,3], n=9, style="fisher")
  brks <- brks$brks
  pdf(paste(save.dir, "/", title, ".pdf", sep=""))
  plot(df[,c(1,2),drop=F], col = colors[findInterval(df[,3],brks,all.inside=TRUE)], axes = T , pch = 16, main=title)#, ylim=c(0.0,1.5))
  legend(x=6, y=0.5, legend=leglabs(round(brks)), fill=colors, bty="n",x.intersp = .5, y.intersp = 0.75)
  par(new = T)
  #points(refcoords, col = "red", pch=18)
  par(new = F)
  dev.off()
  
  if(top > 0) {
    df <- df[order(-1.0*df[,3]),]
    df <- df[1:top,]
    print("New data frame:")
    print(df[1:10,])
    new.title <- paste(title, " Top ", top, " Docking Scores", sep="")
    pdf(paste(save.dir, "/", new.title, ".pdf", sep=""))
    plot(df[,c(1,2),drop=F], col = colors[findInterval(df[,3],brks,all.inside=TRUE)], axes = T , pch = 16, main=new.title, ylim=c(0.0,1.5))
    par(new = T)
    points(refcoords, col = "red", pch=18)
    par(new = F)
    dev.off()
    
  }
  
}

convert.docking.to.rel.probs <- function(docking) {
  docking <- 2.718 ^ (docking/0.58)
  docking <- docking / max(docking)
  return(docking)
}

calc.ranks <- function(df) {
  num.unique.residues <- length(unique(df[1,]))
  ranks.per.tIC <- as.data.frame(matrix(0,nrow=num.unique.residues,ncol=dim(df)[2]))
  print(dim(ranks.per.tIC))
  rownames(ranks.per.tIC) <- sort(unique(df[1,]))
  
  ranks <- rep(1:(dim(df)[1]/2),each=2)
  
  for(i in 1:dim(df)[2]) {
    average.ranks <- aggregate(ranks,by=list(df[,i]),FUN = mean)
    average.ranks <- average.ranks[order(average.ranks[,1]),]
    print(average.ranks[1:10,])
    print(dim(average.ranks))
    print(dim(ranks.per.tIC))
    ranks.per.tIC[,i] <- average.ranks[,2]
  }
  return(ranks.per.tIC)
  
}

get.tic.residues <- function(tic.residue.ranks.csv = "", tic.residues.csv = "", tic.duplicated.residues.csv = "", analysis.dir = "", features.dir = "", tica = "") {
  #tic.residues.csv <- paste(tica, "/top_residues_per_tIC.csv", sep="")
  #tic.residues <- read.csv(tic.residues.csv,stringsAsFactors=T,sep=",",header=F)[,-1]
  #tic.residues <- t(tic.residues)
  #write.table(tic.residues, paste(tica, "/top_residues_per_tIC_r.csv", sep=""))
  
  #tic.residues.csv <- paste(analysis.dir, "/feature_coefs.csv", sep="")
  #tic.residues <- read.csv(tic.residues.csv,stringsAsFactors=T,sep=",",header=F)[,-1]
  #tic.residues <- t(tic.residues)
  #colnames(tic.residues) <- rep(NA, dim(tic.residues)[2])
  #for(i in 1:(dim(tic.residues)[2]-2)) {
  #  colnames(tic.residues)[i] <- paste("tIC.",i,sep="")
  #}
  #colnames(tic.residues)[dim(tic.residues)[2]-1] <- "residue.1"
  #colnames(tic.residues)[dim(tic.residues)[2]] <- "residue.2"
  #write.table(tic.residues, paste(tica, "/feature_coefs_r.csv", sep=""))
  
  tic.residues.csv <- paste(tica, "/duplicated_feature_coefs.csv", sep="")
  tic.residues <- read.csv(tic.residues.csv,stringsAsFactors=T,sep=",",header=F)[,-1]
  tic.residues <- t(tic.residues)
  colnames(tic.residues) <- rep(NA, dim(tic.residues)[2])
  for(i in 1:(dim(tic.residues)[2]-1)) {
    colnames(tic.residues)[i] <- paste("tIC.",i,sep="")
  }
  colnames(tic.residues)[dim(tic.residues)[2]] <- "residue"
  rownames(tic.residues) <- tic.residues[,"residue"]
  tic.residues <- tic.residues[,-dim(tic.residues)[2]]
  write.table(tic.residues, paste(tica, "/duplicated_feature_coefs_r.csv", sep=""))
  
  return(tic.residues)
}

get.tic.residues.pairs <- function(tic.residue.ranks.csv = "", tic.residues.csv = "", tic.duplicated.residues.csv = "", analysis.dir = "", features.dir = "", tica = "") {
  #tic.residues.csv <- paste(tica, "/top_residues_per_tIC.csv", sep="")
  #tic.residues <- read.csv(tic.residues.csv,stringsAsFactors=T,sep=",",header=F)[,-1]
  #tic.residues <- t(tic.residues)
  #write.table(tic.residues, paste(tica, "/top_residues_per_tIC_r.csv", sep=""))
  
  tic.residues.csv <- paste(tica, "/feature_coefs.csv", sep="")
  print(tic.residues.csv)
  tic.residues <- read.csv(tic.residues.csv,stringsAsFactors=T,sep=",",header=F)[,-1]
  tic.residues <- t(tic.residues)
  print(tic.residues[1:3,])
  print(dim(tic.residues))
  colnames(tic.residues) <- rep(NA, dim(tic.residues)[2])
  for(i in 1:(dim(tic.residues)[2]-2)) {
    colnames(tic.residues)[i] <- paste("tIC.",i,sep="")
  }
  colnames(tic.residues)[dim(tic.residues)[2]-1] <- "residue.1"
  colnames(tic.residues)[dim(tic.residues)[2]] <- "residue.2"
  write.table(tic.residues, paste(tica, "/feature_coefs_r.csv", sep=""))
  
  return(tic.residues)
}

