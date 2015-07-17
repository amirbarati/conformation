library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(classInt)
library(maptools)
inverse.agonists <- c("s.atenolol", "s.carazolol")

base <- "/Users/Evan/vsp/b2ar_analysis"
#base <- "/Users/Evan/vsp/b2ar_analysis/exacycle_data"
#tica <- "/tICA_t5_n_components10_skip5_switches_pp_npxx_contact"
tica <- "/tICA_t5_n_components10_skip5_switches_pp_npxx_contact/ktICA_random_specified"
#tica <- "/tICA_t10_n_components10_skip5_switches_pp_npxx_contact"
#tica <- "/tICA_t10_n_components5_switches_npxx_tm6_bp"
#tica <- "tICA_t10_n_components"
tica <- paste(base, tica, sep = "")
clusters <- 1000
exacycle <- ""
method <- "random"
analysis.dir <- paste(tica, "/analysis_n_clusters", clusters, "_", method, sep = "")

pnas.coords.csv <- paste(analysis.dir, "/pnas_coords_new.csv", sep = "")
tica.coords.csv <- paste(analysis.dir, "/tica_coords.csv", sep = "")
docking.csv <- paste(analysis.dir, "/all_docking_combined.csv", sep = "")
docking.aggregated.csv <- paste(analysis.dir, "/aggregate_docking.csv", sep = "")
sasa.csv <- paste(analysis.dir, "/sasa_bp.csv", sep = "")
#pnas.coords.all_csv <- "/Users/evan/scratch_enf/b2ar_analysis/exacycle_data/all_pnas_features/pnas_all_coords.csv"
#pnas.coords.all <- data.frame(read.csv(pnas.coords.all_csv, stringsAsFactors = F, row.names=1))
#badrows <- which(pnas.coords.all[,2] > 4.0 | pnas.coords.all[,3] > 4.0 | pnas.coords.all[,4] > 4.0 | pnas.coords.all[,5] > 4.0)
#if(length(badrows) > 0) {
#  pnas.coords.all <- pnas.coords.all[-badrows,]
#}
#colnames(pnas.coords.all) <- colnames(pnas.coords)
#pnas.coords.all["tm6_tm3_dist"] <- 7.14 * pnas.coords.all["tm6_tm3_dist"]


pnas.coords <- data.frame(read.csv(pnas.coords.csv, stringsAsFactors = F, row.names=1))[,c(1,2,3,4,5)]
pnas.coords["tm6_tm3_dist"] <- 7.14 * pnas.coords["tm6_tm3_dist"]
tica.coords <- data.frame(read.csv(tica.coords.csv, stringsAsFactors = F, row.names=1))
for (i in 1:dim(tica.coords)[2]) {
  colnames(tica.coords)[i] <- paste("tIC.", i, sep="")
}
badrows <- which(pnas.coords[,2] > 4.0 | pnas.coords[,3] > 4.0 | pnas.coords[,4] > 4.0 | pnas.coords[,5] > 4.0)
if(length(badrows) > 0) {
  pnas.coords <- pnas.coords[-badrows,]
  tica.coords <- tica.coords[-badrows,]
}
dim(pnas.coords)
pnas.coords[1:10,]


#colnames(pnas.coords.all) <- colnames(pnas.coords)
#pnas.coords.all["tm6_tm3_dist"] <- 7.14 * pnas.coords.all["tm6_tm3_dist"]

#docking <- data.frame(read.csv(docking.csv, stringsAsFactors = F, row.names=1))
#docking.aggregated <- data.frame(read.csv(docking.aggregated.csv, stringsAsFactors = F, row.names=1))
reference.docking <- data.frame(read.csv("/Users/Evan/vsp/b2ar_analysis/reference_docking/docking_SP/all_docking_combined_manual.csv", stringsAsFactors = F, row.names=1))

#sasa <- data.frame(read.csv(sasa.csv, stringsAsFactors = F, row.names=1))

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
    cluster.averages <- apply(df[rows,],2,mean)
    df.averages[cluster,] <- cluster.averages
  }
  return(df.averages)
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

plot.docking.vs.reaction.coord <- function(docking, pnas.coords, save.dir) {
  d <- data.frame(pnas.coords, docking)
  colnames(d)[dim(d)[2]] <- "aggregate.docking.score"
  d <- data.frame(log(1/pnas.coords.averages[,3]), docking.aggregated[rownames(pnas.coords.averages),1])
  colnames(d) <- c("log_inverse_npxxy_rmsd_active", "aggregate_docking_score")
  pdf(paste(save.dir, "/aggregate_docking_vs_npxxy_rmsd_active", sep=""))
  logEstimate <- lm(aggregate_docking_score ~ log_inverse_npxxy_rmsd_active,data=d)
  plot(d, xlab = "ln(1/npxxy_rmsd_active)", ylab = "aggregate docking score")
  abline(logEstimate)
  dev.off()
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
  plot(df[,c(1,2),drop=F], col = colors[findInterval(df[,3],brks,all.inside=TRUE)], axes = T , pch = 16, main=title, ylim=c(0.0,1.5))
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

refcoords <- read.csv(file = "/Users/evan/vsp/b2ar_analysis/reference_receptors/ref_coords.csv", row.names = 1)
colnames(refcoords) <- colnames(pnas.coords)
plot_coords_and_docking(pnas.coords.docking.sorted[1:50,], pnas.coords.docking.sorted, refcoords)

#plot.docking.vs.reaction.coord(docking.aggregated, pnas.coords.averages, analysis.dir)

#all_pnas_rows <- find_active(pnas.coords.all)
#all_active_rows <- all_pnas_rows[all_pnas_rows == T]

pnas.coords.averages <- cluster_averages(pnas.coords)
tica.coords.averages <- cluster_averages(tica.coords)

pnas.tica <- combine.dfs(pnas.coords.averages, tica.coords.averages)
pnas.tica <- pnas.tica[,c(1,3,5,6,7,8,9,10,11,12,13,14,15)]
pairs(pnas.tica)
pnas_1_tIC_0 <- lm(tIC_0~npxxy_rmsd_active, data=pnas.tica)
summary(pnas_1_tIC_0)
plot(pnas.tica$npxxy_rmsd_active, pnas.tica$tIC_0)
abline(pnas_1_tIC_0)

#docking.subset <- docking
#docking.subset <- docking[,which(colnames(docking) %in% inverse.agonists)]
#docking.subset <- docking[,c("X3p0g_lig"),drop=F]
#aggregate.docking <- as.data.frame(compute.aggregate.docking(docking.subset, inverse.agonists))
#docking.averages <- cluster_averages(aggregate.docking)
#test_method(docking.averages, pnas.rows, "Aggregate Docking Score - Co-Crystallized 3P0G Agonist", analysis.dir)
#plot.docking.vs.reference(docking.averages, reference.docking, analysis.dir)

coords <- c("tm6_tm3_dist", "npxxy_rmsd_active")
plot.colmap(docking.averages, pnas.coords.averages[,coords, drop=F], refcoords[,coords], "Docking Score All Clusters vs TM6_TM3 dist and NPxxY RMSD to Inactive", analysis.dir, top=100)
plot.colmap(tica.coords.averages[,9,drop=F],pnas.coords.averages[,coords,drop=F], refcoords[,coords], "A priori Reaction Coordinates vs. tIC 9", analysis.dir)


#pnas.rows <- find_active(pnas.coords)
active.rows <- find_active(pnas.coords.averages)
write.table(t(as.data.frame(names(active.rows[active.rows==T]))), file = paste(analysis.dir, "/", "active_clusters.csv", sep=""), row.names = F, col.names = F, sep = ",")

intermediate.rows <- find.intermediate(pnas.coords.averages)
write.table(t(as.data.frame(names(intermediate.rows[intermediate.rows==T]))), file = paste(analysis.dir, "/", "intermediate_clusters.csv", sep=""), row.names = F, col.names = F, sep = ",")

inactive.rows <- find.inactive(pnas.coords.averages)
write.table(t(as.data.frame(names(inactive.rows[inactive.rows==T]))), file = paste(analysis.dir, "/", "inactive_clusters.csv", sep=""), row.names = F, col.names = F, sep = ",")
#docking.aggregated <- as.data.frame(apply(docking.averages[,c(1,2,3,4),drop=F],1,mean))

#test_method(docking.averages, pnas.rows)
active.rows <- pnas.rows[pnas.rows == T]
intermediate.rows <- 
inactive.rows <- pnas.rows[pnas.rows == F]

#all_pnas_rows <- find_active(pnas.coords.all)
#all_active_rows <- all_pnas_rows[all_pnas_rows == T]

#docking.csv <- "/Users/evan/biox3/b2ar_analysis_sherlock_all/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking_joined.csv"



#0.779 for cutoff of 0.3
#0.818 for cutoff of 0.25
#0.948 for cutoff of 0.2

#print(length(active_rows))
#print(length(active_rows)/dim(pnas.coords.averages)[1])
#print(length(all_active_rows)/(length(all_pnas_rows)))



#pnas_clusters_rows <- find_active(pnas.coords.averages)
#active_clusters_rows <- pnas_clusters_rows[pnas_clusters_rows == T]
#inactive_clusters_rows <- pnas_clusters_rows[pnas_clusters_rows == F]

#pnas_rows <- find_active(pnas.coords)
#active_rows <- pnas_rows[pnas_rows == T]
#inactive_rows <- pnas_rows[pnas_rows == F]

#print(length(active_rows))
#print(length(active_rows)/dim(pnas.coords)[1])
#print(length(all_active_rows)/(length(all_pnas_rows)))

#print(length(active_clusters_rows))
#print(length(active_clusters_rows)/dim(pnas.coords.averages)[1])
#print(length(all_active_rows)/(length(all_pnas_rows)))

#hist(pnas.coords[,1], breaks=20)
#hist(tica.coords[,2], breaks=20)
