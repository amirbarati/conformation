library(ggplot2)
library(reshape2)


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
  p
  #pdf(paste(title,".pdf"),width=12,height=6)
  #plot(p)
  #dev.off()
}

get_cluster_average <- function(cluster_name, df) {
  cluster_rows <- grep(paste(cluster_name, "_", sep=""), rownames(df))
  cluster_df <- df[cluster_rows,]
  if(is.vector(cluster_df)) {
    cluster_averages <- median(cluster_df)
  } else {
    cluster_averages <- apply(cluster_df, 2, median)
  }
  
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


is_active <- function(row) {
  if(row["tm6_tm3_dist"] > 12.0 &  row["connector_rmsd_active"] < 1.0 & row["npxxy_rmsd_active"] < .7 & row["connector_rmsd_active"] < row["connector_rmsd_inactive"] & row["npxxy_rmsd_inactive"] > 0.8) {    return(TRUE)
  } else {
    return(FALSE)
  }
}

find_active <- function(df) {
  active_rows <- apply(df, 1, is_active)
  return(active_rows)
}

#pnas_coords_csv <- "/Users/evan/Downloads/biox3_home/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"
pnas_coords_csv <- "/Users/evan/Downloads/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_skip5_switches_pp_npxx_contact/analysis_n_clusters3000_random/pnas_coords_new.csv"
pnas_coords <- data.frame(read.csv(pnas_coords_csv, stringsAsFactors = F, row.names=1))
pnas_coords["tm6_tm3_dist"] <- 7.14 * pnas_coords["tm6_tm3_dist"]
badrows <- which(pnas_coords[,2] > 4.0 | pnas_coords[,3] > 4.0 | pnas_coords[,4] > 4.0 | pnas_coords[,5] > 4.0)
if(length(badrows) > 0) {
  pnas_coords <- pnas_coords[-badrows,]
}

dim(pnas_coords)
pnas_coords[1:10,]

#pnas_coords_all_csv <- "/Users/evan/biox3/b2ar_analysis_sherlock_all/b2ar_analysis/all_pnas_features/pnas_all_coords.csv"
#pnas_coords_all_csv <- "/Users/evan/Downloads/vsp/b2ar_analysis/exacycle_data/all_pnas_features/pnas_all_coords.csv"
#pnas_coords_all <- data.frame(read.csv(pnas_coords_all_csv, stringsAsFactors = F, row.names=1))
#badrows <- which(pnas_coords_all[,2] > 4.0 | pnas_coords_all[,3] > 4.0 | pnas_coords_all[,4] > 4.0 | pnas_coords_all[,5] > 4.0)
#if(length(badrows) > 0) {
#  pnas_coords_all <- pnas_coords_all[-badrows,]
#}

#colnames(pnas_coords_all) <- colnames(pnas_coords)
#pnas_coords_all["tm6_tm3_dist"] <- 7.14 * pnas_coords_all["tm6_tm3_dist"]

pnas_coords_averages <- cluster_averages(pnas_coords)

pnas_rows <- find_active(pnas_coords_averages)
active_rows <- pnas_rows[pnas_rows == T]
inactive_rows <- pnas_rows[pnas_rows == F]

#all_pnas_rows <- find_active(pnas_coords_all)
#all_active_rows <- all_pnas_rows[all_pnas_rows == T]

docking_csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/docking_n_clusters1000_n_samples10_dist_SP/s-carazololl/docking_summary.csv"
#docking_csv <- "/Users/evan/Downloads/biox3_home/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/docking_n_clusters1000_n_samples10_dist_SP/docking_summary.csv"
#docking_csv <- "/Users/evan/Downloads/biox3_home/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking_joined.csv"
#docking_csv <- "/Users/evan/Downloads/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_skip5_switches_pp_npxx_contact/analysis_n_clusters3000_random/aggregate_docking_joined.csv"
#docking_csv <- "/Users/evan/Downloads/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_skip5_switches_pp_npxx_contact/docking_n_clusters3000_n_samples20_random_SP/3p0g_lig/docking_summary.csv"
docking_original <- data.frame(read.csv(docking_csv, stringsAsFactors = F, row.names=1))
docking <- cluster_averages(docking_original)

clusters.active.docking <- merge(docking, pnas_rows, by = "row.names", all = TRUE)
rownames(clusters.active.docking) <- clusters.active.docking[,1] 
clusters.active.docking <- clusters.active.docking[,2:dim(clusters.active.docking)[2]]
colnames(clusters.active.docking) <-c("aggregate_docking_score", "is_active")

clusters.active.docking.ordered <- clusters.active.docking[order(-1.0 * clusters.active.docking$aggregate_docking_score),]
good_tm6_rows <- rownames(pnas_coords_averages)[which(pnas_coords_averages$tm6_tm3_dist > 12.0)]
clusters.active.docking.ordered <- clusters.active.docking.ordered[good_tm6_rows,]

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

N <- dim(clusters.active.docking.ordered)[1]
cutoff <- 0.20
index.fpr.tpr <- cbind(seq(1, N, 1), rep(0, N), rep(0, N)) 
colnames(index.fpr.tpr) <- c("Index", "FPR", "TPR")
index.fpr.tpr[,2] <- sapply(index.fpr.tpr[,1], calc_fpr, df = clusters.active.docking.ordered, col=2)
index.fpr.tpr[,3] <- sapply(index.fpr.tpr[,1], calc_tpr, df = clusters.active.docking.ordered, col=2)
fpr.tpr <- index.fpr.tpr[,c(2,3)]
colnames(fpr.tpr) <- colnames(index.fpr.tpr)[c(2,3)]
make_roc_plot(fpr.tpr, paste("ROC Curve Plot"))

pos.scores <- clusters.active.docking.ordered[clusters.active.docking.ordered[,2] == T, 1]
neg.scores <- clusters.active.docking.ordered[clusters.active.docking.ordered[,2] == F, 1]
auc <- mean(sample(pos.scores,1000,replace=T) > sample(neg.scores,1000,replace=T))

#0.779 for cutoff of 0.3
#0.818 for cutoff of 0.25
#0.948 for cutoff of 0.2

print(length(active_rows))
print(length(active_rows)/dim(pnas_coords_averages)[1])
#print(length(all_active_rows)/(length(all_pnas_rows)))
