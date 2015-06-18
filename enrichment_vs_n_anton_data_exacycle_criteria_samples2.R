library(ggplot2)
library(reshape2)


make_roc_plot <- function(data, title) {
  data <- data.frame(data)[1:10000,]
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

#pnas_coords_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"
pnas_coords_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"
pnas_coords <- data.frame(read.csv(pnas_coords_csv, stringsAsFactors = F, row.names=1))
pnas_coords["tm6_tm3_dist"] <- 7.14 * pnas_coords["tm6_tm3_dist"]
pnas_coords[1:10,]

pnas_coords_all_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/all_pnas_features/pnas_all_coords.csv"
pnas_coords_all <- data.frame(read.csv(pnas_coords_all_csv, stringsAsFactors = F, row.names=1))
colnames(pnas_coords_all) <- colnames(pnas_coords)
pnas_coords_all["tm6_tm3_dist"] <- 7.14 * pnas_coords_all["tm6_tm3_dist"]


is_active <- function(row) {
  #print(row)
  #print(rownames(row))
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

is_active(pnas_coords[1,])
pnas_rows <- find_active(pnas_coords)
active_rows <- pnas_rows[pnas_rows == T]
inactive_rows <- pnas_rows[pnas_rows == F]
print(length(active_rows))
print(length(active_rows)/dim(pnas_coords)[1])

#all_rows <- find_active(pnas_coords_all)
#all_active_rows <- all_rows[all_rows == T]
#print(length(all_active_rows))
#print(length(all_active_rows)/length(all_rows))

#docking_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking.csv"
docking_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking.csv"

docking <- data.frame(read.csv(docking_csv, stringsAsFactors = F, row.names=1))

samples.active.docking <- merge(docking, pnas_rows, by = "row.names", all = TRUE)
rownames(samples.active.docking) <- samples.active.docking[,1] 
samples.active.docking <- samples.active.docking[,2:dim(samples.active.docking)[2]]
colnames(samples.active.docking) <-c("aggregate_docking_score", "is_active")

samples.active.docking.ordered <- samples.active.docking[order(-1.0 * samples.active.docking$aggregate_docking_score),]


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

N <- dim(samples.active.docking.ordered)[1]
cutoff <- 0.20
index.fpr.tpr <- cbind(seq(1, N, 1), rep(0, N), rep(0, N)) 
colnames(index.fpr.tpr) <- c("Index", "FPR", "TPR")
index.fpr.tpr[,2] <- sapply(index.fpr.tpr[,1], calc_fpr, df = samples.active.docking.ordered, col=2)
index.fpr.tpr[,3] <- sapply(index.fpr.tpr[,1], calc_tpr, df = samples.active.docking.ordered, col=2)
fpr.tpr <- index.fpr.tpr[,c(2,3)]
colnames(fpr.tpr) <- colnames(index.fpr.tpr)[c(2,3)]
make_roc_plot(fpr.tpr, paste("ROC Curve Plot"))

pos.scores <- samples.active.docking.ordered[samples.active.docking.ordered[,2] == T, 1]
neg.scores <- samples.active.docking.ordered[samples.active.docking.ordered[,2] == F, 1]
auc <- mean(sample(pos.scores,10000,replace=T) > sample(neg.scores,10000,replace=T))

#0.779 for cutoff of 0.3
#0.818 for cutoff of 0.25
#0.948 for cutoff of 0.2
