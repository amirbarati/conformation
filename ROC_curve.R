library(ggplot2)
library(reshape2)


make_roc_plot <- function(data, title) {
  data <- data.frame(data)[1:1000,]
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

pdb_lig_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/docking_n_clusters1000_n_samples10_dist_SP/3p0g_lig/docking_means_pnas_means.csv"
pdb_lig_pnas_active <- data.frame(read.csv(pdb_lig_csv))
pdb_lig_sorted <- pdb_lig_pnas_active[order(-pdb_lig_pnas_active[,2]),]

all_pnas_active_csv <- "/Users/evan/scratch_enf/b2ar_analysis/features_pnas/all_pnas_distances.csv"
all_pnas_active <- data.frame(read.csv(all_pnas_active_csv))


aggregate_docking_pnas_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking_pnas_joined.csv"
aggregate_docking_pnas <- data.frame(read.csv(aggregate_docking_pnas_csv))
aggregate_docking_sorted <- aggregate_docking_pnas[order(-aggregate_docking_pnas[,2]),]

total_active <- function(df, cutoff, col) {
  num_active = length(df[df[,col] < cutoff, col])
  return(num_active)
}

active_under_index <- function(index, df, cutoff, col) {
  under_index = df[1:index,]
  num_active = length(under_index[under_index[,col] < cutoff, col])
  return(num_active)
}

enrichment <- function(percent_active, percent_all) {
  return(percent_active/percent_all)
}

calc_tpr <- function(index, df, cutoff, col) {
  TP <- active_under_index(index, df, cutoff, col)
  TP_FN <- total_active(df, cutoff, col)
  sensitivity <- TP / TP_FN
  return(sensitivity)
}

calc_fpr <- function(index, df, cutoff, col) {
  under_index = df[1:index,]
  TP_FP <- dim(under_index)[1]
  TP <- active_under_index(index, df, cutoff, col)
  FP <- TP_FP - TP
  over_index <- df[index:dim(df)[1],]
  over_index_over_cutoff <- over_index[over_index[col] > cutoff,]
  TN <- dim(over_index_over_cutoff)[1]
  fpr <- FP / (TN + FP)
  return(fpr)
}

N <- dim(aggregate_docking_sorted)[1]
cutoff <- 0.20
index.fpr.tpr <- cbind(seq(1, N, 1), rep(0, N), rep(0, N))
colnames(index.fpr.tpr) <- c("Index", "FPR", "TPR")
index.fpr.tpr[,2] <- sapply(index.fpr.tpr[,1], calc_fpr, df = aggregate_docking_sorted, cutoff = cutoff, col=3)
index.fpr.tpr[,3] <- sapply(index.fpr.tpr[,1], calc_tpr, df = aggregate_docking_sorted, cutoff = cutoff, col=3)
fpr.tpr <- index.fpr.tpr[,c(2,3)]
colnames(fpr.tpr) <- colnames(index.fpr.tpr)[c(2,3)]
make_roc_plot(fpr.tpr, paste("ROC Curve Plot, PNAS Score cutoff = ", cutoff))

pos.scores <- aggregate_docking_sorted[aggregate_docking_sorted[,3] < cutoff, 2]
neg.scores <- aggregate_docking_sorted[aggregate_docking_sorted[,3] >= cutoff, 2]

auc <- mean(sample(pos.scores,1000,replace=T) > sample(neg.scores,1000,replace=T))
#0.779 for cutoff of 0.3
#0.818 for cutoff of 0.25
#0.948 for cutoff of 0.2