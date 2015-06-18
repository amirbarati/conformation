library(ggplot2)
library(reshape2)

pdb_lig_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/docking_n_clusters1000_n_samples10_dist_SP/3p0g_lig/docking_means_pnas_means.csv"
pdb_lig_pnas_active <- data.frame(read.csv(pdb_lig_csv))
pdb_lig_sorted <- pdb_lig_pnas_active[order(-pdb_lig_pnas_active[,2]),]

all_pnas_active_csv <- "/Users/evan/scratch_enf/b2ar_analysis/features_pnas/all_pnas_distances.csv"
all_pnas_active <- data.frame(read.csv(all_pnas_active_csv))


aggregate_docking_pnas_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking_pnas_joined.csv"
aggregate_docking_pnas <- data.frame(read.csv(aggregate_docking_pnas_csv))
aggregate_docking_sorted <- aggregate_docking_pnas[order(-aggregate_docking_pnas[,2]),]

make_plot <- function(data, title) {
  data <- data.frame(data)[1:1000,]
  final.data.melted <- melt(as.data.frame(data), id.vars = colnames(data)[1])
  p <- ggplot(final.data.melted, aes(x = Index, y = value, colour= variable)) + geom_line() #+ geom_hline(yintercept = 1.0, color = "blue", linetype = "longdash") #+ geom_vline(xintercept=50, color = "blue", linetype = "longdash")
  p <- p + ylab("Accuracy") + xlab("Top 'n' clusters scored by docking") + ggtitle(title)
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

total_active <- function(df, cutoff, col) {
  num_active = length(df[df[,col] < cutoff, col])
  print(num_active)
  print(length(df[,1]))
  percent_active = num_active / length(df[,1])
  return num_active
}

active_under_index <- function(index, df, cutoff, col) {
  under_index = df[1:index,]
  num_active = length(under_index[under_index[,col] < cutoff, col])
  print(num_active)
  percent_active = num_active / length(under_index[,1])
  return(num_active)
}

enrichment <- function(percent_active, percent_all) {
  return(percent_active/percent_all)
}

accuracy <- function(percent_active)

cutoff <- 0.15

total_percent_active_1000_clusters <- total_active(pdb_lig_pnas_active, cutoff, 3)
#multiple scoring functions
active_vs_index <- cbind(1:1000, rep(0,1000))
active_vs_cutoff <- cbind(seq(0.1, 0.5, .05), rep(0,9))

active_vs_index[,2] <- sapply(enrichment_vs_index[,1], active_under_index, df = aggregate_docking_sorted, cutoff = cutoff, col = 3)  
  
enrichment_vs_index <- active_vs_index
enrichment_vs_index[,2] <- (enrichment_vs_index[,2] / enrichment_vs_index[,1]) / (total_percent_active_1000_clusters)

accuracy_vs_index <- active_vs_index
accuracy_vs_index[,2] <- active_vs_index[,2]/active_vs_index[,1]
plot(accuracy_vs_index)

#single scoring function w/ 3p0g ligand

active_vs_index_single <- cbind(1:1000, rep(0,1000))
active_vs_cutoff_single <- cbind(seq(0.1, 0.5, .05), rep(0,9))
colnames(active_vs_index_single) <- c("Index", "Enrichment")

active_vs_index_single[,2] <- sapply(active_vs_index_single[,1], active_under_index, df = pdb_lig_sorted, cutoff = cutoff, col = 3)  

enrichment_vs_index_single <- active_vs_index_single
enrichment_vs_index_single[,2] <- (enrichment_vs_index_single[,2] / enrichment_vs_index_single[,1]) / (total_percent_active_1000_clusters)

accuracy_vs_index_single <- active_vs_index_single
colnames(accuracy_vs_index_single)[2] <- "Accuracy"
accuracy_vs_index_single[,2] <- active_vs_index_single[,2]/active_vs_index_single[,1]
plot(accuracy_vs_index_single)

## nicer plots:
accuracy_single_vs_multiple <- cbind(accuracy_vs_index[,c(1,2)], accuracy_vs_index_single[,2])
colnames(accuracy_single_vs_multiple) <- c("Index", "Multiple Scoring Functions", "Single Scoring Function")
make_plot(accuracy_single_vs_multiple, "Accuracy of Top 'n' Clusters Scored by Docking")

enrichment_single_vs_multiple <- accuracy_single_vs_multiple
enrichment_single_vs_multiple[,3] <- enrichment_single_vs_multiple[,3] / total_percent_active_1000_clusters
enrichment_single_vs_multiple[,2] <- enrichment_single_vs_multiple[,2] / total_percent_active_1000_clusters
make_plot(enrichment_single_vs_multiple, "Enrichment of Top 'n' Clusters Scored by Docking")
