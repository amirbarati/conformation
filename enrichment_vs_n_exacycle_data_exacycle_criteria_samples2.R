library(ggplot2)
library(reshape2)

pnas_coords_all_csv <- "/Users/evan/scratch_enf/b2ar_analysis/exacycle_data/all_pnas_features/pnas_all_coords.csv"
pnas_coords_all <- data.frame(read.csv(pnas_coords_all_csv, stringsAsFactors = F, row.names=1))
badrows <- which(pnas_coords_all[,2] > 4.0 | pnas_coords_all[,3] > 4.0 | pnas_coords_all[,4] > 4.0 | pnas_coords_all[,5] > 4.0)
if(length(badrows) > 0) {
  pnas_coords_all <- pnas_coords_all[-badrows,]
}
colnames(pnas_coords_all) <- colnames(pnas_coords)
pnas_coords_all["tm6_tm3_dist"] <- 7.14 * pnas_coords_all["tm6_tm3_dist"]


all_pnas_rows <- find_active(pnas_coords_all)
all_active_rows <- all_pnas_rows[all_pnas_rows == T]


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

get_cluster_average <- function(cluster_name, df) {
  cluster_rows <- grep(paste(cluster_name, "_", sep=""), rownames(df))
  cluster_df <- df[cluster_rows,]
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


is_active <- function(row) {
  if(row["tm6_tm3_dist"] > 12.0 &  row["connector_rmsd_active"] < 1.0 & row["npxxy_rmsd_active"] < .8 & row["connector_rmsd_active"] < row["connector_rmsd_inactive"] & 0.8 < row["npxxy_rmsd_inactive"]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

find_active <- function(df) {
  active_rows <- apply(df, 1, is_active)
  return(active_rows)
}

#pnas_coords_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_dihedrals_switches_npxx_contact/analysis_n_clusters1000/pnas_coords_new.csv"
pnas_coords_csv <- "/Users/evan/scratch_enf/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"

pnas_coords <- data.frame(read.csv(pnas_coords_csv, stringsAsFactors = F, row.names=1))
pnas_coords["tm6_tm3_dist"] <- 7.14 * pnas_coords["tm6_tm3_dist"]
badrows <- which(pnas_coords[,2] > 4.0 | pnas_coords[,3] > 4.0 | pnas_coords[,4] > 4.0 | pnas_coords[,5] > 4.0)
if(length(badrows) > 0) {
  pnas_coords <- pnas_coords[-badrows,]
}
dim(pnas_coords)
pnas_coords[1:10,]


pnas_coords_averages <- cluster_averages(pnas_coords)

pnas_clusters_rows <- find_active(pnas_coords_averages)
active_clusters_rows <- pnas_clusters_rows[pnas_clusters_rows == T]
inactive_clusters_rows <- pnas_clusters_rows[pnas_clusters_rows == F]

pnas_rows <- find_active(pnas_coords)
active_rows <- pnas_rows[pnas_rows == T]
inactive_rows <- pnas_rows[pnas_rows == F]

print(length(active_rows))
print(length(active_rows)/dim(pnas_coords)[1])
print(length(all_active_rows)/(length(all_pnas_rows)))

print(length(active_clusters_rows))
print(length(active_clusters_rows)/dim(pnas_coords_averages)[1])
print(length(all_active_rows)/(length(all_pnas_rows)))
