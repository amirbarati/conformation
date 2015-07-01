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
  if(row["tm6_tm3_dist"] > 12.0 & row["connector_rmsd_active"] < 1.0 & row["npxxy_rmsd_active"] < .7 & row["connector_rmsd_active"] < row["connector_rmsd_inactive"] & row["npxxy_rmsd_inactive"] > 0.8) {
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
#pnas_coords_csv <- "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"
#pnas_coords_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components5_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/pnas_coords_new.csv"
#pnas_coords_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"
pnas.coords.csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/pnas_coords_new.csv"


#tica_coords_csv <- "/Users/evan/vsp/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/tica_coords.csv"

pnas_coords <- data.frame(read.csv(pnas_coords_csv, stringsAsFactors = F, row.names=1))[,c(1,2,3,4,5)]
pnas_coords["tm6_tm3_dist"] <- 7.14 * pnas_coords["tm6_tm3_dist"]
tica_coords <- data.frame(read.csv(tica_coords_csv, stringsAsFactors = F, row.names=1))
badrows <- which(pnas_coords[,2] > 4.0 | pnas_coords[,3] > 4.0 | pnas_coords[,4] > 4.0 | pnas_coords[,5] > 4.0)
if(length(badrows) > 0) {
  pnas_coords <- pnas_coords[-badrows,]
  tica_coords <- tica_coords[-badrows,]
}
dim(pnas_coords)
pnas_coords[1:10,]

docking_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/all_docking_combined.csv"
docking <- data.frame(read.csv(docking_csv, stringsAsFactors = F, row.names=1))
docking.means <- apply(docking, 2, mean)
docking.sds <- apply(docking, 2, sd)
docking.z <- data.frame(matrix(data = 0, nrow = dim(docking)[1], dim(docking)[2]))  
for(i in (1:dim(docking)[2])) {
  docking.z[,i] <- (docking[,i] - docking.means[i])/docking.sds[i]
}
rownames(docking.z) <- rownames(docking)
colnames(docking.z) <- colnames(docking)

pnas_coords_averages <- cluster_averages(pnas_coords)
tica_coords_averages <- cluster_averages(tica_coords)
docking_averages <- cluster_averages(docking.z)

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

docking_csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking.csv"
docking.python <- data.frame(read.csv(docking_csv, stringsAsFactors = F, row.names=1))

docking.agonists <- docking.z[,c(1, 3, 4, 5, 7)]
docking.means <- apply(docking.agonists, 1, mean)
docking.means.sorted <- docking.means[order(-1.0*docking.means)]

samples.active.docking <- merge(docking.means, pnas_rows, by = "row.names", all = TRUE)
rownames(samples.active.docking) <- samples.active.docking[,1] 
samples.active.docking <- samples.active.docking[,2:dim(samples.active.docking)[2]]
colnames(samples.active.docking) <-c("aggregate_docking_score", "is_active")

samples.active.docking.ordered <- samples.active.docking[order(-1.0 * samples.active.docking$aggregate_docking_score),]

