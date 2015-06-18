library(RColorBrewer)
library(classInt)
library(maptools)

pnas.coords <- read.csv(file = "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords.csv", sep = ",", row.names = 1)
pnas.coords[1:10,]
pnas.coords.clusters <- cbind(rep(0,1000), rep(0,1000))
colnames(pnas.coords.clusters) <- colnames(pnas.coords)

sample.zero.rows <- which(grepl("sample0", rownames(pnas.coords)))
pnas.coords.clusters <- pnas.coords[sample.zero.rows,]

remove_sample <- function(my_string) {
  return(strsplit(my_string, "_")[[1]][1])
}

rownames(pnas.coords.clusters) <- sapply(rownames(pnas.coords.clusters), remove_sample)

average.docking.clusters <- read.csv(file = "/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking_joined.csv", row.names=1)
pnas.coords.docking <- merge(pnas.coords.clusters, average.docking.clusters, by = "row.names", all = TRUE)
dfrows <- pnas.coords.docking[,"Row.names"]
pnas.coords.docking <- pnas.coords.docking[,-1]
rownames(pnas.coords.docking) <- dfrows
pnas.coords.docking.sorted <- pnas.coords.docking[order(-pnas.coords.docking[,3]),]
plot(pnas.coords.docking[,1], pnas.coords.docking[,2])

plot_coords_and_docking <- function(df, ori_df, refcoords) {
  colors <- brewer.pal(9, "YlGnBu")
  brks<-classIntervals(ori_df[,3], n=9, style="fisher")
  brks <- brks$brks
  plot(df[,c(1,2)], col=colors[findInterval(df[,3], brks, all.inside=TRUE)], axes = T, pch=16)
  par(new = T)
  points(refcoords, col = "red", pch=18)
  par(new = F)
}

refcoords <- read.csv(file = "/Users/evan/scratch_enf/b2ar_analysis/reference_receptors/ref_coords.csv", row.names = 1)
plot_coords_and_docking(pnas.coords.docking.sorted[1:50,], pnas.coords.docking.sorted, refcoords)

#Make ROC Curves 

pnas.coords.docking.sorted[1:10,]
