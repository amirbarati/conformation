library(RColorBrewer)
library(classInt)
library(maptools)
plot.new()
tica.coords <- read.csv("/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/kmeans_csv.csv", row.names = 1)
tica.coords[1:10,]
append_cluster <- function(index) {
  return(paste("cluster", index, sep=""))
}
rownames(tica.coords) <- sapply(rownames(tica.coords), append_cluster)

active.pnas <- read.csv("/Users/evan/scratch_enf/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/active_pnas_joined.csv", row.names=1)
tica.pnas <- merge(tica.coords, active.pnas, by = "row.names", all = TRUE)
tica.pnas <- tica.pnas[-1]


plot_coords_and_pnas <- function(df, c1, c2) {
  pnas_col = dim(df)[2]
  colors <- rev(brewer.pal(3, "YlGnBu"))
  brks<-classIntervals(df[,pnas_col], n=3, style="fisher")
  brks <- brks$brks
  brks <- c(0.0, 0.3, 0.7, 1.3)
  #brks <- brks[order(-brks)]
  plot(df[,c(c1,c2)], col=colors[findInterval(df[,pnas_col], brks, all.inside=TRUE)], axes = T, pch=16)
  legend(x=1.25, y=-1, legend=leglabs(round(brks,3)), fill=colors, bty="n",x.intersp = .5, y.intersp = 1.0)
}

tica.pnas <- tica.pnas[order(tica.pnas[,dim(tica.pnas)[2]]),]
tica.pnas[1:10,]
plot_coords_and_pnas(tica.pnas, 1,2)
