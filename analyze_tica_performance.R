library(Hmisc)

#pnas.coords.csv <- "/Users/evan/vsp/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/pnas_coords_new.csv"
#tica.coords.csv <- "/Users/evan/vsp/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/tica_coords.csv"

#pnas.coords.csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/pnas_coords_new.csv"
#tica.coords.csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t20_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/tica_coords.csv"

pnas.coords.csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/pnas_coords_new.csv"
tica.coords.csv <- "/Users/evan/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact/analysis_n_clusters1000/tica_coords.csv"


pnas.coords <- data.frame(read.csv(pnas.coords.csv, stringsAsFactors = F, row.names=1))[,c(1,2,3,4,5)]
tica.coords <- data.frame(read.csv(tica.coords.csv, stringsAsFactors = F, row.names=1))

badrows <- which(pnas_coords_all[,2] > 4.0 | pnas_coords_all[,3] > 4.0 | pnas_coords_all[,4] > 4.0 | pnas_coords_all[,5] > 4.0)
if(length(badrows) > 0) {
  pnas.coords <- pnas.coords[-badrows,]
  tica.coords <- tica.coords[-badrows,]
}

pnas.coords["tm6_tm3_dist"] <- 7.14 * pnas.coords["tm6_tm3_dist"]

rcorr(cbind(pnas.coords[,5], tica.coords[,2]),type=c("pearson"))
plot(pnas.coords[,3], tica.coords[,9])

#for(i in (1:20)) {
#  print(rcorr(cbind(pnas.coords[,3], tica.coords[,i]),type=c("pearson")))
#}
