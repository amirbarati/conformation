#Compare SASA distributions of exacycle vs anton data
exacycle.csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/exacycle_data/sasa_bp.csv"
exacycle <- data.frame(read.csv2(exacycle.csv, stringsAsFactors = F))
exacycle <- sapply(exacycle, as.numeric)

anton.csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/sasa_bp.csv"
anton <- read.csv2(anton.csv, stringsAsFactors = FALSE)
anton <- sapply(anton, as.numeric)

par(mfrow = c(2,1))
hist(anton[,1], breaks = seq(5, 25,1))
hist(exacycle[,1], breaks = seq(5, 25,1))

####compare docking

exacycle.csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_skip5_switches_pp_npxx_contact/docking_n_clusters3000_n_samples20_random_SP/3p0g_lig/docking_summary.csv"
exacycle <- read.csv2(exacycle.csv, stringsAsFactors = F, row.names=1, sep=",")
exacycle <- sapply(exacycle, as.numeric)

anton.csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/docking_n_clusters1000_n_samples10_dist_SP/3p0g_lig/docking_summary.csv"
anton <- read.csv2(anton.csv, stringsAsFactors = F, row.names=1, sep=",")
anton <- sapply(anton, as.numeric)

par(mfrow = c(2,1))
hist(anton, breaks = seq(0, 16, 0.2))
abline(v = 9.99, col="blue")
abline(v = 11.13, col = "red")
abline(v = 8.04, col = "red")
hist(exacycle, breaks = seq(0, 16, 0.3))
abline(v = 11.43, col = "blue")
abline(v = 11.13, col = "red")
abline(v = 8.04, col = "red")

exacycle.csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/exacycle_data/tICA_t10_n_components10_skip5_switches_pp_npxx_contact/docking_n_clusters3000_n_samples20_random_SP/s_carazolol/docking_summary.csv"
exacycle <- read.csv2(exacycle.csv, stringsAsFactors = F, row.names=1, sep=",")
exacycle <- sapply(exacycle, as.numeric)
#carazolol

anton.csv <- "/Users/Evan/downloads/vsp/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/docking_n_clusters1000_n_samples10_dist_SP/s-carazolol/docking_summary.csv"
anton <- read.csv2(anton.csv, stringsAsFactors = F, row.names=1, sep=",")
anton <- sapply(anton, as.numeric)

par(mfrow = c(2,1))
hist(anton, breaks = seq(0, 16, 0.2))
abline(v = 9.99, col="blue")
abline(v = 11.13, col = "red")
abline(v = 8.04, col = "red")
hist(exacycle, breaks = seq(0, 16, 0.3))
abline(v = 11.43, col = "blue")
abline(v = 11.13, col = "red")
abline(v = 8.04, col = "red")
