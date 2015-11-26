from PDB_Order_Fixer import PDB_Order_Fixer
import mdtraj as md
import os
import numpy as np
import h5py
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import KMeans

from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.cluster import KMedoids
import datetime
import multiprocessing as mp
import glob
import copy
import gc
from functools import partial 
import itertools
import operator
from mdtraj.geometry import dihedral as ManualDihedral
import time
import fileinput
from msmbuilder.cluster import MiniBatchKMedoids

import random 
import subprocess
from subprocess import Popen
import sys
from io_functions import *
from custom_clusterer import *
from custom_tica import *
from custom_featurizer import *
from pdb_editing import *
from analysis import *
from io_functions import *
#from topology_fixing import *
from subsampling import *
from conversions import *
from custom_msm import *
from grids import *
from landmark_kernel_tica import *
#from pymol_color_tICA import *

import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
from rpy2.robjects import r
import rpy2.robjects.numpy2ri as numpy2ri
numpy2ri.activate()
from detect_intermediates import *

from feature_types import *



base = get_base()

R_functions = "%s/conformation/analysis.R" %base
R_analysis = "%s/conformation/b2ar_analysis.R" %base
ro.r.source(R_functions)
ro.r.source(R_analysis)
#r.assign("'R.scripts", R_scripts)
#ro.r('source(R.scripts)')

#traj_dir = "%s/reference_receptors" %base


sasa_file = "%s/sasa_bp.csv" %base
#r.assign('sasa.csv', sasa_file)
#mae_dir = "%s/docking_test" %tica_dir
grid_center = "64.4, 16.9, 11.99"

#inverse_ligands = get_ligands(inverse_agonist_dir)
#agonist_ligands = get_ligands(agonist_dir)

#inactive = md.load_frame(inactive_ref_dir, top = inactive_ref_dir, index = 0)
#active = md.load_frame(active_ref_dir, top = active_ref_dir, index = 0)

#compute_pnas_distance(traj_dir, inactive, active)

'''
precision = "SP"
tica_dir = "/scratch/users/enf/b2ar_analysis/reference_docking"
docking_dir = "%s/docking_%s" %(tica_dir, precision)
grid_dir = "%s/reference_grids" %tica_dir
ligands_dir = agonist_dir
docking_multiple_ligands = "%s/all_docking_combined.csv" %docking_dir
dock_ligands_and_receptors(grid_dir, docking_dir, ligands_dir, precision = precision, ext = "-out.maegz", chosen_ligands = False, chosen_receptors = False, parallel = True)
analyze_docking_results_multiple(docking_dir, precision = precision)
compute_aggregate_docking_scores(docking_multiple_ligands, docking_dir)
'''
'''
precision = "SP"
tica_dir = "/scratch/users/enf/b2ar_analysis/reference_docking"
docking_dir = "%s/docking_%s" %(tica_dir, precision)
grid_dir = "%s/reference_grids" %tica_dir
ligands_dir = inverse_agonist_dir
docking_multiple_ligands = "%s/all_docking_combined.csv" %docking_dir
dock_ligands_and_receptors(grid_dir, docking_dir, ligands_dir, precision = precision, ext = "-out.maegz", chosen_ligands = False, chosen_receptors = False, parallel = True)
analyze_docking_results_multiple(docking_dir, precision = precision, summary = docking_multiple_ligands)
#compute_aggregate_docking_scores(docking_multiple_ligands, docking_dir)
'''

'''
top_dock = "%s/top50_docking.txt" %tica_dir
to_dock = []
dockfile = open(top_dock, "rb")
for line in dockfile.readlines():
	line = line.split('\n')[0]
	to_dock.append(line)
print to_dock
'''


def make_extreme_tIC_barplots(tica_extremes_dir, feature_residues_csv, n_components):
	feature_files = [f for f in get_trajectory_files(tica_extremes_dir, ext=".csv") if "standardized" in f]
	feature_files = [f for f in feature_files if "standardized" in f]
	for i in range(1, n_components+1):
		low_file = [f for f in feature_files if "tIC.%d_" %i in f and "low" in f][0]
		print low_file
		high_file = [f for f in feature_files if "tIC.%d_" %i in f and "high" in f][0]
		r['analyze.extreme.tic.values'](low_file, high_file, feature_residues_csv, i, tica_extremes_dir)
	return

#reimage_traj_new("%s/A-00.h5" %traj_dir, base, "", ".h5")
#featurize_pnas_distance(base, base, "-00.h5", inactive_ref_dir, active_ref_dir, "%s/pnas_inactive_dist_test.csv" %base, "%s/pnas_active_dist_test.csv" %base, "%s/pnas_coords_dir.csv" %base, None, "%s/pnas_active_dist_test.csv" %base, "%s/pnas_all_dist_test.csv" %base, scale = 7.14, residues_map = None)


from tica_variables import *

####Featurize with PNAS distances and coords, 2D####
#featurize_pnas_distance(traj_dir, pnas_features_dir, ".h5", inactive_ref_dir, active_ref_dir, inactive_pnas_distances_dir, active_pnas_distances_dir, pnas_coords_dir, scale = 1.0)
####
#featurize_pnas_distance_pdbs(reimaged_dir, "%s/combined.h5" %reimaged_dir, features_dir, inactive_ref_dir, active_ref_dir, inactive_pnas_distances_dir, active_pnas_distances_dir, pnas_coords_dir, scale = 7.14)

#featurize_pnas_distance(traj_dir, whole_trajectory_pnas, ".h5", inactive_ref_dir, active_ref_dir, inactive_pnas_distances_dir, active_pnas_distances_dir, pnas_coords_dir, None, active_pnas_all_distances_dir, pnas_all_coords_csv, scale = 7.14, residues_map = None)
#featurize_pnas_distance(ref_receptors_dir, pnas_features_dir, ".pdb", inactive_ref_dir, active_ref_dir, inactive_pnas_distances_dir, active_pnas_distances_dir, pnas_coords_dir, None, active_pnas_all_distances_dir, pnas_all_coords_csv, scale = 7.14)
#featurize_sasa(traj_dir = reimaged_dir, traj_ext = ".pdb", bp_residues = bp_residues, sasa_file = sasa_file, residues_map = None, anton = False, skip = 1, stride = 1)
#plot_columns(whole_trajectory_pnas, pnas_coords_dir, titles = pnas_titles, tICA = False, scale = 7.14, refcoords_file = "%s/ref_coords.h5" %ref_receptors_dir)

#plot_hex(pnas_coords_dir, pnas_coords_hexbin_dir)
#plot_col(pnas_coords_dir, pnas_coords_active_colors_dir, active_pnas_distances_dir)
#to_dock = ["cluster0_sample1", "cluster0_sample2", "cluster0_sample3"]

residues_map = generate_residues_map(residues_map_csv)
new_residues_map = {}
for k,v in residues_map.iteritems():
    new_residues_map[(0, k)] = (0, v)
residues_map = new_residues_map
contact_residues = [res for res in contact_residues if res in residues_map.keys()]
 
#Featurize reference receptors
#featurize_contacts_custom(traj_dir, features_dir = features_dir, traj_ext = traj_ext, contact_residue_pairs_csv = feature_residues_csv, featurized_traj_and_frame = None, dihedral_residues =  [], dihedral_types = ["phi", "psi", "chi1", "chi2"], contact_residues =  contact_residues, residues_map = residues_map, contact_cutoff = cutoff, parallel = parallel, exacycle = exacycle)
#save_features_to_residues_map(get_trajectory_files(traj_dir, ext = traj_ext)[0], contact_residues, feature_residues_csv, cutoff=cutoff,residues_map = residues_map, exacycle=True)
#featurize_contacts_custom(ref_receptors_dir, features_dir = ref_features_dir, traj_ext = ".pdb", featurized_traj_and_frame = [get_trajectory_files(traj_dir, traj_ext)[0], 0], contact_residue_pairs_csv = feature_residues_csv, dihedral_residues =  [], dihedral_types = ["phi", "psi", "chi1", "chi2"], contact_residues =  contact_residues, residues_map = residues_map, contact_cutoff = cutoff, exacycle = False)
#if exacycle:
#    while not os.path.exists("%s/trj39998.dataset" %features_dir): 
#        print("Waiting for featurization to complete")
#        time.sleep(60*10)
#standardize_features(features_dir, ".h5", standardized_features_dir)

'''
fit_and_transform(features_directory = features_dir, model_dir = tica_dir, stride=5, lag_time = lag_time, n_components = n_components, sparse = sparse, wolf = wolf, rho = rho, shrinkage = shrinkage, parallel=feature_parallel, traj_ext = traj_ext)

plot_pnas_vs_tics(pnas_coords_dir, projected_features_dir, ["tm6_tm3_dist", "rmsd_npxxy_inactive", "rmsd_npxxy_active", "rmsd_connector_inactive", "rmsd_connector_active"], tica_dir)

transform(existing_model = projection_operator_dir, features_directory = ref_features_dir, tica_dir = ref_tica_dir)
#plot_columns(tica_dir, projected_features_dir, titles = None, tICA = True, scale = 1.0, refcoords_file = ref_tica_coords)
#plot_columns_3d(tica_dir, projected_features_dir, titles = None, tICA = True, scale = 1.0, refcoords_file = ref_tica_coords, columns = [1,6, 8])
#plot_columns_3d_contour(tica_dir, projected_features_dir, titles = None, tICA = True, scale = 1.0, refcoords_file = ref_tica_coords, columns = [2,7,9])
#tica_coords = np.concatenate(load_file(projected_features_dir))
#print(np.shape(tica_coords))
#r['tIC.mixture.models'](tica_coords, tica_dir)
#sample_tIC_extremes(projected_features_dir, features_dir, standardized_features_dir, tica_extremes_dir, ".h5", percentile)
#make_extreme_tIC_barplots(tica_extremes_dir, feature_residues_csv, n_components)
#plot_all_tics(tica_dir, projected_features_dir, lag_time)

cluster_minikmeans(tica_dir, projected_features_dir, traj_dir, n_clusters, lag_time)
#cluster_kmeans(tica_dir, projected_features_dir, traj_dir, n_clusters, lag_time)
#find_missing_features(traj_dir, features_dir)

sample_clusters(clusterer_dir, projected_features_dir, traj_dir, traj_ext, save_dir, n_samples, method = sampling_method, clusters_map_file = clusters_map_file)
#dist_to_means(clusterer_dir, projected_features_dir, n_samples = n_samples, n_components = n_components, tica_coords_csv = tica_coords_csv, kmeans_csv = kmeans_csv)
#reverse_sign_csv(docking_joined)
#plot_all_tics_samples(kmeans_csv, analysis_dir, docking_csv = docking_joined, specific_clusters = [49, 353, 994, 397, 456, 517, 51])

cluster_pnas_distances(clusterer_dir, features_dir, active_pnas_distances_dir, pnas_coords_dir, projected_features_dir, traj_dir, traj_ext, active_pnas_distances_new_csv, pnas_coords_csv, tica_coords_csv, features_csv, n_samples, sampling_method, clusters_map_file = clusters_map_file)

find_most_important_residues_in_tIC(get_trajectory_files(traj_dir, traj_ext)[0], projection_operator_dir, feature_residues_csv, contact_residues, tic_residue_csv, feature_coefs_csv, duplicated_feature_coefs_csv, cutoff = cutoff)
#r['compute.decision.trees']("", tica_coords_csv, tica_classes_csv, ro.IntVector([2,3,4,7,9,11,24]), "", features_csv, feature_residues_csv, analysis_dir)
#r['compute.rf']("", tica_coords_csv, tica_classes_csv, ro.IntVector([2,3,4,7,9,11,24]), "", features_csv, feature_residues_csv, analysis_dir)
#for filename in [f for f in get_trajectory_files(analysis_dir, ".rda") if "rf" in f] : r['write.rf.gini'](filename, "%s.csv" %filename.split(".")[0])
#interpret_tIC(ref_receptors_dir, inactive_ref_dir, "%s/tIC7rf_importance.csv" %analysis_dir, analysis_dir, 7)
#r['sample.tic.ranges'](tica_coords_csv, tica_classes_csv, ro.IntVector([2,7,9]), tica_samples_csv)

r['do.analysis'](tica_dir, analysis_dir, pnas_coords_csv, tica_coords_csv, features_dir, docking_multiple_ligands)


with open(active_clusters_csv, 'rb') as f:
    reader = csv.reader(f)
    active_clusters = list(reader)[0]
active_clusters = [int(c[7:]) for c in active_clusters]
with open(intermediate_clusters_csv, 'rb') as f:
    reader = csv.reader(f)
    intermediate_clusters = list(reader)[0]
intermediate_clusters = [int(c[7:]) for c in intermediate_clusters]
print(intermediate_clusters[0:10])
with open(inactive_clusters_csv, 'rb') as f:
    reader = csv.reader(f)
    inactive_clusters = list(reader)[0]
inactive_clusters = [int(c[7:]) for c in inactive_clusters]

plot_all_tics_and_clusters(tica_dir, projected_features_dir, clusterer_dir, lag_time, label = "dot", active_cluster_ids = active_clusters, intermediate_cluster_ids = intermediate_clusters, inactive_cluster_ids = inactive_clusters)


find_correlation(features_dir, projected_features_dir, mutual_information_csv, pearson_csv, bins=50, exacycle = exacycle)
r['analyze.tic.feature.correlations'](pearson_csv, feature_residues_csv, tica_dir, "pearson_tic_feature_coefficients", " :: Feature Pearson Correlation Coefficients")
#r['analyze.tic.feature.correlations'](mutual_information_csv, feature_residues_csv, tica_dir, "MI_tic_feature_coefficients", " :: Feature Mutual Information Score")

#check_tica_vs_features(projected_features_dir, features_dir)

#plot_timescales(clusterer_dir, n_clusters, tica_dir)
#build_msm(clusterer_dir, msm_lag_time)
#macrostate_pcca(msm_model_dir, clusterer_dir, n_macrostates, macrostate_dir)
#construct_graph(msm_model_dir, clusterer_dir, n_clusters, lag_time, msm_lag_time, graph_file, inactive = None,active = active_pnas_joined, pnas_clusters_averages = pnas_clusters_averages, tica_clusters_averages = tica_clusters_averages, docking = None, macrostate = macrostate_dir) #docking=aggregate_docking_joined)
'''

compute_gmms_R(projected_features_dir, max_components=10, save_dir=gmm_dir, max_j=n_components)
#plot_tics_gmm(gmm_dir, projected_features_dir, gmm_dir, R=True, titles = None, tICA = False, scale = 1.0, refcoords_file = ref_tica_coords)

#compute_one_vs_all_rf_models(features_dir, projected_features_dir, gmm_dir, rf_dir, n_trees = 10, n_tica_components=25)
compute_overall_rf_models(features_dir, projected_features_dir, gmm_dir, rf_dir, n_trees = 500, n_tica_components=10)


#LANDMARK Kernel tICA


from ktica_variables import *
'''
residues_map = generate_residues_map(residues_map_csv)
contact_residues = [res for res in contact_residues if res in residues_map.keys()]

if ktica_ticaTraj:
    landmark_ktica_ticaTraj(projected_features_dir, ori_clusterer_dir, tica_dir, clusters_map_file = clusters_map_file, landmarks_dir = landmarks_dir, nystroem_components=1000, n_components = n_components, lag_time=5, nystroem_data_filename = nystroem_data_filename, fit_model_filename = ktica_fit_model_filename, projected_data_filename = ktica_projected_data_filename, landmark_subsample=landmark_subsample, sparse = sparse, wolf = wolf, rho = rho, shrinkage = shrinkage)
else:
    landmark_ktica(features_dir, combined_features_dir, tica_dir, clusters_map_file = clusters_map_file, landmarks_dir = landmarks_dir, nystroem_components=1000, n_components=k_tica_components, lag_time=5, nystroem_data_filename = nystroem_data_filename, fit_model_filename = ktica_fit_model_filename, projected_data_filename = ktica_projected_data_filename, landmark_subsample=landmark_subsample, sparse = sparse, wolf = wolf, rho = rho, shrinkage = shrinkage)
plot_pnas_vs_tics(pnas_coords_dir, ktica_projected_data_filename, ["tm6_tm3_dist", "rmsd_npxxy_inactive", "rmsd_npxxy_active", "rmsd_connector_inactive", "rmsd_connector_active"], tica_dir)

#ktica_test(features_dir, tica_dir, landmark_indices = None, nystroem_components=1000, tica_components=10, lag_time=5, nystroem_data_filename = nystroem_data_filename, fit_model_filename = ktica_fit_model_filename, projected_data_filename = ktica_projected_data_filename)
#plot_columns(tica_dir, projected_features_dir, titles = None, tICA = True, scale = 1.0, refcoords_file = ref_tica_coords)
##plot_all_tics(tica_dir, ktica_projected_data_filename, lag_time)
cluster_minikmeans(tica_dir, ktica_projected_data_filename, traj_dir, n_clusters, lag_time)

sample_clusters(clusterer_dir, ktica_projected_data_filename, traj_dir, traj_ext, save_dir, n_samples, method = sampling_method, clusters_map_file = ktica_clusters_map_file)

cluster_pnas_distances(clusterer_dir, None, active_pnas_distances_dir, pnas_coords_dir, ktica_projected_data_filename, traj_dir, traj_ext, active_pnas_distances_new_csv, pnas_coords_csv, tica_coords_csv, features_csv, n_samples, sampling_method, clusters_map_file = ktica_clusters_map_file)

r['do.analysis'](tica_dir, analysis_dir, pnas_coords_csv, tica_coords_csv, features_dir)

with open(active_clusters_csv, 'rb') as f:
    reader = csv.reader(f)
    active_clusters = list(reader)[0]
active_clusters = [int(c[7:]) for c in active_clusters]
with open(intermediate_clusters_csv, 'rb') as f:
    reader = csv.reader(f)
    intermediate_clusters = list(reader)[0]
intermediate_clusters = [int(c[7:]) for c in intermediate_clusters]
print(intermediate_clusters[0:10])
with open(inactive_clusters_csv, 'rb') as f:
    reader = csv.reader(f)
    inactive_clusters = list(reader)[0]
inactive_clusters = [int(c[7:]) for c in inactive_clusters]

plot_all_tics_and_clusters(tica_dir, ktica_projected_data_filename, clusterer_dir, lag_time, label = "dot", active_cluster_ids = active_clusters, intermediate_cluster_ids = intermediate_clusters, inactive_cluster_ids = inactive_clusters)

find_correlation(features_dir, ktica_projected_data_filename, mutual_information_csv, pearson_csv, bins=50, exacycle = exacycle)
r['analyze.tic.feature.correlations'](pearson_csv, feature_residues_csv, tica_dir, "pearson_tic_feature_coefficients", " :: Feature Pearson Correlation Coefficients")
#r['analyze.tic.feature.correlations'](mutual_information_csv, feature_residues_csv, tica_dir, "MI_tic_feature_coefficients", " :: Feature Mutual Information Score")
'''

#compute_gmms(projected_features_dir, gmm_dir, gmm_max_components)
#plot_tics_gmm(gmm_dir, projected_features_dir, gmm_dir, titles = None, tICA = False, scale = 1.0, refcoords_file = ref_tica_coords)
#compute_one_vs_all_rf_models(features_dir, projected_features_dir, gmm_dir, rf_dir, n_trees = 10, n_tica_components=25)

#compute_overall_rf_models(features_dir, projected_features_dir, gmm_dir, rf_dir, n_trees = n_trees, n_tica_components=25)


#while(not os.path.exists(pearson_csv)):
#	print("Waiting")
#	time.sleep(100)

#time.sleep(60*130)

#featurize_pnas_distance("%s/reference_receptors" %base, "%s/reference_receptors" %base, ".pdb", inactive_ref_dir, active_ref_dir, "%s/reference_receptors/inactive_pnas_distances_ref.h5" %base, "%s/reference_receptors/active_pnas_distances_ref.h5" %base, "%s/reference_receptors/ref_coords.h5" %base, "%s/reference_receptors/inactive_distances.csv" %base, "%s/reference_receptors/active_distances.csv" %base, "%s/reference_receptors/ref_coords.csv" %base, scale = 1.0)
#convert_matrix_to_map("%s/reference_receptors/ref_coords.h5" %base, "%s/reference_receptors" %base, ".pdb", "%s/reference_receptors/ref_coords.csv" %base)

#featurize_pnas_distance(reimaged_dir, features_dir, ".pdb", inactive_ref_dir, active_ref_dir, inactive_pnas_distances_dir, active_pnas_distances_dir, pnas_coords_dir, scale = 7.14)
#convert_matrix_to_map(active_pnas_distances_dir, reimaged_dir, ".pdb",active_pnas_distances_new_csv)

#pymol_fixpdb(save_dir, pymol_fixpdb_dir)
#reorder(save_dir)
#reimage_trajs(save_dir, ext = ".pdb")
reimaged_dir = save_dir
mae_dir = reimaged_dir
#remove_ter(reimaged_dir)
#reorder(reimaged_dir)
#rmsd_pymol(reimaged_dir, inactive_ref_dir, script_dir, inactive_rmsd_dir)
#rmsd_pymol(reimaged_dir, active_ref_dir, script_dir, active_rmsd_dir)
#analyze_docking_results(docking_dir)
#pnas_distances(reimaged_dir, inactive_ref_dir, active_ref_dir)
#analyze_rmsds(inactive_rmsd_dir, active_rmsd_dir, inactive_pnas_dir, active_pnas_dir, combined_file, analysis_file)
#plot_pnas_vs_docking(docking_summary, pnas_coords, "%s/pnas_vs_docking.pdf" %docking_dir)

#pnas_distance(traj_file, inactive_file, active_file)


#plot_all_tics_and_clusters(tica_dir, projected_features_dir, clusterer_dir, lag_time)
indices = [800,1000]
chosen_receptors = []
for i in range(indices[0], indices[1]):
	for j in range(0, n_samples):
		chosen_receptors.append("cluster%d_sample%d" %(i, j))
#pprep(mae_dir, ref = active_ref_dir, chosen_receptors = chosen_receptors)
#generate_grids(mae_dir, grid_center, grid_dir, remove_lig = "BIA", chosen_receptors = chosen_receptors)
#dock_conformations(grid_dir = grid_dir, docking_dir = docking_dir, ligand_dir = ligand_dir, chosen_jobs = False, precision = precision)
#analyze_docking_results(docking_dir, "BI", "SP", docking_summary)
#combine_csv_list([docking_summary, active_pnas_dir], docking_distances_file)

#docking_joined_map = convert_csv_to_joined_map(docking_summary, docking_joined)[0]
#pnas_joined_map = convert_csv_to_joined_map(active_pnas_distances_new_csv, active_pnas_joined)[0]
#inactive_pnas_joined_map = convert_csv_to_joined_map(inactive_pnas_dir, )[0]

#docking_averages = calc_mean(docking_joined_map)
#pnas_averages = calc_mean(pnas_joined_map)

#write_map_to_csv(docking_joined, docking_averages, ["cluster", "mean_docking_score"])
#write_map_to_csv(active_pnas_means, pnas_averages, ["cluster", "pnas_averages"])

#top_n =  top_n_scoring_samples(active_pnas_means, score_type = "pnas_averages", n = 10, n_samples = 1)
#print top_n
#plot_all_tics_and_clusters(tica_dir, projected_features_dir, clusterer_dir, lag_time, cluster_ids = [692, 349, 311, 356, 705, 866, 527, 763, 0, 132, 799, 685, 754])


#combine_csv_list([docking_joined, active_pnas_joined], docking_pnas_joined)

#top_n =  top_n_scoring_samples(aggregate_docking_joined, score_type = "mean_aggregate_docking_z_score", n = n_mmgbsa, n_samples = n_samples)
#print top_n



#rmsd_pymol(reimaged_dir, inactive_ref_dir, script_dir, inactive_rmsd_dir)
#rmsd_pymol(reimaged_dir, active_ref_dir, script_dir, active_rmsd_dir)
#active_pnas_map = convert_csv_to_map(active_pnas_dir)
#active_pnas_stats = calc_mean_and_stdev(active_pnas_map)
#write_map_to_csv("%s/active_pnas_stats.csv" %reimaged_dir, active_pnas_stats, ["cluster, average_pnas, stdev_pnas"])
#analyze_rmsds(inactive_rmsd_dir, active_rmsd_dir, inactive_pnas_dir, active_pnas_dir, combined_file, analysis_file)
#combine_docking_distances(docking_summary, combined_file, docking_dir)

#pnas_values = convert_matrix_to_map(active_pnas_distances_dir, reimaged_dir, ".pdb", active_pnas_distances_new_csv)
#compute_means_ligands(docking_dir, active_pnas_joined, inverse_ligands + agonist_ligands)
#convert_matrix_list_to_list("%s/features_pnas/active_pnas_distances.h5" %base, "%s/features_pnas/all_pnas_distances.csv" %base)
#unzip_receptors(grid_dir, chosen_receptors)
inverse_ligands = get_ligands(inverse_agonist_dir)
agonist_ligands = get_ligands(agonist_dir)
#prepare_ligands(inverse_agonist_dir, ext = ".sdf")

#dock_ligands_and_receptors(grid_dir, docking_dir, agonist_dir, precision = precision, ext = "-out.maegz", chosen_ligands = agonist_ligands, chosen_receptors = chosen_receptors, parallel = "receptor", grid_ext = ".grd")
#dock_ligands_and_receptors(grid_dir, docking_dir, inverse_agonist_dir, precision = precision, ext = "-out.maegz", chosen_ligands = inverse_ligands, chosen_receptors = chosen_receptors, parallel = "receptor", grid_ext = ".grd")

#analyze_docking_results_multiple(docking_dir, precision = "SP", ligands = (inverse_ligands + agonist_ligands), summary = docking_multiple_ligands, redo = True)
#compute_aggregate_scores(docking_multiple_ligands, inverse_agonists = inverse_ligands, summary = aggregate_docking, z_scores_csv = docking_z_scores_csv)
#combine_csv_list([docking_joined, active_pnas_joined], aggregate_docking_pnas_joined)

#aggregate_docking_joined_map = convert_csv_to_joined_map(aggregate_docking, aggregate_docking_joined)[0]
#aggregate_docking_means = calc_mean(aggregate_docking_joined_map)
#write_map_to_csv(aggregate_docking_joined, aggregate_docking_means, ["cluster", "mean_aggregate_docking_z_score"])
#combine_csv_list([aggregate_docking_joined, active_pnas_joined], aggregate_docking_pnas_joined)

#mmgbsa(docking_dir, mmgbsa_dir, chosen_jobs = top_n)
#mmgbsa_ligands_and_receptors(docking_dir, mmgbsa_dir, inverse_ligands, chosen_receptors = top_n)
#analyze_mmgbsa_results_multiple(mmgbsa_dir, summary = mmgbsa_csv , ligands = inverse_ligands + agonist_ligands, chosen_receptors = top_n)
#compute_aggregate_scores(mmgbsa_csv, inverse_agonists = inverse_ligands, summary = aggregate_mmgbsa, z_scores_csv = mmgbsa_z_scores_csv)
#aggregate_mmgbsa_joined_map = convert_csv_to_joined_map(aggregate_mmgbsa, aggregate_mmgbsa_joined)[0]
#aggregate_mmgbsa_means = calc_mean(aggregate_mmgbsa_joined_map)
#write_map_to_csv(aggregate_mmgbsa_joined, aggregate_mmgbsa_means, ["cluster", "mean_aggregate_mmgbsa_z_score"])
#combine_csv_list([aggregate_mmgbsa_joined, active_pnas_joined], aggregate_mmgbsa_pnas_joined)

#analyze_mmgbsa_results(mmgbsa_dir, mmgbsa_csv)
#combine_docking_mmgbsa(docking_distances_file, mmgbsa_csv, mmgbsa_dir, mmgbsa_docking_distances)
#plot_pnas_vs_docking(docking_summary, pnas_coords, docking_pdf, selected = to_dock)

#plot_pnas_vs_docking(docking_summary, pnas_coords, "%s/pnas_vs_docking.pdf" %docking_dir)

#reimage_trajs(save_dir)

#featurize_known(traj_dir, inactive_ref_dir)

#plot_hex("%s/features_known/A-00.h5")


#### Dock SP, Find top n_s and do Dock XP on those, Find top n_x and MM_GBSA on those #####
#generate_grids(mae_dir, grid_center, tica_dir, n_clusters, n_samples)
#precision = "SP"
#to_dock = False
#dock_conformations(grid_dir, tica_dir, precision = precision, to_dock = to_dock)
#analyze_docking_results_multiple()

####


'''important residues for GPCRs:
Asp3.49	-- Asp130
Arg3.50	--	Arg131
Phe5.47	--	Phe208
Pro5.50	--	Pro211
Tyr5.58	--	Tyr219
Glu6.30	--	Glu268
Thr6.34	--	Thr272
Trp6.48	--	Trp286
Pro6.50	--	Pro288
Lys7.43	--	Lys316
Asn7.49	--	Asn322
Pro7.50	--	Pro323
Tyr7.53	--	Tyr326	
'''



'''
B2AR binding pocket residues:

Met82
Val86
His93
Cyx106
Trp109
Thr110
Asp113
Val114
Val117
Thr118
Thr164
Cyx191
Asp192
Phe193
Thr195
Tyr199
Ala200
Ser203
Val206
Ser208
Trp286
Phe289
Phe290
Asn293
Lys305
Tyr308
Ile309
Asn312
Tyr316



'''

'''
http://pubs.acs.org/doi/pdf/10.1021/jm800710x
Full aognist list:
procaterol
R-isoproterenol
R-epinephrine
TA-2005
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2693451/
http://scicurve.com/paper/9249239

inverse agonist list:
carazolol 
atenolol
timolol
alprenolol

http://www.ncbi.nlm.nih.gov/pubmed/7908406

'''



'''
def extract_centers(clusterer_dir, traj_dir, lag_time):
	clusterer = verboseload(clusterer_dir)
	n_clusters = clusterer.n_clusters
	labels = clusterer.labels_
	sample_conformations = []
	visited_clusters = set()
	all_visited = False

	for i in range(0, len(labels)):
		trajectory = labels[i]
		for j in range(0, len(trajectory)):
			label = trajectory[j]
			if label not in visited_clusters:
				sample_conformations.append((label,i,j))
				visited_clusters.add(label)
				if len(visited_clusters) == n_clusters:
					print("sampled all clusters")
					all_visited = True
					break
		if all_visited == True: break

	trajectories = get_trajectory_files(traj_dir)

	for cluster in sample_conformations:
		cluster_id = cluster[0]
		traj = trajectories[cluster[1]]
		frame = cluster[2]

		conformation = md.load_frame(traj,index=frame)

		save_dir = "%s/%d_clusters_t%d" %(n_clusters, lag_time)
		if not os.path.exists(save_dir): os.makedirs(save_dir)
		conformation.save_pdb("%s/%d.pdb" %(save_dir, cluster_id))
'''

'''
For each subdirectory i in DESRES:
	For each subdir j in i: 
		get list of all dcd files
		split subdirectory to find the correct directory to save the file
		apply function read_and_save_traj on each dcd file
'''

'''

traj_dir = "/home/harrigan/data/gpcr/DESRES/DESRES-Trajectory_pnas2011b-H-05-all/pnas2011b-H-05-all"
traj_files = get_trajectory_files(traj_dir)[0:5]

a = datetime.datetime.now().replace(microsecond=0)
trajectories = map(read_trajectory, traj_files)
#trajectories = []
#for traj_file in traj_files:
#	trajectories.append(read_trajectory(traj_file))
b = datetime.datetime.now().replace(microsecond=0)
print(b-a)
'''


'''
dataset = []
trajs = []

traj_dir = "/home/harrigan/data/gpcr/DESRES/DESRES-Trajectory_pnas2011b-H-05-all/pnas2011b-H-05-all"
traj_files = []

if not (os.path.isfile("%s/H-05/%s" %("combined_traj_stride10.h5"))):
	print("traj not loaded yet")
	for traj in os.listdir(traj_dir):
		if traj.endswith(".dcd"):
			traj_files.append("%s/%s" %(traj_dir,traj))
	traj_files.sort()
	traj = md.load(traj_files, top = "/home/harrigan/compute/wetmsm/gpcr/des/system_mae_to_pdb/des_trajs/DESRES-Trajectory_pnas2011b-H-05-all/system.pdb", stride=10)
	traj = traj[0].join(traj[1:])
	traj.save("%s/H-05/%s" %("combined_traj_stride10.h5"))
else:
	print("loading h5 traj")
	traj = md.load("%s/H-05/%s" %("combined_traj_stride10.h5"))

'''
'''
if not (os.path.isfile("phi_psi_chi2_features_vd_stride10.h5")):
	print("featurizing")
	phi_psi_chi2 = DihedralFeaturizer(types=['phi','psi','chi2'])
	features = phi_psi_chi2.transform(traj_list = traj)
	print("finished featurizing")
	verbosedump(features, "phi_psi_chi2_features_vd_stride10.h5")
else:
	print("loading existing features")
	features = verboseload("phi_psi_chi2_features_vd_stride10.h5")

if not (os.path.isfile("reduced_phi_psi_chi_stride10.h5")):
	print("Fitting tICA model")
	tica_model = tICA(n_components=4)
	fitted_model = tica_model.fit(features)
	reduced_data = fitted_model.transform(features)
	verbosedump(reduced_data, "reduced_phi_psi_chi_stride10.h5")
	print(tica_model.summarize())
else:
	reduced_data = verboseload("reduced_phi_psi_chi_stride10.h5")

clusterer = KMedoids(n_clusters=9)

clusters = clusterer.fit_transform(reduced_data)[0]

center_locations = []

for i in range(0, len(clusters)):
	print i
	for j in range(0, len(clusterer.cluster_centers_)):
		if np.linalg.norm(reduced_data[0][i] - clusterer.cluster_centers_[j]) < 0.001:
			print("found match")
			center_locations.append(i)

print(center_locations)

for center in center_locations:
	frame = md.load_frame("combined_traj_stride10.h5", index=center)
	frame.save_pdb("frame_%d.pdb" %(center))






#trajs = np.concatenate(reduced_data)
#plt.hexbin(trajs[:,0], trajs[:,1], bins='log', mincnt=1)
#plt.show()
'''


