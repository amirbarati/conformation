from PDB_Order_Fixer import PDB_Order_Fixer
import mdtraj as md
import os
import numpy as np
import h5py
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import KMeans
from msmbuilder.msm import MarkovStateModel
from sklearn.pipeline import Pipeline
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.utils import verbosedump, verboseload
from msmbuilder.dataset import dataset
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
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
from msmbuilder.cluster import MiniBatchKMeans
from matplotlib.backends.backend_pdf import PdfPages
from msmbuilder.msm import implied_timescales
import networkx as nx
import pytraj.io as mdio
from pytraj import adict
import random 
import subprocess
from subprocess import Popen
import sys
from grids import *
from io_functions import *
from custom_clusterer import *
from custom_tica import *
from custom_featurizer import *
from pdb_editing import *
from analysis import *
from io_functions import *
from grids import *
from topology_fixing import *
from subsampling import *
from conversions import *


def get_trajectory_files(traj_dir):
	traj_files = []
	for traj in os.listdir(traj_dir):
			if traj.endswith(".pdb"):
				traj_files.append("%s/%s" %(traj_dir,traj))
	return sorted(traj_files)

n_clusters = 500
lag_time = 10
msm_lag_time = 20
n_components = 5
n_samples = 20
#feature_types = "_switches_npxx_tm6_bp"
feature_types = ""

switch_residues = [130, 131, 208, 211, 219, 268, 272, 286, 288, 316, 322, 323, 326]
switch_npxx = [130, 131, 208, 211, 219, 268, 272, 286, 288, 316] + range(322,328)
tm6_residues = range(269, 299)
bp_residues = [82, 86, 93, 106, 110, 113, 114, 117, 118, 164, 191, 192, 193, 195, 199, 200, 203, 206, 208, 286, 289, 290, 293, 305, 308, 309, 312, 316]
dihedral_residues = list(set(switch_npxx + tm6_residues))
sampling_method = "dist"
precision = "SP"

sherlock_base = "/scratch/users/enf/b2ar_analysis"
biox3_base = "/home/enf/b2ar_analysis_sherlock_all"

if os.path.exists(sherlock_base):
	print("we are operating on sherlock")
	base = sherlock_base
elif os.path.exists(biox3_base):
	print("we are operating on biox3")
	base = biox3_base
else:
	print("WHERE ARE WE?")
	sys.exit()

traj_dir = "%s/subsampled_allprot_combined_reimaged" %base


tica_dir = "%s/tICA_t%d_n_components%d%s" %(base, lag_time, n_components, feature_types)
clusterer_dir = "%s/clusterer_%dclusters.h5" %(tica_dir, n_clusters)
#features_dir = "%s/features_allprot"
features_dir = "%s/features%s" %(base,feature_types)
features_known = "%s/features_known" %base
model_dir = tica_dir
projected_features_dir = "%s/phi_psi_chi2_allprot_projected.h5" %(tica_dir)
save_dir = "%s/clusters%d_n_components%d_n_samples%d_%s" %(tica_dir, n_clusters, n_components, n_samples, sampling_method)
#save_dir = "%s/reorder_test" %tica_dir
reimaged_dir = "%s/clusters%d_n_components%d_n_samples%d_%s_reimaged" %(tica_dir, n_clusters, n_components, n_samples, sampling_method)
active_ref_dir = "%s/3P0G_pymol_prepped.pdb" %base
inactive_ref_dir = "%s/2RH1_prepped.pdb" %base
scripts_dir = "%s/scripts" %base
script_dir = "%s/pymol_rmsd.py" %scripts_dir
pymol_fixpdb_dir = "%s/pymol_fixpdb.py" %scripts_dir
active_rmsd_dir =  "%s/active_rmsds.csv" %reimaged_dir
inactive_rmsd_dir = "%s/inactive_rmsd.csv" %reimaged_dir
active_pnas_dir = "%s/active_pnas_distances.csv" %reimaged_dir
inactive_pnas_dir = "%s/inactive_pnas_distances.csv" %reimaged_dir
analysis_file = "%s/rmsd_analyzed.csv" %reimaged_dir
combined_file = "%s/rmsd_combined.csv" %reimaged_dir
ligand_dir = "%s/ligprep_2/ligprep_2-out.maegz" %base
grid_dir = "%s/grids_n_clusters%d_n_samples%d_%s" %(tica_dir, n_clusters, n_samples, sampling_method)
docking_dir = "%s/docking_n_clusters%d_n_samples%d_%s_%s" %(tica_dir, n_clusters, n_samples, sampling_method, precision)
docking_summary = "%s/docking_summary.csv" %docking_dir
pnas_coords = "%s/pnas_coords.csv" %reimaged_dir
#mae_dir = "%s/docking_test" %tica_dir
mae_dir = reimaged_dir
grid_center = "64.4, 16.9, 11.99"

#featurize_custom(traj_dir, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = ["phi", "psi", "chi1", "chi2"], contact_residues = bp_residues)
#fit_and_transform(features_directory = features_dir, model_dir = tica_dir, stride=5, lag_time = lag_time, n_components = n_components)
#cluster_kmeans(tica_dir, projected_features_dir, traj_dir, n_clusters, lag_time)
#cluster_kmeans(tica_dir, projected_features_dir, traj_dir, n_clusters, lag_time)
#plot_tica_and_clusters(tica_dir, projected_features_dir, clusterer_dir, lag_time)
#sample_clusters(clusterer_dir, projected_features_dir, traj_dir, save_dir, n_samples, method = "dist")

#pymol_fixpdb(save_dir, pymol_fixpdb_dir)
#reorder(save_dir)
#reimage_trajs(save_dir, ext = ".pdb")
#remove_ter(reimaged_dir)
#reorder(reimaged_dir)
#pprep(mae_dir)
#rmsd_pymol(reimaged_dir, inactive_ref_dir, script_dir, inactive_rmsd_dir)
#rmsd_pymol(reimaged_dir, active_ref_dir, script_dir, active_rmsd_dir)
#analyze_docking_results(docking_dir)
#pnas_distances(reimaged_dir, inactive_ref_dir, active_ref_dir)
#analyze_rmsds(inactive_rmsd_dir, active_rmsd_dir, inactive_pnas_dir, active_pnas_dir, combined_file, analysis_file)
#plot_pnas_vs_docking(docking_summary, pnas_coords, "%s/pnas_vs_docking.pdf" %docking_dir)

#pnas_distance(traj_file, inactive_file, active_file)


#plot_all_tics_and_clusters(tica_dir, projected_features_dir, clusterer_dir, lag_time)
#pprep(mae_dir)
#generate_grids(mae_dir, grid_center, tica_dir, n_clusters, n_samples)
#dock_conformations(grid_dir, docking_dir, ligand_dir)
#analyze_docking_results(docking_dir)
#analyze_rmsds(inactive_rmsd_dir, active_rmsd_dir, inactive_pnas_dir, active_pnas_dir, combined_file, analysis_file)

#plot_pnas_vs_docking(docking_summary, pnas_coords, "%s/pnas_vs_docking.pdf" %docking_dir)

#reimage_trajs(save_dir)

#featurize_known(traj_dir, inactive_ref_dir)

#plot_hex("%s/features_known/A-00.h5")



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


