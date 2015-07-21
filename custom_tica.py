import os
from msmbuilder.decomposition import tICA
from io_functions import *
import multiprocessing as mp

def load_features(filename):
	return np.transpose(verboseload(filename))

def fit_and_transform(features_directory, model_dir, stride=5, lag_time=10, n_components = 5, tica_regularization = 0.05):
	if not os.path.exists(model_dir):
		os.makedirs(model_dir)

	projected_data_filename = "%s/phi_psi_chi2_allprot_projected.h5" %model_dir
	fit_model_filename  = "%s/phi_psi_chi2_allprot_tica_coords.h5" %model_dir
	#active_pdb_file = "/scratch/users/enf/b2ar_analysis/renamed_topologies/A-00.pdb"

	tica_model = tICA(n_components = n_components, lag_time = lag_time, gamma = tica_regularization)

	if not os.path.exists(projected_data_filename):
		print("loading feature files")
		feature_files = get_trajectory_files(features_directory, ext = ".h5")
		pool = mp.Pool(mp.cpu_count())
		features = pool.map(load_features, feature_files)
		pool.terminate()
		if np.shape(features[1]) != np.shape(features[1])[1]:
			for i in range(0, len(features)):
				features[i] = np.transpose(features[i])
		print np.shape(features[0])
		print np.shape(features[1])
		if not os.path.exists(fit_model_filename):
			print("fitting data to tICA model")
			fit_model = tica_model.fit(features)
			verbosedump(fit_model, fit_model_filename)
			transformed_data = fit_model.transform(features)
			verbosedump(transformed_data, projected_data_filename)
		else:
			print("loading tICA model")
			fit_model = verboseload(fit_model_filename)
			print("transforming")
			transformed_data = fit_model.transform(features)
			verbosedump(transformed_data, projected_data_filename)
	else:
		fit_model = verboseload(fit_model_filename)
		transformed_data = verboseload(projected_data_filename)

	print fit_model.summarize()

	#active_pdb = md.load(active_pdb_file)
	#top = active_pdb.topology
	#atom_indices = atom_indices = [a.index for a in top.atoms if str(a.residue)[0:3] != "POP" and not a.residue.is_water and str(a.residue)[0:2] != "NA" and str(a.residue)[0:2] != "CL"]
	#active_pdb = md.load(active_pdb_file, atom_indices=atom_indices)
	#read_and_featurize_custom(active_pdb_file, condition = "A-00_custom_features", location = "/scratch/users/enf/b2ar_analysis")
	#active_features = [np.transpose(verboseload("/scratch/users/enf/b2ar_analysis/A-00_custom_features.h5"))]
	#active_pdb_projected = fit_model.transform(active_features)
	#print(active_pdb_projected)
