from msmbuilder.decomposition import tICA, SparseTICA
from sklearn.kernel_approximation import Nystroem as ScikitNystroem
from LandmarkNystroem import Nystroem as LandmarkNystroem
from msmbuilder.decomposition.base import MultiSequenceDecompositionMixin
from msmbuilder.utils import verbosedump, verboseload
import numpy as np
import multiprocessing as mp
from io_functions import *
import json
import inspect
from scipy.linalg import svd
import types
from sklearn.metrics.pairwise import pairwise_kernels
import tables

'''
def landmark_fit(self, X, y=None):
	n_samples = X.shape[0]
	basis = y
	print("Here is landmark fit data:")
	print(np.shape(X))
	print(np.shape(basis))

	basis_kernel = pairwise_kernels(basis, metric=self.kernel, filter_params=True,**self._get_kernel_params())

	# sqrt of kernel matrix on basis vectors
	U, S, V = svd(basis_kernel)
	S = np.maximum(S, 1e-12)
	self.normalization_ = np.dot(U * 1. / np.sqrt(S), V)
	return self

	print(inspect.getargspec(ScikitNystroem.fit))
'''

class Nystroem(MultiSequenceDecompositionMixin, LandmarkNystroem):
	pass





'''
load clusterer
load tica object 
load featurized trajectories 
load clusters map for random or re-generate it if dist/cos

landmarks: use i,j coords to get a list of vectors where each vector is the features in feature space of the landmarks

for each feature vector, compute its kernel distance to each of the landmarks
	save each feature object in a new feature fold

'''
'''
def landmark_kernel_tica(clusterer_dir, cluster_map_file, n_landmarks_per_cluster, sampling_method, features_dir, reimaged_dir, kernel_features_dir):
	if method == "cos":	
		clusters_map = cos_to_means(clusterer_dir, features_dir)
	elif method == "random":
		with open(clusters_map_file) as f:
			clusters_map = json.load(f)
			clusters_map = {int(k):v for k,v in clusters_map.items()}
			print(clusters_map.keys())
	elif method == "dist":
		clusters_map = dist_to_means(clusterer_dir, features_dir)
	else: 
		print("method not recognized")
		return

	landmark_frames = []
	landmark_vectors = []
	for i in range(0,len(clusters_map.keys())):
		frames = clusters_map[i]
		for j in range(0, n_landmarks_per_cluster):
			if j == len(frames): continue

			landmark_frames.append(frames[j])

	feature_files = get_trajectory_files(features_dir, ext = ".h5")
	for frame in landmark_frames:
		features = verboseload(feature_files[frame[0]])
		feature = features[frame[1]]
		landmark_vectors.append(feature)g

	landmark_vectors = np.concatenate(landmark_vectors)
'''



def landmark_ktica(features_dir, combined_features_dir, tica_dir, clusters_map_file = "", landmarks_dir = "", nystroem_components=1000, n_components=10, lag_time=5, nystroem_data_filename = "", fit_model_filename = "", projected_data_filename = "", landmark_subsample=1, sparse = False, shrinkage = 0.05, wolf = True, rho = 0.01):
	with open(clusters_map_file) as f:
		clusters_map = json.load(f)
		clusters_map = {int(k):v for k,v in clusters_map.items()}
	i = 0
	for cluster_id,sample_list in clusters_map.items():
		for sample in sample_list:
			i+= 1

	print i

	if not sparse:
		if shrinkage is None:
			tica_model = tICA(n_components = n_components, lag_time = lag_time)
		else:
			tica_model = tICA(n_components = n_components, lag_time = lag_time, shrinkage = shrinkage)
		
	else:
		if shrinkage is None:
			tica_model = SparseTICA(n_components = n_components, lag_time = lag_time, rho = rho)
		else:
			tica_model = SparseTICA(n_components = n_components, lag_time = lag_time, rho = rho, shrinkage = shrinkage)

	if not os.path.exists(nystroem_data_filename):
		features = load_file_list(get_trajectory_files(features_dir, ext = ".dataset"))
		if os.path.exists(landmarks_dir):
			landmarks = verboseload(landmarks_dir)
			print(np.shape(landmarks))
		else:
			landmarks = []
			for cluster_id,sample_list in clusters_map.items():
				for sample in sample_list:
					traj = sample[0]
					frame = sample[1]
					landmark = features[traj][frame]
					landmarks.append(landmark)
			verbosedump(landmarks, landmarks_dir)

		#if os.path.exists(nystroem_data_filename):
		#	nyx = verboseload(nystroem_data_filename)
		#else:
			#features = verboseload(combined_features_dir)
		print("here's what goes into the combined class:")
		landmarks = [landmarks[i] for i in range(0,np.shape(landmarks)[0]) if i%landmark_subsample==0] #%landmark_subsample == 0]
		#print(np.shape(features))
		print(np.shape(landmarks))
		print(type(landmarks))
		#print(landmarks[9930:np.shape(landmarks)[0]])
		#features = [f[0:100,0:100] for f in features]
		nys = Nystroem(n_components = np.shape(landmarks)[0], basis = landmarks)#np.shape(landmarks)[0])# basis=landmarks)
		nyx = nys.fit_transform(features)
		del features
		del landmarks
		try:
			save_dataset(nyx, nystroem_data_filename)
		except:
			os.system("rm -rf %s" %nystroem_data_filename)
			save_dataset(nyx, nystroem_data_filename)
		#f = tables.openFile(nystroem_data_filename, 'w')
		#atom = tables.Atom.from_dtype(nyx.dtype)
		#ds = f.createCArray(f.root, 'somename', atom, nyx.shape)
		#ds[:] = nyx
		#f.close()
		#np.savez_compressed(nystroem_data_filename, nyx)
	else:
		nyx = load_dataset(nystroem_data_filename)

	print(np.shape(nyx))
	print(dir(nyx))

	if not os.path.exists(projected_data_filename):
		fit_model = tica_model.fit(nyx)
		verbosedump(fit_model, fit_model_filename)
		transformed_data = fit_model.transform(nyx)
		del(nyx)
		try:
			save_dataset(transformed_data, projected_data_filename)
		except:
			os.system("rm -rf %s" %projected_data_filename)
			save_dataset(transformed_data, projected_data_filename)
	else:
		print("Already performed landmark kernel tICA.")

def landmark_ktica_ticaTraj(tica_dir, clusterer_dir, ktica_dir, clusters_map_file = "", landmarks_dir = "", nystroem_components=1000, n_components=10, lag_time=5, nystroem_data_filename = "", fit_model_filename = "", projected_data_filename = "", landmark_subsample=1, sparse = False, wolf = True, rho = 0.01, shrinkage = None):
	if not os.path.exists(ktica_dir): os.makedirs(ktica_dir)
	
	if not sparse:
		if shrinkage is None:
			tica_model = tICA(n_components = n_components, lag_time = lag_time)
		else:
			tica_model = tICA(n_components = n_components, lag_time = lag_time, shrinkage = shrinkage)
		
	else:
		if shrinkage is None:
			tica_model = SparseTICA(n_components = n_components, lag_time = lag_time, rho = rho)
		else:
			tica_model = SparseTICA(n_components = n_components, lag_time = lag_time, rho = rho, shrinkage = shrinkage)

	if not os.path.exists(nystroem_data_filename):
		clusterer = verboseload(clusterer_dir)
		tica = verboseload(tica_dir)
		features = tica
		clusters = clusterer.cluster_centers_
		landmarks = clusters

		print("here's what goes into the combined class:")
		#print(np.shape(features))
		print(np.shape(landmarks))
		print(type(landmarks))
		nys = Nystroem(n_components = np.shape(landmarks)[0], basis = landmarks)#np.shape(landmarks)[0])# basis=landmarks)
		nyx = nys.fit_transform(features)
		del features
		del landmarks
		try:
			save_dataset(nyx, nystroem_data_filename)
		except:
			os.system("rm -rf %s" %nystroem_data_filename)
			save_dataset(nyx, nystroem_data_filename)
	else:
		nyx = load_dataset(nystroem_data_filename)

	print(np.shape(nyx))
	print(dir(nyx))

	if not os.path.exists(projected_data_filename):
		fit_model = tica_model.fit(nyx)
		verbosedump(fit_model, fit_model_filename)
		transformed_data = fit_model.transform(nyx)
		del(nyx)
		try:
			save_dataset(transformed_data, projected_data_filename)
		except:
			os.system("rm -rf %s" %projected_data_filename)
			save_dataset(transformed_data, projected_data_filename)
	else:
		print("Already performed landmark kernel tICA.")
	

def ktica_test(features_dir, tica_dir, landmark_indices = None, nystroem_components=1000, tica_components=10, lag_time=5, nystroem_data_filename = "", fit_model_filename = "", projected_data_filename = ""):
	nys = Nystroem(n_components=nystroem_components)
	tica_model = tICA(n_components = tica_components, lag_time = lag_time)
	feature_files = get_trajectory_files(features_dir, ext = ".h5")[0:3]

	#if os.path.exists(nystroem_data_filename):
	#	nyx = verboseload(nystroem_data_filename)
	#else:
	features = load_file_list(feature_files)
	nyx = nys.fit_transform(features)
	verbosedump(nyx, nystroem_data_filename)

	print(np.shape(nyx))
	print(dir(nyx))

	fit_model = tica_model.fit(nyx)
	verbosedump(fit_model, fit_model_filename)
	transformed_data = fit_model.transform(nyx)
	verbosedump(transformed_data, projected_data_filename)

	return
	
	#features = 





