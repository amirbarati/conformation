from msmbuilder.utils import verbosedump, verboseload
from io_functions import *
from functools import partial
from analysis import *
import multiprocessing as mp
import mdtraj as md
#from msmbuilder.cluster import KMeans
#from msmbuilder.cluster import KCenters
from msmbuilder.cluster import MiniBatchKMeans
import random
import json
from sklearn import mixture
from msmbuilder.cluster import GMM

def cluster(data_dir, traj_dir, n_clusters, lag_time):
	clusterer_dir = "/scratch/users/enf/b2ar_analysis/clusterer_%d_t%d.h5" %(n_clusters, lag_time)
	if (os.path.exists(clusterer_dir)):
		print("Already clustered")
	else:
		try:
			reduced_data = verboseload(data_dir)
		except:
			reduced_data = load_dataset(data_dir)
		trajs = np.concatenate(reduced_data)
		clusterer = MiniBatchKMedoids(n_clusters = n_clusters)
		clusterer.fit_transform(reduced_data)
		verbosedump(clusterer, "/scratch/users/enf/b2ar_analysis/clusterer_%d_t%d.h5" %(n_clusters, lag_time))	

def recompute_cluster_means(means, tICs):
	return

def cluster_minikmeans(tica_dir, data_dir, traj_dir, n_clusters, clusterer_dir=None,tICs=None):
	if (os.path.exists(clusterer_dir)):
		print("Already clustered")
	else:
		print("Clustering by KMeans")
		try:
			reduced_data = verboseload(data_dir)
		except:
			reduced_data = load_dataset(data_dir)
		if tICs is not None:
			X = []
			for traj in reduced_data:
				X.append(traj[:,tICs])
		else:
			X = reduced_data

		clusterer = MiniBatchKMeans(n_clusters = n_clusters, n_init=10)
		clusterer.fit_transform(X)
		verbosedump(clusterer, clusterer_dir)

def cluster_gmm(projected_features_file, model_file, tICs=None, n_components=25):
	try:
		projected_features = verboseload(projected_features_file)
	except:
		projected_features = load_dataset(projected_features_file)
	X = np.concatenate(projected_features)

	X = []
	for traj in projected_features:
		X.append(traj[:,tICs])

	gmm = GMM(n_components=n_components, covariance_type='diag')
	print("Now fitting GMM model")
	gmm.fit(X)
	labels = gmm.predict(X)
	print("Completed GMM model. Saving now.")

	msmb_gmm = MSMB_GMM(labels, n_components, gmm.means_)
	verbosedump(msmb_gmm, model_file)

class MSMB_GMM(object):
	def __init__(self, labels, n_clusters, centers):
		self.labels_ = labels 
		self.n_clusters = n_clusters
		self.cluster_centers_ = centers

def cluster_kcenters(tica_dir, data_dir, traj_dir, n_clusters, lag_time):
	clusterer_dir = "%s/kcenters_clusterer_%dclusters.h5" %(tica_dir, n_clusters)
	if (os.path.exists(clusterer_dir)):
		print("Already clustered")
	else:
		print("Clustering by KMeans")
		reduced_data = verboseload(data_dir)
		trajs = np.concatenate(reduced_data)
		clusterer = KClusters(n_clusters = n_clusters)
		clusterer.fit_transform(reduced_data)
		verbosedump(clusterer, clusterer_dir)	

def make_clusters_map(clusterer):
	n_clusters = clusterer.n_clusters
	labels = clusterer.labels_
	clusters_map = {}

	for i in range(0,n_clusters):
		clusters_map[i] = set()

	for i in range(0, len(labels)):
		trajectory = labels[i]
		for j in range(0, len(trajectory)):
			cluster = trajectory[j]
			clusters_map[cluster].add((i,j))

	for cluster in list(clusters_map.keys()):
		print(len(clusters_map[cluster]))

	return clusters_map

def cos_to_means(clusterer_dir, features_dir):
	clusterer = verboseload(clusterer_dir)
	clusters_map = make_clusters_map(clusterer)

	features = verboseload(features_dir)
	feature_distances = {}

	for i in range(0, len(list(clusters_map.keys()))):
		indices = clusters_map[i]
		k_mean = clusterer.cluster_centers_[i]
		print(k_mean)
		find_cos_partial = partial(find_cos, k_mean=k_mean, features = features)
		feature_distances_i = list(map(find_cos_partial, indices))
		feature_distances[i] = feature_distances_i

	print((feature_distances[0][0:10]))
	sorted_map = {}

	print((list(feature_distances.keys())))
	print((len(list(feature_distances.keys()))))

	for i in range(0, len(list(feature_distances.keys()))):
		sorted_features = sorted(feature_distances[i], key = lambda x: x[2], reverse = True)
		sorted_map[i] = sorted_features

	print(sorted_map[0][0:10])
	return sorted_map

def find_dist(index, k_mean, features, tICs=None):
		traj = index[0]
		frame = index[1]
		conformation = features[traj][frame]
		a = conformation
		if tICs is not None:
			a = a[tICs]
		b = k_mean
		return (traj, frame, np.linalg.norm(b-a))

def dist_to_means(clusterer_dir, features_dir, n_samples = False, n_components = False, tica_coords_csv = False, kmeans_csv = False, tICs=None):
	clusterer = verboseload(clusterer_dir)
	clusters_map = make_clusters_map(clusterer)

	try: 
		features = verboseload(features_dir)
	except:
		features = load_dataset(features_dir)
	feature_distances = {}

	for i in range(0, len(list(clusters_map.keys()))):
		indices = clusters_map[i]
		k_mean = clusterer.cluster_centers_[i]
		print(k_mean)
		find_dist_partial = partial(find_dist, k_mean=k_mean, features = features, tICs=tICs)
		feature_distances_i = list(map(find_dist_partial, indices))
		feature_distances[i] = feature_distances_i

	print((feature_distances[0][0:10]))
	sorted_map = {}

	print((list(feature_distances.keys())))
	print((len(list(feature_distances.keys()))))

	for i in range(0, len(list(feature_distances.keys()))):
		sorted_features = sorted(feature_distances[i], key = lambda x: x[2], reverse = False)
		sorted_map[i] = sorted_features

	if n_samples is not False and n_components is not False and tica_coords_csv is not False:
		tica_coords_map = {}
		for cluster_id in list(sorted_map.keys()):
			for j in range(0, n_samples):
				sample = "cluster%d_sample%d" %(cluster_id, j)
				sample_tuple = sorted_map[cluster_id][j][0:2]
				sample_coords = features[sample_tuple[0]][sample_tuple[1]]
				tica_coords_map[sample] = sample_coords
		titles = ["sample"]
		for k in range(0, n_components):
			titles.append("component_%d" %k)
		print((list(tica_coords_map.keys())[0]))
		print((tica_coords_map[list(tica_coords_map.keys())[0]]))
		write_map_to_csv(tica_coords_csv, tica_coords_map, titles)

	if kmeans_csv is not False:
		kmeans_map = {}
		for cluster in range(0,clusterer.n_clusters):
			k_mean = clusterer.cluster_centers_[cluster]
			cluster_id = "cluster%d" %cluster
			kmeans_map[cluster_id] = k_mean
		titles = ["cluster"]
		for k in range(0, n_components):
			titles.append("component_%d" %k)
		write_map_to_csv(kmeans_csv, kmeans_map, titles)			


	print(sorted_map[0][0:10]) 
	return sorted_map


def get_samples(cluster, trajectories, clusters_map, clusterer_dir, features_dir, traj_dir, save_dir, n_samples, method, structure=None, residue_cutoff=10000):
	num_configurations = len(clusters_map[cluster])
	if method == "random":
		try:
			indices = random.sample(list(range(num_configurations)), n_samples)
		except:
			return(list(range(0, min(n_samples, num_configurations))))
		#print indices
	else:
		indices = list(range(0, min(n_samples, num_configurations)))
	
	for s in range(0, n_samples):
		if s == len(clusters_map[cluster]): return(indices[0:s])
		if method != "random":
			k = s
		else:
			k = indices[s]
		sample = clusters_map[cluster][k]
		traj_id = sample[0]
		frame = sample[1]
		traj = trajectories[traj_id]
		print(("cluster %d sample %d" %(cluster, k)))
		#print traj

		#traj_obj = md.load(traj)
		#print traj_obj
		#print frame

		if structure is None:
			top = md.load_frame(traj, index=frame).topology
		else:
			top = md.load_frame(traj, index=frame, top=structure).topology

		atom_indices = [a.index for a in top.atoms if str(a.residue)[0:3] != "SOD" and str(a.residue)[0:3] != "CLA" and a.residue.resSeq < residue_cutoff and str(a.residue)[0:3] != "POP" and not a.residue.is_water]
		#print indices
		if structure is None:
			conformation = md.load_frame(traj, index=frame, atom_indices=sorted(atom_indices))
		else:
			conformation = md.load_frame(traj, index=frame, atom_indices=sorted(atom_indices), top=structure)

		conformation.save_pdb("%s/cluster%d_sample%d.pdb" %(save_dir, cluster, s))
	
	print(cluster)
	#print(indices)
	#print(len(indices))
	return indices


def sample_clusters(clusterer_dir, features_dir, traj_dir, traj_ext, save_dir, n_samples, method, clusters_map_file = "", tICs=None, structure=None, worker_pool=None):
	if method == "cos":	
		clusters_map = cos_to_means(clusterer_dir, features_dir)
	elif method == "dist":
		clusters_map = dist_to_means(clusterer_dir, features_dir)
	elif method == "random":
		clusters_map = dist_to_means(clusterer_dir, features_dir, tICs=tICs)
	clusters = list(range(0, len(list(clusters_map.keys()))))
	if not os.path.exists(save_dir): os.makedirs(save_dir)
	
	trajectories = get_trajectory_files(traj_dir, traj_ext)
	
	sampler = partial(get_samples, trajectories = trajectories, clusters_map = clusters_map, clusterer_dir = clusterer_dir, features_dir = features_dir, traj_dir = traj_dir, save_dir = save_dir, n_samples = n_samples, method = method, structure=structure)
	if worker_pool is not None:
		list_of_indices = worker_pool.map_sync(sampler, clusters)
	else:
		num_workers = mp.cpu_count()
		pool = mp.Pool(num_workers)
		list_of_indices = pool.map(sampler, clusters)
		pool.terminate()
	print("Done sampling, now saving clusters map")
	for i in clusters:
		print(i)
		indices = list_of_indices[i]
		if method == "random":
			#print(len(indices))
			#print(indices[0:5])
			clusters_map[i] = [clusters_map[i][j] for j in indices]
	if method == "random":
		with open(clusters_map_file, 'w') as f:
			json.dump(clusters_map, f)


	
	#for cluster in clusters:
	#	if cluster != 118: continue
	#	sampler(cluster)
	#	print cluster





def get_pnas(cluster, clusters_map, pnas_coords, tica_coords, feature_coords, n_samples):
		distances = []
		coords = []
		tica_coords_list = []
		feature_coords_list = []
		print("cluster = %d" %cluster)
		for s in range(0, n_samples):
			print("sample  = %d" %s)
			if s == len(clusters_map[cluster]): return
			sample = clusters_map[cluster][s]
			print(sample)
			traj_id = sample[0]
			frame = sample[1]
			pnas_coord = pnas_coords[traj_id][frame]
			tica_coord = tica_coords[traj_id][frame]
			if feature_coords is not None: 
				feature_coord = feature_coords[traj_id][frame]
				feature_coords_list.append(feature_coord)
			coords.append(pnas_coord)
			tica_coords_list.append(tica_coord)
		return [coords, tica_coords_list, feature_coords_list]



def cluster_pnas_distances(clusterer_dir, features_dir, pnas_coords_dir, projected_features_dir, traj_dir, traj_ext, pnas_coords_csv, tica_coords_csv, feature_coords_csv, n_samples, method, coord_names, clusters_map_file = None):
	if method == "cos":	
		clusters_map = cos_to_means(clusterer_dir, features_dir)
	elif method == "random":
		with open(clusters_map_file) as f:
			clusters_map = json.load(f)
			clusters_map = {int(k):v for k,v in list(clusters_map.items())}
			print((list(clusters_map.keys())))
	elif method == "dist":
		clusters_map = dist_to_means(clusterer_dir, features_dir)
	else: 
		print("method not recognized")
		return
	clusters = list(range(0, len(list(clusters_map.keys()))))

	trajectories = get_trajectory_files(traj_dir, traj_ext)

	pnas_coords = verboseload(pnas_coords_dir)
	try:
		tica_coords = verboseload(projected_features_dir)
	except:
		tica_coords = load_dataset(projected_features_dir)
	feature_coords = None
	if features_dir is not None: feature_coords = load_file_list(None, features_dir, ".dataset")

	sampler = partial(get_pnas, clusters_map = clusters_map, pnas_coords = pnas_coords, tica_coords = tica_coords, feature_coords = feature_coords, n_samples = n_samples)
	num_workers = mp.cpu_count()
	#pool = mp.Pool(num_workers)
	pnas_feature = []
	for cluster in clusters: 
		pnas_feature.append(sampler(cluster))
	#pnas_feature = pool.map(sampler, clusters)
	#pool.terminate()

	pnas_distance_map = {}
	pnas_coords_map = {}
	tica_coords_map = {}
	feature_coords_map = {}

	for i in range(0, len(list(clusters_map.keys()))):
		try:
			pnas_coord = pnas_feature[i][0]
			print(pnas_coord)
			tica_coord = pnas_feature[i][1]
			if features_dir is not None: feature_coord = pnas_feature[i][2]
		except:
			continue
		for j in range(0, len(pnas_coord)):
			pnas_coords_map["cluster%d_sample%d" %(i,j)] = pnas_coord[j]
			tica_coords_map["cluster%d_sample%d" %(i,j)] = tica_coord[j]
			if features_dir is not None: feature_coords_map["cluster%d_sample%d" %(i,j)] = feature_coord[j]

	n_components = len(tica_coords_map[list(tica_coords_map.keys())[0]])

	write_map_to_csv(pnas_coords_csv, pnas_coords_map, ["sample"] + coord_names)
	tic_names = []
	for i in range(0, n_components):
		tic_names.append("tIC_%d" %i)
	write_map_to_csv(tica_coords_csv, tica_coords_map, ["sample"] + tic_names)
	if features_dir is not None: write_map_to_csv(feature_coords_csv, feature_coords_map, [])

def sample_features(clusterer_dir, features_dir, features_ext, n_samples, method, clusters_map_file = None):
	if method == "cos":	
		clusters_map = cos_to_means(clusterer_dir, features_dir)
	elif method == "random":
		with open(clusters_map_file) as f:
			clusters_map = json.load(f)
			clusters_map = {int(k):v for k,v in list(clusters_map.items())}
			print((list(clusters_map.keys())))
	elif method == "dist":
		clusters_map = dist_to_means(clusterer_dir, features_dir)
	else: 
		print("method not recognized")
		return
	clusters = list(range(0, len(list(clusters_map.keys()))))

	features = load_file_list(None, features_dir, features_ext)

	sampler = partial(get_pnas, clusters_map = clusters_map, pnas_coords = pnas_coords, tica_coords = tica_coords, n_samples = n_samples)
	num_workers = mp.cpu_count()
	#pool = mp.Pool(num_workers)
	pnas_feature = []
	for cluster in clusters: 
		pnas_feature.append(sampler(cluster))
	#pnas_feature = pool.map(sampler, clusters)
	#pool.terminate()

	pnas_distance_map = {}
	pnas_coords_map = {}
	tica_coords_map = {}

	for i in range(0, len(list(clusters_map.keys()))):
		try:
			pnas_distance = pnas_feature[i][0]
			print(pnas_distance)
			pnas_coord = pnas_feature[i][1]
			print(pnas_coord)
			tica_coord = pnas_feature[i][2]
		except:
			continue
		for j in range(0, len(pnas_distance)):
			pnas_distance_map["cluster%d_sample%d" %(i, j)] = [pnas_distance[j]]
			pnas_coords_map["cluster%d_sample%d" %(i,j)] = pnas_coord[j]
			tica_coords_map["cluster%d_sample%d" %(i,j)] = tica_coord[j]

	n_components = len(tica_coords_map[list(tica_coords_map.keys())[0]])

	write_map_to_csv(active_pnas_csv, pnas_distance_map, ["sample", "active_pnas_distance"])
	write_map_to_csv(pnas_coords_csv, pnas_coords_map, ["sample", "tm6_tm3_dist", "npxxy_rmsd_inactive", "npxxy_rmsd_active", "connector_rmsd_inactive", "connector_rmsd_active"])
	tic_names = []
	for i in range(0, n_components):
		tic_names.append("tIC_%d" %i)
	write_map_to_csv(tica_coords_csv, tica_coords_map, ["sample"] + tic_names)






