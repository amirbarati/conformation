from msmbuilder.utils import verbosedump, verboseload
from io_functions import *
from functools import partial
from analysis import *
import multiprocessing as mp
import mdtraj as md
from msmbuilder.cluster import KMeans
from msmbuilder.cluster import MiniBatchKMeans

def cluster(data_dir, traj_dir, n_clusters, lag_time):
	clusterer_dir = "/scratch/users/enf/b2ar_analysis/clusterer_%d_t%d.h5" %(n_clusters, lag_time)
	if (os.path.exists(clusterer_dir)):
		print "Already clustered"
	else:
		reduced_data = verboseload(data_dir)
		trajs = np.concatenate(reduced_data)
		clusterer = MiniBatchKMedoids(n_clusters = n_clusters)
		clusterer.fit_transform(reduced_data)
		verbosedump(clusterer, "/scratch/users/enf/b2ar_analysis/clusterer_%d_t%d.h5" %(n_clusters, lag_time))	

def cluster_kmeans(tica_dir, data_dir, traj_dir, n_clusters, lag_time):
	clusterer_dir = "%s/clusterer_%dclusters.h5" %(tica_dir, n_clusters)
	if (os.path.exists(clusterer_dir)):
		print "Already clustered"
	else:
		print "Clustering by KMeans"
		reduced_data = verboseload(data_dir)
		trajs = np.concatenate(reduced_data)
		clusterer = KMeans(n_clusters = n_clusters, n_jobs = -1)
		clusterer.fit_transform(reduced_data)
		verbosedump(clusterer, clusterer_dir)	

def cluster_minikmeans(tica_dir, data_dir, traj_dir, n_clusters, lag_time):
	clusterer_dir = "%s/clusterer_%dclusters.h5" %(tica_dir, n_clusters)
	if (os.path.exists(clusterer_dir)):
		print "Already clustered"
	else:
		print "Clustering by KMeans"
		reduced_data = verboseload(data_dir)
		trajs = np.concatenate(reduced_data)
		clusterer = MiniBatchKMeans(n_clusters = n_clusters)
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

	for cluster in clusters_map.keys():
		print len(clusters_map[cluster])

	return clusters_map

def cos_to_means(clusterer_dir, features_dir):
	clusterer = verboseload(clusterer_dir)
	clusters_map = make_clusters_map(clusterer)

	features = verboseload(features_dir)
	feature_distances = {}

	for i in range(0, len(clusters_map.keys())):
		indices = clusters_map[i]
		k_mean = clusterer.cluster_centers_[i]
		print k_mean
		find_cos_partial = partial(find_cos, k_mean=k_mean, features = features)
		feature_distances_i = map(find_cos_partial, indices)
		feature_distances[i] = feature_distances_i

	print(feature_distances[0][0:10])
	sorted_map = {}

	print(feature_distances.keys())
	print(len(feature_distances.keys()))

	for i in range(0, len(feature_distances.keys())):
		sorted_features = sorted(feature_distances[i], key = lambda x: x[2], reverse = True)
		sorted_map[i] = sorted_features

	print sorted_map[0][0:10]
	return sorted_map

def find_dist(index, k_mean, features):
		traj = index[0]
		frame = index[1]
		conformation = features[traj][frame]
		a = conformation
		b = k_mean
		return (traj, frame, np.linalg.norm(b-a))

def dist_to_means(clusterer_dir, features_dir, n_samples = False, n_components = False, tica_coords_csv = False, kmeans_csv = False):
	clusterer = verboseload(clusterer_dir)
	clusters_map = make_clusters_map(clusterer)

	features = verboseload(features_dir)
	feature_distances = {}

	for i in range(0, len(clusters_map.keys())):
		indices = clusters_map[i]
		k_mean = clusterer.cluster_centers_[i]
		print k_mean
		find_dist_partial = partial(find_dist, k_mean=k_mean, features = features)
		feature_distances_i = map(find_dist_partial, indices)
		feature_distances[i] = feature_distances_i

	print(feature_distances[0][0:10])
	sorted_map = {}

	print(feature_distances.keys())
	print(len(feature_distances.keys()))

	for i in range(0, len(feature_distances.keys())):
		sorted_features = sorted(feature_distances[i], key = lambda x: x[2], reverse = False)
		sorted_map[i] = sorted_features

	if n_samples is not False and n_components is not False and tica_coords_csv is not False:
		tica_coords_map = {}
		for cluster_id in sorted_map.keys():
			for j in range(0, n_samples):
				sample = "cluster%d_sample%d" %(cluster_id, j)
				sample_tuple = sorted_map[cluster_id][j][0:2]
				sample_coords = features[sample_tuple[0]][sample_tuple[1]]
				tica_coords_map[sample] = sample_coords
		titles = ["sample"]
		for k in range(0, n_components):
			titles.append("component_%d" %k)
		print(tica_coords_map.keys()[0])
		print(tica_coords_map[tica_coords_map.keys()[0]])
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


	print sorted_map[0][0:10] 
	return sorted_map


def get_samples(cluster, trajectories, clusters_map, clusterer_dir, features_dir, traj_dir, save_dir, n_samples, method):
		for s in range(0, n_samples):
			if s == len(clusters_map[cluster]): return
			sample = clusters_map[cluster][s]
			traj_id = sample[0]
			frame = sample[1]
			traj = trajectories[traj_id]
			print traj

			traj_obj = md.load(traj)
			print traj_obj
			print frame

			top = md.load_frame(traj, index=frame).topology
			indices = [a.index for a in top.atoms if str(a.residue)[0:3] != "SOD" and str(a.residue)[0:3] != "CLA" and a.residue.resSeq < 341]
			#print indices

			conformation = md.load_frame(traj, index=frame, atom_indices=sorted(indices))
			conformation.save_pdb("%s/cluster%d_sample%d.pdb" %(save_dir, cluster, s))



def sample_clusters(clusterer_dir, features_dir, traj_dir, save_dir, n_samples, method):
	if method == "cos":	
		clusters_map = cos_to_means(clusterer_dir, features_dir)
	else:
		clusters_map = dist_to_means(clusterer_dir, features_dir)
	clusters = range(0, len(clusters_map.keys()))
	if not os.path.exists(save_dir): os.makedirs(save_dir)
	
	trajectories = get_trajectory_files(traj_dir, ".lh5")
	
	sampler = partial(get_samples, trajectories = trajectories, clusters_map = clusters_map, clusterer_dir = clusterer_dir, features_dir = features_dir, traj_dir = traj_dir, save_dir = save_dir, n_samples = n_samples, method = method)
	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(sampler, clusters)
	pool.terminate()
	#for cluster in clusters:
	#	if cluster != 118: continue
	#	sampler(cluster)
	#	print cluster

def get_pnas(cluster, clusters_map, pnas_active_distances, pnas_coords, n_samples):
		distances = []
		coords = []
		print "cluster = %d" %cluster
		for s in range(0, n_samples):
			print "sample  = %d" %s
			if s == len(clusters_map[cluster]): return
			sample = clusters_map[cluster][s]
			print sample
			traj_id = sample[0]
			frame = sample[1]
			print(pnas_active_distances[traj_id])
			print(np.shape(pnas_active_distances[traj_id]))
			active_pnas_distance = pnas_active_distances[traj_id][frame]
			pnas_coord = pnas_coords[traj_id][frame]
			distances.append(active_pnas_distance)
			coords.append(pnas_coord)
		return [distances, coords]

def cluster_pnas_distances(clusterer_dir, features_dir, active_pnas_dir, pnas_coords_dir, traj_dir, active_pnas_csv, pnas_coords_csv, n_samples, method):
	if method == "cos":	
		clusters_map = cos_to_means(clusterer_dir, features_dir)
	else:
		clusters_map = dist_to_means(clusterer_dir, features_dir)
	clusters = range(0, len(clusters_map.keys()))

	trajectories = get_trajectory_files(traj_dir, ".lh5")

	active_pnas_distances = verboseload(active_pnas_dir)
	pnas_coords = verboseload(pnas_coords_dir)

	sampler = partial(get_pnas, clusters_map = clusters_map, pnas_active_distances = active_pnas_distances, pnas_coords = pnas_coords, n_samples = n_samples)
	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pnas_feature = []
	for cluster in clusters: 
		pnas_feature.append(sampler(cluster))
	#pnas_feature = pool.map(sampler, clusters)
	pool.terminate()

	pnas_distance_map = {}
	pnas_coords_map = {}

	for i in range(0, len(clusters_map.keys())):
		try:
			pnas_distance = pnas_feature[i][0]
			print pnas_distance
			pnas_coord = pnas_feature[i][1]
			print pnas_coord
		except:
			continue
		for j in range(0, len(pnas_distance)):
			pnas_distance_map["cluster%d_sample%d" %(i, j)] = [pnas_distance[j]]
			pnas_coords_map["cluster%d_sample%d" %(i,j)] = pnas_coord[j]

	write_map_to_csv(active_pnas_csv, pnas_distance_map, ["sample", "active_pnas_distance"])
	write_map_to_csv(pnas_coords_csv, pnas_coords_map, ["sample", "tm6_tm3_dist", "npxxy_rmsd_inactive", "npxxy_rmsd_active", "connector_rmsd_inactive", "connector_rmsd_active"])





