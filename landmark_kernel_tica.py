'''
load clusterer
load tica object 
load featurized trajectories 
load clusters map for random or re-generate it if dist/cos

landmarks: use i,j coords to get a list of vectors where each vector is the features in feature space of the landmarks

for each feature vector, compute its kernel distance to each of the landmarks
	save each feature object in a new feature fold

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
		landmark_vectors.append(feature)

	landmark_vectors = np.concatenate(landmark_vectors)
	





