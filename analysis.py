
def calc_mean_and_stdev(rmsd_map):
	stats_map = {}
	for key in rmsd_map.keys():
		rmsds = np.array(rmsd_map[key])
		mean = np.mean(rmsds, axis = 0)
		stdev = np.std(rmsds, axis = 0)
		stats_map[key] = (mean, stdev)
	return stats_map

def find_cos(index, k_mean, features):
		traj = index[0]
		frame = index[1]
		conformation = features[traj][frame]
		a = conformation
		b = k_mean
		return (traj, frame, np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)))

def rmsd_npxxy(traj, inactive):
	npxxy_atoms = [a.index for a in traj.topology.atoms if a.residue.resSeq in range(322,328) and a.is_backbone]
	traj_stripped = traj.atom_slice(npxxy_atoms)

	npxxy_atoms_target = [a.index for a in inactive.topology.atoms if a.residue.resSeq in range(322,328) and a.is_backbone]
	inactive_stripped = inactive.atom_slice(npxxy_atoms_target)

	traj_stripped_aligned = traj_stripped.superpose(inactive_stripped)
	rmsds = md.rmsd(traj_stripped, inactive_stripped) * 10.0
	return rmsds

def helix6_helix3_dist(traj):
	atom_3 = [a.index for a in traj.topology.atoms if a.residue.resSeq == 131 and a.name == "CA"][0]
	atom_6 = [a.index for a in traj.topology.atoms if a.residue.resSeq == 272 and a.name == "CA"][0]
	indices = np.empty([1,2])
	indices[0] = [atom_3, atom_6]
	dist = md.compute_distances(traj, indices) * 10.0
	return dist



def plot_pnas_vs_docking(docking_dir, pnas_dir, save_dir):
	dock_scores = convert_csv_to_map_nocombine(docking_dir)
	pnas_vectors = convert_csv_to_map_nocombine(pnas_dir)
	x = []
	y = []
	c = []
	for key in pnas_vectors.keys():
		#print "PNAS"
		#print key
		#print pnas_vectors[key]
		if key in dock_scores.keys():
			#print pnas_vectors[key]
			x.append(pnas_vectors[key][0])
			y.append(pnas_vectors[key][1])
			c.append(abs(dock_scores[key][0]))
	print dock_scores.keys()
	plt.scatter(x, y, c=c, s=50, cmap = mpl.cm.RdYlBu, norm = LogNorm())
	pp = PdfPages(save_dir)
	pp.savefig()
	pp.close()


def pnas_distance(traj_file, inactive_file, active_file):
	traj = md.load(traj_file)
	inactive = md.load(inactive_file)
	active = md.load(active_file)
	scale = 7.14
	inactive_tuple = np.array([helix6_helix3_dist(inactive) / scale, rmsd_npxxy(inactive, inactive)])
	active_tuple = np.array([helix6_helix3_dist(active) / scale, rmsd_npxxy(active, inactive)])
	traj_coords = [helix6_helix3_dist(traj) / scale, rmsd_npxxy(traj, inactive)]
	traj_tuple = np.array(traj_coords)
	active_dist = np.linalg.norm(traj_tuple - active_tuple)
	inactive_dist = np.linalg.norm(traj_tuple - inactive_tuple)
	distances = [inactive_dist, active_dist]
	print distances[1]
	return [traj_coords, distances]
	

def pnas_distances(traj_dir, inactive_file, active_file):
	scale = 7.14
	save_file_i = "%s/inactive_pnas_distances.csv" %traj_dir
	save_file_a = "%s/active_pnas_distances.csv" %traj_dir
	coord_file = "%s/pnas_coords.csv" %traj_dir
	distances_i = open(save_file_i, "wb")
	distances_i.write("conformation, inactive_pnas_distance, a \n")
	distances_a = open(save_file_a, "wb")
	distances_a.write("conformation, active_pnas_distance, a \n")
	coordinates = open(coord_file, "wb")
	coordinates.write("cluster, tm6_tm3_dist, npxxy_rmsd_inactive \n")

	pnas_partial = partial(pnas_distance, inactive_file = inactive_file, active_file = active_file)

	trajs = get_trajectory_files(traj_dir, ext = ".pdb")
	pool = mp.Pool(mp.cpu_count())
	distances = pool.map(pnas_partial, trajs)
	pool.terminate()

	for i in range(0, len(trajs)):
		traj = trajs[i]
		pnas_coord = distances[i][0]
		pnas_dist = distances[i][1]
		traj_name = traj.split("/")[len(traj.split("/"))-1].split(".")[0]
		distances_i.write("%s, %f \n" %(traj_name, pnas_dist[0]))
		distances_a.write("%s, %f \n" %(traj_name, pnas_dist[1]))
		coordinates.write("%s, %f, %f \n" %(traj_name, (pnas_coord[0]*scale), pnas_coord[1]))

	distances_i.close()
	distances_a.close()
	coordinates.close()
	return [save_file_i, save_file_a]



def plot_hex(transformed_data_file):
	transformed_data = verboseload(transformed_data_file)
	trajs = transformed_data
	#trajs = np.concatenate(transformed_data)
	plt.hexbin(trajs[:,1], trajs[:,0], bins='log', mincnt=1)
	pp = PdfPages("%s.pdf" %transformed_data_file)
	pp.savefig()
	pp.close()
	return


def save_pdb(traj_dir, clusterer, i):
	location = clusterer.cluster_ids_[i,:]
	traj = get_trajectory_files(traj_dir)[location[0]]
	print("traj = %s, frame = %d" %(traj, location[1]))
	conformation = md.load_frame(traj, location[1])
	conformation.save_pdb("/scratch/users/enf/b2ar_analysis/clusters_1000_allprot/%d.pdb" %i)
	return None

def get_cluster_centers(clusterer_dir, traj_dir):
	clusterer = verboseload(clusterer_dir)

	centers = clusterer.cluster_centers_

	save_pdb_partial = partial(save_pdb, traj_dir = traj_dir, clusterer = clusterer)
	indices = range(0, np.shape(centers)[0])
	pool = mp.Pool(mp.cpu_count())
	pool.map(save_pdb_partial, indices)
	pool.terminate()
	return

def plot_tica(transformed_data_dir, lag_time):
	transformed_data = verboseload(transformed_data_dir)
	trajs = np.concatenate(transformed_data)
	plt.hexbin(trajs[:,0], trajs[:,1], bins='log', mincnt=1)
	pp = PdfPages("/scratch/users/enf/b2ar_analysis/tica_phi_psi_chi2_t%d.pdf" %lag_time)
	pp.savefig()
	pp.close()


def plot_tica_and_clusters(tica_dir, transformed_data_dir, clusterer_dir, lag_time, component_i = 0, component_j = 1, cluster_ids = False):
	transformed_data = verboseload(transformed_data_dir)
	clusterer = verboseload(clusterer_dir)

	trajs = np.concatenate(transformed_data)
	plt.hexbin(trajs[:,component_i], trajs[:,component_j], bins='log', mincnt=1)

	centers = clusterer.cluster_centers_
	if cluster_ids is False:
		for i in range(0, np.shape(centers)[0]):
			center = centers[i,:]
			plt.annotate('%d' %i, xy=(center[component_i],center[component_j]), xytext=(center[component_i], center[component_j]),size=6)
	else:
		for i in cluster_ids:
			center = centers[i,:]
			plt.annotate('%d' %i, xy=(center[component_i],center[component_j]), xytext=(center[component_i], center[component_j]),size=6)

	pp = PdfPages("%s/c%d_c%d_clusters%d.pdf" %(tica_dir, component_i, component_j, np.shape(centers)[0]))
	pp.savefig()
	pp.close()

def plot_all_tics_and_clusters(tica_dir, transformed_data_dir, clusterer_dir, lag_time, cluster_ids = False):
	transformed_data = verboseload(transformed_data_dir)
	num_tics = np.shape(transformed_data[0])[1]
	print "Looking at %d tICS" %num_tics
	for i in range(0,num_tics):
		for j in range(i+1,num_tics):
			plot_tica_and_clusters(tica_dir, transformed_data_dir, clusterer_dir, lag_time, component_i = i, component_j = j, cluster_ids = cluster_ids)
	print "Printed all tICA coords and all requested clusters"



def plot_timescales(clusterer_dir, n_clusters, lag_time):
	clusterer = verboseload(clusterer_dir)
	sequences = clusterer.labels_
	lag_times = list(np.arange(1,150,5))
	n_timescales = 5

	msm_timescales = implied_timescales(sequences, lag_times, n_timescales=n_timescales, msm=MarkovStateModel(verbose=False))
	print msm_timescales

	for i in range(n_timescales):
		plt.plot(lag_times, msm_timescales[:,i])
	plt.semilogy()
	pp = PdfPages("/scratch/users/enf/b2ar_analysis/kmeans_%d_%d_implied_timescales.pdf" %(n_clusters, lag_time))
	pp.savefig()
	pp.close()

def rmsd_to_structure(clusters_dir, ref_dir, text):
	pdbs = get_trajectory_files(clusters_dir)

	ref = md.load_frame(ref_dir, index=0)
	rmsds = np.zeros(shape=(len(pdbs),2))

	for i in range(0,len(pdbs)):
		print i 
		pdb_file = pdbs[i]
		pdb = md.load_frame(pdb_file, index=0)
		rmsd = md.rmsd(pdb, ref, 0)
		rmsds[i,0] = i
		rmsds[i,1] = rmsd[0]

	rmsd_file = "%s/%s_rmsds.csv" %(clusters_dir, text)
	np.savetxt(rmsd_file, rmsds, delimiter=",")

def rmsd_pymol(pdb_dir, ref_dir, script_dir, rmsd_dir):
	script = open(script_dir, "rb")
	lines = script.readlines()

	new_script = open(script_dir, "wb")

	for line in lines:
		if line[0:7] == "pdb_dir": 
			print ("found pdb line")
			line = "pdb_dir = '%s'\n" %pdb_dir
		elif line[0:7] == "ref_dir": line = "ref_dir = '%s'\n" %ref_dir
		elif line[0:8] == "rmsd_dir": line = "rmsd_dir = '%s'\n" %rmsd_dir
		new_script.write(line)

	new_script.close()
	command = "/scratch/users/enf/pymol/pymol %s" %script_dir
	print command
	os.system(command)

def analyze_rmsds(inactive_rmsd_file, active_rmsd_file, pnas_i, pnas_a, combined_file, analysis_file):
	inactive_rmsd_map = convert_csv_to_map(inactive_rmsd_file)
	inactive_stats_map = calc_mean_and_stdev(inactive_rmsd_map)

	active_rmsd_map = convert_csv_to_map(active_rmsd_file)
	active_stats_map = calc_mean_and_stdev(active_rmsd_map)

	pnas_map_i = convert_csv_to_map(pnas_i)
	pnas_stats_i = calc_mean_and_stdev(pnas_map_i)

	pnas_map_a = convert_csv_to_map(pnas_a)
	pnas_stats_a = calc_mean_and_stdev(pnas_map_a)

	rmsd_i_map = convert_csv_to_map_nocombine(inactive_rmsd_file)
	rmsd_a_map = convert_csv_to_map_nocombine(active_rmsd_file)
	dist_i_map = convert_csv_to_map_nocombine(pnas_i)
	dist_a_map = convert_csv_to_map_nocombine(pnas_a)

	new_file = open(analysis_file, "wb")
	new_file.write("cluster, inactive_rmsd, inactive_stdev, inactive_pnas_dist, inactive_pnas_stdev, active_rmsd, active_stdev, active_pnas_dist, active_pnas_stdev \n")

	combined = open(combined_file, "wb")
	combined.write("cluster, inactive_rmsd, inactive_pnas_dist, active_rmsd, active_pnas_dist \n")

	for key in sorted(inactive_rmsd_map.keys()):
		new_file.write("%s, %f, %f, %f, %f, %f, %f, %f, %f \n" %(key, inactive_stats_map[key][0], inactive_stats_map[key][1], pnas_stats_i[key][0], pnas_stats_i[key][1], active_stats_map[key][0], active_stats_map[key][1], pnas_stats_a[key][0], pnas_stats_a[key][1]))
	for key in sorted(rmsd_i_map.keys()):
		combined.write("%s, %f, %f, %f, %f \n" %(key, rmsd_i_map[key][0], dist_i_map[key][0], rmsd_a_map[key][0], dist_a_map[key][0]))
	new_file.close()
	combined.close()
	
	return [inactive_stats_map, active_stats_map]

def merge_samples(results_map):
	merged_results = {}
	for key in results_map.keys():
		cluster = key.split("_")[0]
		if cluster not in merged_results.keys():
			merged_results[cluster] = [results_map[key]]
		else:
			merged_results[cluster].append(results_map[key])
	stats_map = calc_mean_and_stdev(merged_results)
	return stats_map

def print_stats_map(merged_results, directory):
	mapfile = "%s/stats_map.csv" %directory
	mapcsv = open(mapfile, "wb")
	mapcsv.write("cluster, mean_score, stdev \n")
	for key in sorted(merged_results.keys()):
		mapcsv.write("%s, %f, %f \n" %(key, merged_results[key][0], merged_results[key][1]))
	mapcsv.close()
	return

def analyze_docking_results(docking_dir):
	results_file = "%s/docking_summary.csv" %docking_dir
	results = open(results_file, "wb")
	logs = get_trajectory_files(docking_dir, ext = ".log")
	scores = {}

	for log_file in logs:
		log = open(log_file, "rb")
		conformation = log_file.rsplit(".", 1)[0]
		conformation = conformation.split("/")[len(conformation.split("/"))-1 ]
		score = 1000.0
		xp_score = None
		lines = log.readlines()
		for line in lines:
			line = line.split()
			if len(line) >= 3:
				if (line[0] == "Best" and line[1] == "XP" and line[2] == "pose:"):
					xp_score = float(line[6])
					print "%f, %f" %(xp_score, score)
					if xp_score < score: score = xp_score
				elif  (line[0] == "Best" and line[1] == "Emodel="):
					xp_score = float(line[8])
					print "%f, %f" %(xp_score, score)
					if xp_score < score: score = xp_score
		if score == 1000.0: continue
		scores[conformation] = score

	for conf in sorted(scores.keys()):
		score = scores[conf]
		results.write("%s, %f \n" %(conf, score))

	merged_results = merge_samples(scores)
	print_stats_map(merged_results, docking_dir)

	results.close()

def combine_docking_distances(docking_csv, distances_csv):
	docking_map = convert_csv_to_map_nocombine(docking_csv)
	distances_map = convert_csv_to_map_nocombine()
	combined_map = copy.deepcopy(distances_map)

	for key in distances_map.keys():
		if key in docking_map.keys():
			combined_map[key].append(docking_map[key][0])
	
def combine_rmsd_docking_maps(rmsd_csv, docking_csv):
	rmsd_map = {}
	docking_map = {}

	rmsd_file = open(rmsd_csv, "rb")
	docking_file = open(docking_csv, "rb")

	rmsd_lines = rmsd_file.readlines()
	docking_lines = docking_file.readlines()

	i = 0
	for line in rmsd_lines:
		if i == 0: continue
		line = line.split()
		rmsd_map[line[0]] = []
		j = 0
		for value in line:
			if j == 0: continue
			rmsd_map.append(line[j])

	i = 0
	for line in docking_lines:
		if i == 0: continue
		line = line.split()
		j = 0
		for value in line:
			if j == 0: continue
			docking_map.append(line[j])

	print docking_map[docking_map.keys()[0]]