from io_functions import *
import numpy as np
import copy
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial
import multiprocessing as mp
import mdtraj as md
import csv
import operator
from msmbuilder.utils import verbosedump, verboseload


def calc_mean_and_stdev(rmsd_map):
	stats_map = {}
	for key in rmsd_map.keys():
		rmsds = np.array(rmsd_map[key])
		mean = np.mean(rmsds, axis = 0)
		stdev = np.std(rmsds, axis = 0)
		stats_map[key] = (mean, stdev)
	return stats_map

def calc_mean(rmsd_map):
	stats_map = {}
	for key in rmsd_map.keys():
		rmsds = np.array(rmsd_map[key])
		mean = np.mean(rmsds, axis = 0)
		stats_map[key] = [mean]
	return stats_map

def find_cos(index, k_mean, features):
		traj = index[0]
		frame = index[1]
		conformation = features[traj][frame]
		a = conformation
		b = k_mean
		return (traj, frame, np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)))

def rmsd_connector(traj, inactive, residues_map = None):
	residues = [121, 282]
	if residues_map is not None:
		residues = map_residues(residues_map, residues)

	
	npxxy_atoms = [a.index for a in traj.topology.atoms if a.residue.resSeq in residues and a.is_backbone]
	#print npxxy_atoms
	traj_stripped = traj.atom_slice(npxxy_atoms)

	npxxy_atoms_target = [a.index for a in inactive.topology.atoms if a.residue.resSeq in [121, 282] and a.is_backbone]
	#print npxxy_atoms_target
	inactive_stripped = inactive.atom_slice(npxxy_atoms_target)

	traj_stripped_aligned = traj_stripped.superpose(inactive_stripped)
	rmsds = md.rmsd(traj_stripped, inactive_stripped) * 10.0
	return rmsds

def rmsd_npxxy(traj, inactive, residues_map = None):
	residues = range(322,328)
	if residues_map is not None:
		residues = map_residues(residues_map, residues)

	
	npxxy_atoms = [a.index for a in traj.topology.atoms if a.residue.resSeq in residues and a.is_backbone]
	#print npxxy_atoms
	traj_stripped = traj.atom_slice(npxxy_atoms)

	npxxy_atoms_target = [a.index for a in inactive.topology.atoms if a.residue.resSeq in range(322,328) and a.is_backbone]
	#print npxxy_atoms_target
	inactive_stripped = inactive.atom_slice(npxxy_atoms_target)

	traj_stripped_aligned = traj_stripped.superpose(inactive_stripped)
	rmsds = md.rmsd(traj_stripped, inactive_stripped) * 10.0
	return rmsds

def helix6_helix3_dist(traj, residues_map = None):
	residues = [131, 272]
	if residues_map is not None:
		residues = map_residues(residues_map, residues)

	print "New residues = "
	print [r for r in traj.topology.residues if r.resSeq == residues[1]]

	atom_3 = [a.index for a in traj.topology.atoms if a.residue.resSeq == residues[0] and a.name == "CA"][0]
	atom_6 = [a.index for a in traj.topology.atoms if a.residue.resSeq == residues[1] and a.name == "CA"][0]
	indices = np.empty([1,2])
	indices[0] = [atom_3, atom_6]
	dist = md.compute_distances(traj, indices) * 10.0
	return np.concatenate(dist)



def plot_pnas_vs_docking(docking_dir, pnas_dir, save_dir, selected = False):
	dock_scores = convert_csv_to_map_nocombine(docking_dir)
	pnas_vectors = convert_csv_to_map_nocombine(pnas_dir)
	x = []
	y = []
	c = []
	for key in dock_scores.keys():
		if selected is not False:
			if key not in selected: continue
		#print "PNAS"
		#print key
		#print pnas_vectors[key]
		if key in pnas_vectors.keys():
			#print pnas_vectors[key]
			x.append(pnas_vectors[key][0])
			y.append(pnas_vectors[key][1])
			c.append(abs(dock_scores[key][0]))
			#plt.annotate("%s" %key, xy=(pnas_vectors[key][0], pnas_vectors[key][1]), xytext=(pnas_vectors[key][0], pnas_vectors[key][1]),size=6)

	#print dock_scores.keys()
	plt.scatter(x, y, c=c, s=50, cmap = mpl.cm.RdYlBu_r)
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
	traj_coords = [helix6_helix3_dist(traj, traj_file = traj_file) / scale, rmsd_npxxy(traj, inactive, traj_file = traj_file)]
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



def plot_hex(transformed_data_file, figure_directory, colors = None, scale = 1.0):
	transformed_data = verboseload(transformed_data_file)
	trajs = np.concatenate(transformed_data)
	print trajs
	plt.hexbin(trajs[:,0] * scale, trajs[:,1], bins='log', mincnt=1)
	pp = PdfPages(figure_directory)
	pp.savefig()
	pp.close()
	return

def plot_col(transformed_data_file, figure_directory, colors_file):
	transformed_data = verboseload(transformed_data_file)
	trajs = np.concatenate(transformed_data)
	colors = np.concatenate(verboseload(colors_file))
	sc = plt.scatter(trajs[:,0], trajs[:,1], c=colors, s=50, cmap = mpl.cm.RdYlBu_r)
	plt.colorbar(sc)
	plt.show()
	pp = PdfPages(figure_directory)
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


def plot_tica_and_clusters(tica_dir, transformed_data_dir, clusterer_dir, lag_time, component_i, component_j, cluster_ids = []):
	transformed_data = verboseload(transformed_data_dir)
	clusterer = verboseload(clusterer_dir)

	trajs = np.concatenate(transformed_data)
	plt.hexbin(trajs[:,component_i], trajs[:,component_j], bins='log', mincnt=1)

	centers = clusterer.cluster_centers_
	for i in cluster_ids:
		center = centers[i,:]
		plt.annotate('%d' %i, xy=(center[component_i],center[component_j]), xytext=(center[component_i], center[component_j]),size=6)

	pp = PdfPages("%s/c%d_c%d_clusters%d.pdf" %(tica_dir, component_i, component_j, np.shape(centers)[0]))
	pp.savefig()
	pp.close()
	plt.clf()

def plot_all_tics_and_clusters(tica_dir, transformed_data_dir, clusterer_dir, lag_time, cluster_ids = []):
	transformed_data = verboseload(transformed_data_dir)
	num_tics = np.shape(transformed_data[0])[1]
	print "Looking at %d tICS" %num_tics
	for i in range(0,num_tics):
		for j in range(i+1,num_tics):
			plot_tica_and_clusters(tica_dir = tica_dir, transformed_data_dir = transformed_data_dir, clusterer_dir = clusterer_dir, lag_time = lag_time, component_i = i, component_j = j, cluster_ids = cluster_ids)
	print "Printed all tICA coords and all requested clusters"

def plot_tica_component_i_j(tica_dir, transformed_data_dir, lag_time, component_i = 0, component_j = 1):
	transformed_data = verboseload(transformed_data_dir)

	trajs = np.concatenate(transformed_data)
	plt.hexbin(trajs[:,component_i], trajs[:,component_j], bins='log', mincnt=1)

	pp = PdfPages("%s/c%d_c%d.pdf" %(tica_dir, component_i, component_j))
	pp.savefig()
	pp.close()
	plt.clf()

def plot_all_tics(tica_dir, transformed_data_dir, lag_time):
	transformed_data = verboseload(transformed_data_dir)
	num_tics = np.shape(transformed_data[0])[1]
	print "Looking at %d tICS" %num_tics
	for i in range(0,num_tics):
		for j in range(i+1,num_tics):
			plot_tica_component_i_j(tica_dir, transformed_data_dir, lag_time, component_i = i, component_j = j)
	print "Printed all tICA coords"

#Add way to plot location of specific clusters as well
def plot_all_tics_samples(tica_coords_csv, save_dir, docking_csv = False, specific_clusters = False):
	tica_coords_map = convert_csv_to_map_nocombine(tica_coords_csv)
	n_samples = len(tica_coords_map.keys())
	if docking_csv is not False: docking_map = convert_csv_to_map_nocombine(docking_csv)

	num_tics = len(tica_coords_map[tica_coords_map.keys()[0]])
	for i in range(0, num_tics):
		for j in range(i + 1, num_tics):
			print("plotting tICS %d %d" %(i, j))
			x = []
			y = []
			if docking_csv is not False: c = []
			for sample in tica_coords_map.keys():
				x.append(tica_coords_map[sample][i])
				y.append(tica_coords_map[sample][j])
				if docking_csv is not False: 
					sample_id = "cluster%s" %sample
					print docking_map[sample_id]
					c.append(abs(docking_map[sample_id][0]))
			if docking_csv is not False:
				plt.scatter(x, y, c=c, s=50, cmap = mpl.cm.RdYlBu_r)
			else:
				plt.scatter(x, y, s=50, color = 'red')

			if specific_clusters is not False:
				for cluster in specific_clusters:
					cluster = str(cluster) 
					cluster_x = float(tica_coords_map[cluster][i])
					cluster_y = float(tica_coords_map[cluster][j])
					plt.annotate('%s' %cluster, xy=(cluster_x,cluster_y), xytext=(cluster_x, cluster_y),size=15)

			plot_dir = "%s/samples_%d_tICS_%d_%d.pdf" %(save_dir, n_samples, i, j)
			pp = PdfPages(plot_dir)
			pp.savefig()
			pp.close()
			plt.clf()



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

def analyze_docking_results(docking_dir, ligand, precision, docking_summary):
	results_file = docking_summary
	results = open(results_file, "wb")
	logs = get_trajectory_files(docking_dir, ext = ".log")
	scores = {}

	log_num = 1
	for log_file in logs:
		if log_num % 1000 == 0: print log_num
		log = open(log_file, "rb")
		conformation = log_file.rsplit(".", 1)[0]
		conformation = conformation.split("/")[len(conformation.split("/"))-1 ]
		score = 0.0
		xp_score = None
		lines = log.readlines()
		for line in lines:
			line = line.split()
			if len(line) >= 3:
				if (line[0] == "Best" and line[1] == "XP" and line[2] == "pose:"):
					xp_score = float(line[6])
					#print "%f, %f" %(xp_score, score)
					if xp_score < score: score = xp_score
				elif  (line[0] == "Best" and line[1] == "Emodel="):
					xp_score = float(line[8])
					#print "%f, %f" %(xp_score, score)
					if xp_score < score: score = xp_score
		scores[conformation] = [-1.0*score]
		log_num += 1

	titles = ["sample", "%s" %("%s_%s_score" %(ligand, precision))]

	write_map_to_csv(results_file, scores, titles)

	merged_results = merge_samples(scores)
	return scores 

def analyze_docking_results_wrapper(args):
	return analyze_docking_results(*args)

def get_lig_names(docking_dir):
	subdirs = [x[0] for x in os.walk(docking_dir)]
	lig_names = []

	for subdir in subdirs:
		lig_name = subdir.split("/")[len(subdir.split("/"))-1]
		lig_names.append(lig_name)

	return lig_names

def analyze_docking_results_multiple(docking_dir, precision, ligands, summary):
	subdirs = [x[0] for x in os.walk(docking_dir)]
	subdirs = subdirs[1:len(subdirs)]
	#print subdirs
	results_list = []
	lig_names = []
	arg_tuples = []

	for subdir in subdirs:
		lig_name = subdir.split("/")[len(subdir.split("/"))-1]
		if lig_name not in ligands: continue
		lig_names.append(lig_name)
		docking_summary = "%s/docking_summary.csv" %subdir
		arg_tuples.append([subdir, lig_name, precision, docking_summary])

	pool = mp.Pool(mp.cpu_count())
	results_list = pool.map(analyze_docking_results_wrapper, arg_tuples)
	pool.terminate()

	combined_map = combine_maps(results_list)
	combined_filename = summary
	write_map_to_csv(combined_filename, combined_map, ["sample"] + lig_names)

def compute_means(docking_csv, joined_csv, means_csv):
	print "analyzing %s" %docking_csv
	titles = get_titles(docking_csv)
	docking_scores = convert_csv_to_joined_map(docking_csv, joined_csv)[0]
	docking_averages = calc_mean(docking_scores)
	write_map_to_csv(means_csv, docking_averages, titles)
	return docking_averages

def compute_means_ligands(docking_dir, pnas_means, ligands):
	subdirs = [x[0] for x in os.walk(docking_dir)]
	subdirs = subdirs[1:len(subdirs)]

	docking_csv_files = []
	docking_means_files = []
	docking_joined_files = []

	for subdir in subdirs:
		lig = subdir.split("/")[len(subdir.split("/"))-1]
		if lig not in ligands: continue
		docking_csv = "%s/docking_summary.csv" %subdir
		docking_means = "%s/docking_means.csv" %subdir
		docking_joined = "%s/docking_joined.csv" %subdir
		compute_means(docking_csv, docking_joined, docking_means)
		pnas_joined_csv = "%s/docking_means_pnas_means.csv" %subdir
		combine_csv_list([docking_means, pnas_means], pnas_joined_csv)

	return 

def compute_z(value, mean, stdev):
	return (value - mean) / stdev

def compute_aggregate_scores(docking_csv, inverse_agonists = [], summary = "", z_scores_csv = ""):
	scores_map = convert_csv_to_map_nocombine(docking_csv)
	docking_titles = get_titles(docking_csv)
	lig_names = docking_titles[1:len(docking_titles)]
	print lig_names
	scores_per_ligand = {}
	for lig_name in lig_names:
		scores_per_ligand[lig_name] = []
	for receptor in scores_map.keys():
		receptor_scores = scores_map[receptor]
		for i in range(0, len(receptor_scores)):
			lig_name = lig_names[i]
			lig_score = receptor_scores[i]
			scores_per_ligand[lig_name].append(lig_score)
	ligand_stats = calc_mean_and_stdev(scores_per_ligand)

	z_scores_per_receptor = {}

	for receptor in scores_map.keys():
		receptor_scores = scores_map[receptor]
		z_scores_per_receptor[receptor] = []
		for i in range(0, len(receptor_scores)):
			lig_score = receptor_scores[i]
			lig_name = lig_names[i]
			lig_mean = ligand_stats[lig_name][0]
			lig_stdev = ligand_stats[lig_name][1]
			z_score = compute_z(lig_score, lig_mean, lig_stdev)
			if abs(lig_score - 0.0) < 0.001: z_score = -3.0
			if lig_name in inverse_agonists:
				#print lig_name
				z_score = -1.0 * z_score

			z_scores_per_receptor[receptor].append(z_score)

	aggregate_scores = calc_mean(z_scores_per_receptor)

	titles = ["sample", "average_z_score"]

	write_map_to_csv(z_scores_csv, z_scores_per_receptor, docking_titles)
	write_map_to_csv(summary, aggregate_scores, titles)

def combine_docking_distances(docking_csv, distances_csv, docking_dir):
	docking_map = convert_csv_to_map_nocombine(docking_csv)
	distances_map = convert_csv_to_map_nocombine(distances_csv)
	combined_map = copy.deepcopy(distances_map)

	for key in distances_map.keys():
		if key in docking_map.keys():
			combined_map[key].append(-1.0 * docking_map[key][0])
		else:
			combined_map.pop(key, None)
	
	firstline = "cluster, inactive_rmsd, inactive_pnas, active_rmsd, active_pnas, docking \n"
	filename = "%s/distances_docking.csv" %docking_dir
	write_map_to_csv(filename, combined_map, firstline)

def top_n_scoring_clusters(docking_csv, score_type = 1, n = 50):
	docking_list = convert_csv_to_list(docking_csv)

	docking_list_sorted = sorted(docking_list, key=operator.itemgetter(score_type), reverse = True)

	top_n = []
	for i in range(0, n):
		top_n.append(docking_list_sorted[i][0])

	return top_n

def top_n_scoring_samples(docking_csv, score_type = "mean_docking_score", n = 100, n_samples = 10):
	docking_list = convert_csv_to_list(docking_csv)
	titles = docking_list[0]
	for i in range(0, len(titles)):
		if titles[i] == score_type: break
	#print i
	#print titles[i]
	docking_list_sorted = sorted(docking_list, key=operator.itemgetter(i), reverse = True)
	#print docking_list_sorted
	top_n = []
	for i in range(0, n):
		top_n.append(docking_list_sorted[i][0])

	if "sample" not in top_n[0]:
		top_n_new = []
		for cluster in top_n:
			for i in range(0,n_samples):
				top_n_new.append("%s_sample%d" %(cluster, i))
		top_n = top_n_new

	return top_n

def combine_docking_mmgbsa(combined_csv, mmgbsa_csv, mmgbsa_dir, filename):
	combined_map = convert_csv_to_map_nocombine(combined_csv)
	mmgbsa_map = convert_csv_to_map_nocombine(mmgbsa_csv)
	new_map = copy.deepcopy(combined_map)

	for key in mmgbsa_map.keys():
		if key in combined_map.keys():
			new_map[key].append(-1.0 * mmgbsa_map[key][0])
		else:
			new_map.pop(key, None)

	firstline = "sample, inactive_rmsd, inactive_pnas, active_rmsd, active_pnas, docking, mmgbsa \n"
	write_map_to_csv(filename, new_map, firstline)

def analyze_mmgbsa_results(mmgbsa_dir, ligand, chosen_receptors):
	analysis_csv = "%s/mmgbsa_summary.csv" %mmgbsa_dir
	csvfile = open(analysis_csv, "wb")
	csvfile.write("sample, %s_mmgbsa \n" %ligand)
	outputs = get_trajectory_files(mmgbsa_dir, ".csv")
	analyzed_samples = set()

	for output in outputs:
		output_name = output.split("/")[len(output.split("/"))-1]
		sample = output_name.split("-out.csv")[0]
		if chosen_receptors is not False:
			if sample not in chosen_receptors: continue
		analyzed_samples.add(sample)
		results = open(output, "rt")
		reader = csv.reader(results)
		score = 10000.0
		i = 0
		for row in reader:
			if i == 0:
				i += 1
				continue
			print row[0]
			temp = float(row[1])
			if temp < score:
				score = temp
		score = -1.0 * score
		csvfile.write("%s, %f \n" %(sample, score))
		results.close()

	not_analyzed_samples = set(chosen_receptors) - analyzed_samples
	for sample in not_analyzed_samples:
		if sample not in analyzed_samples:
			csvfile.write("%s, %f \n" %(sample, 0.00))


	csvfile.close()

	results_map = convert_csv_to_map_nocombine(analysis_csv)
	return results_map


def analyze_mmgbsa_results_multiple(mmgbsa_dir, summary, ligands, chosen_receptors):
	subdirs = [x[0] for x in os.walk(mmgbsa_dir)]
	subdirs = subdirs[1:len(subdirs)]
	#print subdirs
	results_list = []
	lig_names = []

	for subdir in subdirs:
		lig_name = subdir.split("/")[len(subdir.split("/"))-1]
		if lig_name not in ligands: continue
		lig_names.append(lig_name)
		mmgbsa_summary = "%s/mmgbsa_summary.csv" %subdir
		results = analyze_mmgbsa_results(subdir, ligand = lig_name, chosen_receptors = chosen_receptors)
		#print results
		#print subdir
		results_list.append(results)
		
	combined_map = combine_maps(results_list)
	combined_filename = summary
	write_map_to_csv(combined_filename, combined_map, ["sample"] + lig_names)

def keep_only_sample0(original_csv, output_csv):
	original = open(original_csv, "rb")
	output = open(output_csv, "wb")
	lines = original.readlines()
	i = 0
	for line in lines:
		if i == 0:
			output.write(line)
			i += 1
		else:
			if "sample0" in line:
				output.write(line)

	original.close()
	output.close()


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