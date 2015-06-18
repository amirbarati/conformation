import fileinput
import os
import copy
import csv
from msmbuilder.utils import verbosedump, verboseload
import numpy as np

def get_trajectory_files(traj_dir, ext = ".pdb"):
	traj_files = []
	for traj in os.listdir(traj_dir):
			if traj.endswith(ext):
				traj_files.append("%s/%s" %(traj_dir,traj))
	return sorted(traj_files)

def get_trajectory_files_conditions(traj_dir, ext, condition_1, condition_2):
	trajs = get_trajectory_files(traj_dir, ext)
	traj_1 = [t for t in trajs if condition_1 in t]
	traj_2 = [t for t in trajs if condition_2 in t]
	traj_1 = sorted(traj_1)
	traj_2 = sorted(traj_2)
	print len(traj_1)
	print(len(traj_2))
	print(len(traj_1+traj_2))
	return(traj_1 + traj_2)

def get_ligands(lig_dir, ext= ".sdf"):
	ligands = get_trajectory_files(lig_dir, ext = ext)
	ligs = []
	for ligand in ligands:
		lig_last_name = ligand.split("/")[len(ligand.split("/"))-1]
		lig = lig_last_name.split(".")[0]
		ligs.append(lig)
	return ligs

def write_map_to_csv(filename, data_map, titles):
	csvfile = open(filename, "wb")
	i = 0
	for title in titles:
		if i < (len(titles) - 1):
			csvfile.write("%s, " %title)
		else:
			csvfile.write("%s \n" %title)
		i += 1

	for key in sorted(data_map.keys()):
		csvfile.write("%s, " %key)
		i = 0
		for value in data_map[key]:
			if i < (len(data_map[key])-1):
				csvfile.write("%f, " %value)
			else:
				csvfile.write("%f \n" %value)
			i += 1
	return

def generateData(files):
	for data in files:
		print data
		yield verboseload(data)

def generateTraj(files, top=None):
	for traj in files:
		if top == None:
			yield md.load(file)
		else:
			yield md.load(file, top=top)

def get_trajectory_files(traj_dir, ext = ".h5"):
	traj_files = []
	for traj in os.listdir(traj_dir):
			if traj.endswith(ext):
				traj_files.append("%s/%s" %(traj_dir,traj))
	return sorted(traj_files)

def get_traj_no_palm(traj_dir):
	trajs = get_trajectory_files(traj_dir)
	non_palm = []

	for i in range(0, len(trajs)):
		traj = trajs[i]
		traj_name = traj.split("/")[len(traj.split("/"))-1]
		if traj_name[0] not in ['C', 'F', 'H', 'J']:
			non_palm.append(i)

	return non_palm

def convert_csv_to_list(filename):
	with open(filename, "rb") as f:
		reader = csv.reader(f)
		lines = list(reader)
		i = 0
		lines.pop(0)
		for line in lines:
			for i in range(1, len(line)):
				line[i] = float(line[i])
		return lines

def convert_csv_to_map_nocombine(filename):
	print filename
	rmsds = open(filename, "rb")
	lines = [line.strip() for line in rmsds.readlines()]
	#if "docking" in filename:
	#	print lines
	rmsd_map = {}
	i = 0
	for line in lines:
		if i == 0:
			i += 1
			continue
		if ';' in line:
			line = line.split(';')
		elif ',' in line:
			line = line.split(',')
		else:
			line = line.split()

		#print line
		cluster = line[0].split('.')[0]
		#print cluster
		#print line
		for i in range(1,len(line)):
			try:
				rmsd = float(line[i].split('\\')[0])
			except:
				continue
			#print rmsd
			#if rmsd > 3.0: print "%s %f" %(cluster, rmsd)

			if cluster in rmsd_map.keys():
				rmsd_map[cluster].append(rmsd)
			else:
				rmsd_map[cluster] = [rmsd]
	return rmsd_map

def get_titles(filename):
	csv_file = open(filename, "rb")
	firstline = csv_file.readlines()[0].split("\n")[0].split(",")
	firstline = [f.strip() for f in firstline]
	return firstline

def convert_csv_to_joined_map(filename, new_filename = False):
	rmsds = open(filename, "rb")
	lines = rmsds.readlines()
	rmsd_map = {}
	i = 0
	for line in lines:
		if i == 0:
			i += 1
			continue
		if ';' in line:
			line = line.split(';')
		elif ',' in line:
			line = line.split(',')
		else:
			line = line.split()



		#print line
		cluster = line[0].split('_')[0]
		rmsd = float(line[1].split('\\')[0])

		if cluster in rmsd_map.keys():
			rmsd_map[cluster].append(rmsd)
		else:
			rmsd_map[cluster] = [rmsd]

	titles = []
	num_entries = len(rmsd_map[rmsd_map.keys()[0]])
	titles.append("cluster")
	for i in range(0, num_entries):
		titles.append("sample%d" %i)

	if new_filename is not False:
		write_map_to_csv(new_filename, rmsd_map, titles)

	return [rmsd_map, titles]

def combine_map(map_i, map_j):
	map_j_val_length = len(map_j[map_j.keys()[0]])
	placeholder = []
	for i in range(0, map_j_val_length):
		placeholder.append(0.0)

	for key in map_i.keys():
		if key in map_j.keys():
			map_i[key] += map_j[key]
		else:
			#map_i[key] += placeholder
			map_i.pop(key)
	return map_i

def combine_maps(map_list):
	combined_map = copy.deepcopy(map_list[0])
	for i in range(1, len(map_list)):
		combined_map = combine_map(combined_map, map_list[i])
	return combined_map

def combine_csv_list(csv_list, new_csv_filename):
	map_list = []
	combined_map = {}
	combined_map_titles = ["cluster"]

	for csv in csv_list:
		map_list.append(convert_csv_to_map_nocombine(csv))
		csv_file = open(csv, "rb")
		csv_titles = csv_file.readlines()[0].split("\n")[0].split(",")
		csv_titles = csv_titles[1:len(csv_titles)]
		print csv_titles
		combined_map_titles += csv_titles

	combined_map = combine_maps(map_list)
	print combined_map_titles
	write_map_to_csv(new_csv_filename, combined_map, combined_map_titles)
	return combined_map

def get_condition(filename):
	filename = filename.split('/')[len(filename.split('/'))-1]
	pieces = filename.split('-')
	condition = "%s-%s" %(pieces[0], pieces[1])
	return condition 

def read_trajectory(directory, filename, stride=10):
	print("reading %s" %(filename))
	traj = md.load(filename, stride=stride, top="/home/harrigan/compute/wetmsm/gpcr/des/system_mae_to_pdb/des_trajs/DESRES-Trajectory_pnas2011b-H-05-all/system.pdb")
	return traj

def reverse_sign_csv(csv_file):
	with open(csv_file, 'rb') as f:
		reader = csv.reader(f)
		lines = list(reader)

	new_csv = open(csv_file, "wb")
	line_num = 0
	for line in lines:
		if line_num > 0:
			for i in range(1, len(line)):
				line[i] = str(-1.0 * (float(line[i])))
			print line
			new_csv.write(",".join(line))
			new_csv.write(" \n")
		else:
			new_csv.write(",".join(line))
			new_csv.write(" \n")
			line_num += 1
	new_csv.close()

def convert_matrix_to_map(matrix_file, traj_dir, ext, header, csv_file):
	trajs = get_trajectory_files(traj_dir, ext = ext)
	matrix = np.vstack(verboseload(matrix_file))
	values_map = {}
	for i in range(0, np.shape(matrix)[0]):
		traj = trajs[i]
		traj_lastname = traj.split("/")[len(traj.split("/"))-1]
		traj_name = traj_lastname.split(".")[0]
		values = tuple(matrix[i,:])
		values_map[traj_name] = values
		write_map_to_csv(csv_file, values_map, header)
	return values_map


def convert_matrix_list_to_list(np_file, csv_file):
	matrix_list = verboseload(np_file)
	all_values = np.concatenate(matrix_list)
	np.savetxt(csv_file, all_values, delimiter=",")
	return all_values	 

def find_missing_features(traj_dir, features_dir):
	trajs = get_trajectory_files(traj_dir, ".lh5")
	trajs = [t.split("/")[len(t.split("/"))-1].split(".")[0] for t in trajs]
	trajs = set(trajs)
	features = get_trajectory_files(features_dir, ".h5")
	features = [f.split("/")[len(f.split("/"))-1].split(".")[0] for f in features]
	features = set(features)
	print(trajs - features)

def generate_residues_map(csv_map):
	reader = csv.reader(open(csv_map, "rb"))
	residues_map = {}
	for line in reader:
		residues_map[int(line[0])] = int(line[1])
	return residues_map

def map_residues(residues_map, residues):
	new_residues = []
	for residue in residues:
		new_residues.append(residues_map[residue])
	return new_residues 

def map_residues_universal(residues, save):
	pdb_file = "/home/enf/b2ar_analysis_sherlock_all/exacycle_data/alignment_universal_pdb.txt"
	lh5_file = "/home/enf/b2ar_analysis_sherlock_all/exacycle_data/alignment_universal_lh5.txt"

	pdb_lines = []
	lh5_lines = []

	pdb = open(pdb_file, "rb")
	lh5 = open(lh5_file, "rb")
	for line in pdb.readlines():
		if "TER" not in line and "END" not in line:
			pdb_lines.append(line.split())

	for line in lh5.readlines():
		if "TER" not in line and "END" not in line:
			lh5_lines.append(line.split())

	if len(pdb_lines) != len(lh5_lines):
		print("Alignemnt no goood")
		sys.exit()

	residues_map = {}
	new_residues = []

	for i in range(0,len(pdb_lines)):
		if pdb_lines[i][4] == lh5_lines[i][4]:
			residues_map[int(pdb_lines[i][5])] = int(lh5_lines[i][5])

	for residue in residues:
		new_residues.append(residues_map[residue])

	if 281 in residues_map.keys():
		print("Residue 281 is now Residue %d" %residues_map[281])

	print("residue 272 is now residue %d" %residues_map[272])

	pdb.close()
	lh5.close()

	if save != False:
		writer = csv.writer(open(save, "wb"))
		for key, value in residues_map.items():
			writer.writerow([key,value])

	return new_residues

def map_residues_condition(residues, condition):
	if condition == "2rh1":
		pdb_file = "/home/enf/b2ar_analysis_sherlock_all/exacycle_data/align_2rh1_pdb.txt"
		lh5_file = "/home/enf/b2ar_analysis_sherlock_all/exacycle_data/align_2rh1_lh5.txt"
	elif condition == "3p0g":
		pdb_file = "/home/enf/b2ar_analysis_sherlock_all/exacycle_data/align_3p0g_pdb.txt"
		lh5_file = "/home/enf/b2ar_analysis_sherlock_all/exacycle_data/align_3p0g_lh5.txt"

	pdb_lines = []
	lh5_lines = []

	pdb = open(pdb_file, "rb")
	lh5 = open(lh5_file, "rb")
	for line in pdb.readlines():
		if "TER" not in line and "END" not in line:
			pdb_lines.append(line.split())

	for line in lh5.readlines():
		if "TER" not in line and "END" not in line:
			lh5_lines.append(line.split())

	if len(pdb_lines) != len(lh5_lines):
		print("Alignemnt no goood")
		sys.exit()

	residues_map = {}
	new_residues = []

	for i in range(0,len(pdb_lines)):
		if pdb_lines[i][4] == lh5_lines[i][4]:
			residues_map[int(pdb_lines[i][5])] = int(lh5_lines[i][5])

	for residue in residues:
		new_residues.append(residues_map[residue])

	if 281 in residues_map.keys():
		print("Residue 281 is now Residue %d" %residues_map[281])

	print("residue 272 is now residue %d" %residues_map[272])

	pdb.close()
	lh5.close()

	if condition == "2rh1":
		writer = csv.writer(open("/home/enf/b2ar_analysis_sherlock_all/exacycle_data/residues_map_2rh1.csv", "wb"))
		for key, value in residues_map.items():
			writer.writerow([key,value])
	if condition == "3p0g":
		writer = csv.writer(open("/home/enf/b2ar_analysis_sherlock_all/exacycle_data/residues_map_3p0g.csv", "wb"))
		for key, value in residues_map.items():
			writer.writerow([key,value])

	return new_residues





