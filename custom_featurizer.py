import mdtraj as md
import numpy as np
from analysis import *
from msmbuilder.utils import verbosedump, verboseload
import time
from pdb_editing import *
from mdtraj.geometry import dihedral as ManualDihedral
import sys
from msmbuilder.featurizer import DihedralFeaturizer
import itertools
from numpy import random
import json


def fix_topology(topology):
	
	new_top = topology.copy()

	residues = {}
	for chain in new_top.chains:
		#print chain
		for residue in chain.residues:
			resname = str(residue)
			if resname in residues.keys():
				residues[resname].append(residue)
			else:
				residues[resname] = [residue]

	for resname in residues.keys():
		fragments = residues[resname]
		if len(fragments) > 1:
			main_fragment = fragments[0]
			new_atom_list = []
			new_atom_list += main_fragment._atoms
			for i in range(1,len(fragments)):
				fragment = fragments[i]
				for atom in fragment.atoms:
					atom.residue = main_fragment
				new_atom_list += fragment._atoms
				fragment._atoms = []
				fragment.chain = main_fragment.chain
			main_fragment._atoms = new_atom_list

	return new_top

def fix_traj(traj):
	time0 = time.time()
	new_traj = copy.deepcopy(traj)
	topology = new_traj.topology 

	new_top = fix_topology(topology)
	topology = new_top
	new_traj.topology = new_top 

	new_atom_sequence = [a for a in topology.atoms]
	new_index_sequence = [a.index for a in topology.atoms]
	
	for i in range(0, np.shape(traj.xyz)[0]):
		new_traj.xyz[i] = new_traj.xyz[i][new_index_sequence]

	for i in range(0, len(new_index_sequence)):
		new_atom_sequence[i].index = i

	time1 = time.time()
	print time1 - time0
	return new_traj

def phi_indices(top, residues = None):
	residues = copy.deepcopy(residues)
	graph = top.to_bondgraph()

	if residues is None:
		c_atoms = [(a, a.residue.resSeq) for a in top.atoms if a.name == "C"]
	else:
		for i in range(0,len(residues)):
			residues[i] -= 1
		c_atoms = [(a, a.residue.resSeq) for a in top.atoms if a.name == "C" and a.residue.resSeq in residues]

	c_atoms.sort(key=operator.itemgetter(1))
	c_atoms = [c_atom[0] for c_atom in c_atoms]
	#print("%d C atoms" %len(c_atoms))

	phi_tuples = []

	for c in c_atoms:
		n = None
		ca = None
		next_c = None

		c_index = c.index
		c_neighbors = graph.edge[c].keys()
		for c_neighbor in c_neighbors:
			if c_neighbor.name == "N":
				n = c_neighbor
				break
		
		if n != None:
			n_neighbors = graph.edge[n].keys()
			for n_neighbor in n_neighbors:
				if n_neighbor.name == "CA":
						ca = n_neighbor
						break
		if ca != None:
			ca_neighbors = graph.edge[ca].keys()
			for ca_neighbor in ca_neighbors:
				if ca_neighbor.name == "C":
					next_c = ca_neighbor
					break
		if n != None and ca != None and next_c != None:
			phi_tuples.append((c.index, n.index, ca.index, next_c.index))
		else:
			print "No phi found for %s " %c.name

	#print("phi angles = %d" %len(phi_tuples))
	return phi_tuples




def psi_indices(top, residues = None):

	graph = top.to_bondgraph()
	if residues is None:
		n_atoms = [(a, a.residue.resSeq) for a in top.atoms if a.name == "N"]
	else:
		n_atoms = [(a, a.residue.resSeq) for a in top.atoms if a.name == "N" and a.residue.resSeq in residues]

	n_atoms.sort(key=operator.itemgetter(1))
	n_atoms = [n_atom[0] for n_atom in n_atoms]

	psi_tuples = []

	for n in n_atoms:
		c = None
		ca = None
		next_n = None

		n_index = n.index
		n_neighbors = graph.edge[n].keys()
		for n_neighbor in n_neighbors:
			if n_neighbor.name == "CA":
				ca = n_neighbor
				break
		
		if ca != None:	
			ca_neighbors = graph.edge[ca].keys()
			for ca_neighbor in ca_neighbors:
				if ca_neighbor.name == "C":
					c = ca_neighbor
					break

		if c != None:
			c_neighbors = graph.edge[c].keys()
			for c_neighbor in c_neighbors:
				if c_neighbor.name == "N":
					next_n = c_neighbor
					break

		if c != None and ca != None and next_n != None:
			psi_tuples.append((n.index, c.index, ca.index, next_n.index))
		else:
			print "No psis found for %s " %c.residue

	#print("psi angles = %d " %len(psi_tuples))
	return psi_tuples

def phi_indices_resSeq(top):
	'''
	for i in residues
		residue_i = residues[i]
		residue_ip1 = residues[i+1]
		if residue_i.resSeq == residue_ip1.resSeq - 1:
			N = bla
			C = bla
			CA = 
			N_next 
	'''
	return

def chi1_indices(top, specified_residues = None):
	term_4 = ('CG', 'CG1', 'OG1', 'SG', 'OG')
	chi1_residues = ["Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]
	chi1_residues = [a.upper() for a in chi1_residues]

	#top = fix_topology(top)
	if specified_residues is None:
		residues = [(res, res.resSeq) for res in top.residues]
	else:
		residues = [(res, res.resSeq) for res in top.residues if res.resSeq in specified_residues]

	residues.sort(key=operator.itemgetter(1))
	residues = [res[0] for res in residues]
	chi1_tuples = []

	#print "CHI1: \n"
	for residue in residues:
		dihedral = [None, None, None, None]
		for atom in residue.atoms:
			if atom.name == 'N': dihedral[0] = atom.index
			if atom.name == 'CA': dihedral[1] = atom.index
			if atom.name == 'CB': dihedral[2] = atom.index
			if atom.name in term_4: dihedral[3] = atom.index
		if None not in dihedral:
			dihedral = tuple(dihedral)
			chi1_tuples.append(dihedral)
			#print residue.resSeq
		elif dihedral != [None, None, None, None] and str(residue.name)[0:3] in chi1_residues:
			print "no chi1 found for %s" %str(residue)	
	return chi1_tuples



def chi2_indices(top, specified_residues = None):
	seq1 = ('CA', 'CB', 'CG', 'CD')
	seq2 = ('CA', 'CB', 'CG', 'OD1')
	seq3 = ('CA', 'CB', 'CG', 'ND1')
	seq4 = ('CA', 'CB', 'CG1', 'CD1')
	seq5 = ('CA', 'CB', 'CG,' 'SD')

	chi2_residues = ["Arg", "Asn", "Asp", "Gln", "Glu", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Trp", "Tyr"]
	chi2_residues = [a.upper() for a in chi2_residues]

	term_4 = ('CD', 'OD1', 'ND1', 'CD1', 'SD')

	#top = fix_topology(top)
	if specified_residues is None:
		residues = [(res, res.resSeq) for res in top.residues]
	else:
		residues = [(res, res.resSeq) for res in top.residues if res.resSeq in specified_residues]

	residues.sort(key=operator.itemgetter(1))
	residues = [res[0] for res in residues]
	chi2_tuples = []

	#print "CHI2: \n"
	for residue in residues:
		dihedral = [None, None, None, None]
		for atom in residue.atoms:
			if atom.name == 'CA': dihedral[0] = atom.index
			if atom.name == 'CB': dihedral[1] = atom.index
			if atom.name == 'CG' or atom.name == 'CG1': dihedral[2] = atom.index
			if atom.name in term_4: dihedral[3] = atom.index
		if (None not in dihedral) and (str(residue.name)[0:3] in chi2_residues):
			dihedral = tuple(dihedral)
			chi2_tuples.append(dihedral)
			#print residue
		elif dihedral != [None, None, None, None] and str(residue.name)[0:3] in chi2_residues:
			print "no chi2 found for %s" %str(residue)		


	return chi2_tuples

def convert_resSeq_to_resIndex(top, residue_pairs):
	resIndices = []
	for pair in residue_pairs:
		resSeq0 = pair[0][1]
		chain0 = pair[0][0]
		resSeq1 = pair[1][1]
		chain1 = pair[1][0]
		indices = [(r.index,r) for r in top.residues if r.resSeq == resSeq0 and not r.is_water and r.chain.index == chain0]
		if len(indices) == 0:
			print("FATAL: No residues in trajectory for residue %d" %resSeq0)
			return None
		else:
			ind0 = indices[0][0]
			for j in indices:
				if j[0] != ind0: 
					#print("Warning: multiple res objects for residue %d " %resSeq0)
					if "CB" in [str(a) for a in r.atoms for r in top.residues if r.index == ind0]:
						ind0 = j[0]
		indices = [(r.index,r) for r in top.residues if r.resSeq == resSeq1 and not r.is_water and r.chain.index == chain1]
		if len(indices) == 0:
			print("FATAL: No residues in trajectory for residue %d" %resSeq1)
			return None
		else:
			ind1 = indices[0][0]
			for j in indices:
				if j[0] != ind1: 
					#print("Warning: multiple res objects for residue %d " %resSeq1)
					if "CB" in [str(a) for a in r.atoms for r in top.residues if r.index == ind1]:
						ind1 = j[0]
		resIndices.append((ind0, ind1))
	print("looking at %d pairs for trajectory" %len(resIndices))
	return(resIndices)



def read_and_featurize(traj_file, features_dir = None, condition=None, dihedral_types = ["phi", "psi", "chi1", "chi2"], dihedral_residues = None, resSeq_pairs = None, iterative = True):

	a = time.time()
	dihedral_indices = []
	residue_order = []
	if len(dihedral_residues) > 0:
		for dihedral_type in dihedral_types:
			if dihedral_type == "phi": dihedral_indices.append(phi_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "psi": dihedral_indices.append(psi_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "chi1": dihedral_indices.append(chi1_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "chi2": dihedral_indices.append(chi2_indices(fix_topology(top), dihedral_residues))

		#print("new features has dim %d" %(2*len(phi_tuples) + 2*len(psi_tuples) + 2*len(chi2_tuples)))

		#print("feauturizing manually:")
		dihedral_angles = []

		for dihedral_type in dihedral_indices:
			angles = np.transpose(ManualDihedral.compute_dihedrals(traj=traj,indices=dihedral_type))
			dihedral_angles.append(np.sin(angles))
			dihedral_angles.append(np.cos(angles))

		manual_features = np.transpose(np.concatenate(dihedral_angles))

	if len(resSeq_pairs) > 0:
		top = md.load_frame(traj_file, index=0).topology
		resIndex_pairs = convert_resSeq_to_resIndex(top, resSeq_pairs)
		contact_features = []
		if iterative:
			try:
				for chunk in md.iterload(traj_file, chunk = 1000):
				#	chunk = fix_traj(chunk)
				#chunk = md.load(traj_file,stride=1000)
				#print(resIndex_pairs[0:10])
					chunk_features = md.compute_contacts(chunk, contacts = resIndex_pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
					print(np.shape(chunk_features))
					contact_features.append(chunk_features)
				contact_features = np.concatenate(contact_features)
			except Exception,e:
				print str(e)
				print("Failed")
				return
				#traj = md.load(traj_file)
				#contact_features = md.compute_contacts(chunk, contacts = contact_residue_pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
		else:
			try:
				traj = md.load(traj_file)
				contact_features =  md.compute_contacts(traj, contacts = resIndex_pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
			except Exception,e:
				print str(e)
				print("Failed for traj")
				return
		if len(dihedral_residues) > 0: 
			manual_features = np.column_stack((manual_features, contact_features))
		else:
			manual_features = contact_features


	b = time.time()

	print("new features %s has shape: " %traj_file)
	print(np.shape(manual_features))

	if condition is None:
		condition = traj_file.split("/")[len(traj_file.split("/"))-1]
		condition = condition.split(".")[0]

	save_dataset(manual_features, "%s/%s.dataset" %(features_dir, condition))

def compute_contacts_below_cutoff(traj_file_frame, cutoff = 100000.0, contact_residues = [], anton = False):
	traj_file = traj_file_frame[0]
	frame = md.load_frame(traj_file, index = 0)
	#frame = fix_traj(frame)
	top = frame.topology
	
	distance_residues = []
	res_indices = []
	resSeq_to_resIndex = {}
	residue_full_infos = []

	for i in range(0, len(contact_residues)):
		residue = contact_residues[i]
		indices = [r.index for r in top.residues if r.resSeq == residue[1] and r.chainid == residue[0] and not r.is_water]
		if len(indices) == 0:
			print("No residues in trajectory for residue %d" %residue)
			continue
		else:
			ind = indices[0]
			for j in indices:
				if j != ind: 
					#print("Warning: multiple res objects for residue %d " %residue)
					if "CB" in [str(a) for a in r.atoms for r in top.residues if r.index == ind]:
						ind = j
			res_indices.append(ind)
			distance_residues.append(residue)
			resSeq_to_resIndex[residue] = ind
	
	resSeq_combinations = itertools.combinations(distance_residues, 2)
	res_index_combinations = []
	resSeq_pairs = [c for c in resSeq_combinations]
	for combination in resSeq_pairs:
		res0 = combination[0]
		res1 = combination[1]
		res_index0 = resSeq_to_resIndex[res0]
		res_index1 = resSeq_to_resIndex[res1]
		res_index_combinations.append((res_index0, res_index1))


	final_resSeq_pairs = []
	final_resIndex_pairs = []

	distances = md.compute_contacts(frame, contacts = res_index_combinations, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
	#print(distances)
	print(np.shape(distances))
	for i in range(0, len(distances[0])):
		distance = distances[0][i]
		#print(distance)
		if distance < cutoff:
			final_resIndex_pairs.append(res_index_combinations[i])
			final_resSeq_pairs.append(resSeq_pairs[i])

	for pair in final_resIndex_pairs:
		info0 = [(r.resSeq, r.name, r.chain.index) for r in top.residues if r.index == pair[0]]
		info1 = [(r.resSeq, r.name, r.chain.index) for r in top.residues if r.index == pair[1]]
		residue_full_infos.append((info0, info1))

	print(len(final_resSeq_pairs))
	print(len(final_resIndex_pairs))
	return((final_resSeq_pairs, residue_full_infos))


def which_trajs_to_featurize(traj_dir, traj_ext):
	all_trajs = get_trajectory_files(traj_dir, traj_ext)[0:3]
	trajs = []
	for fulltraj in all_trajs:
		#if "H-05" not in fulltraj and "A-00" not in fulltraj: continue
		traj = fulltraj.split("/")
		filename = traj[len(traj)-1]
		#if agonist_bound is not False and filename[0] not in agonist_bound: continue
		filename_noext = filename.split(".")[0]
		if os.path.exists("%s/%s.dataset" %(features_dir, filename_noext)):
			print("already featurized")	
		else:
			trajs.append(fulltraj)
	return(trajs)	

def featurize_contacts_custom(traj_dir, features_dir, traj_ext, featurized_traj_and_frame, contact_residue_pairs_csv = None, dihedral_residues = None, dihedral_types = None, contact_residues = None, agonist_bound = False, residues_map = None, contact_cutoff = None, parallel = False, exacycle = False):
	'''
	residues to be used for featurization need to be input as a list of tuples. Each tuple will be of form (chain, residue_number)

	you can also input a residue_map, a dictionary that maps (chain, residue_number) --> (chain, residue_number). The reason for this is that there
		is not a consensus residue numbering for the same protein.
	'''

	if not os.path.exists(features_dir): os.makedirs(features_dir)
	trajs = which_trajs_to_featurize(traj_dir, traj_ext)

	if residues_map is not None:
		contact_residues = [r for r in contact_residues if r in residues_map.keys()]
		if exacycle: contact_residues = [residues_map[key] for key in contact_residues]	

	print(featurized_traj_and_frame)
	if featurized_traj_and_frame is None: featurized_traj_and_frame = [all_trajs[0], 0]
	if contact_residue_pairs_csv == "" or (not os.path.exists(contact_residue_pairs_csv)):
		contact_residue_pairs, residue_full_infos = compute_contacts_below_cutoff(featurized_traj_and_frame, cutoff = contact_cutoff, contact_residues = contact_residues, anton = False)
	else:
		print("Features already computed")
		time.sleep(10)
		contact_residue_pairs = generate_features(contact_residue_pairs_csv)
		if exacycle: contact_residue_pairs = [residues_map[key] for key in contact_residue_pairs]
	print(contact_residue_pairs)
	print("Number of contact pairs = %d" %len(contact_residue_pairs))


	featurize_partial = partial(read_and_featurize, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = dihedral_types, resSeq_pairs = contact_residue_pairs, iterative = True)
	if parallel:
		pool = mp.Pool(mp.cpu_count())
		pool.map(featurize_partial, trajs)
		pool.terminate()
	else:
		for traj in trajs:
			print traj
			featurize_partial(traj)

	print("Completed featurizing")



def featurize_known_traj(traj_dir, inactive, features_dir):
	print("currently featurizing %s" %traj_dir.split("/")[len(traj_dir.split("/"))-1])
	traj = md.load(traj_dir)
	rmsds = rmsd_npxxy(traj, inactive)
	helix6_helix3_distances = helix6_helix3_dist(traj)
	features = np.transpose(np.concatenate([[rmsds], [np.concatenate(helix6_helix3_distances)]]))
	print np.shape(features)

	filename = "%s/%s" %(features_dir, traj_dir.split("/")[len(traj_dir.split("/"))-1])
	verbosedump(features, filename)

def featurize_known(directory, inactive_dir, active_dir):
	features_dir = "/scratch/users/enf/b2ar_analysis/features_known"
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	ianctive = md.load(inactive_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	all_trajs = get_trajectory_files(directory)
	trajs = []
	for fulltraj in all_trajs:
		traj = fulltraj.split("/")
		filename = traj[len(traj)-1]
		if filename[0] in agonist_bound:
			condition = get_condition(fulltraj).split(".h5")[0]
			if os.path.exists("%s/%s" %(features_dir, condition)):
				print("already featurized")
				trajs.append(fulltraj)
			else:
				trajs.append(fulltraj)
	
	featurize_partial = partial(featurize_known_traj, inactive_dir = inactive_dir, features_dir = features_dir)
	#pool = mp.Pool(mp.cpu_count()-1)
	#pool.map(featurize_partial, trajs)
	#pool.terminate()
	featurize_partial(trajs[0])

	print("Completed featurizing")

def compute_pnas_coords_and_distance(traj_file, inactive, active, scale = 7.14, residues_map = None):
	print "featurizing %s" %traj_file
	traj = md.load(traj_file)
	inactive_tuple = np.array([helix6_helix3_dist(inactive) / scale, rmsd_npxxy(inactive, inactive)])
	active_tuple = np.array([helix6_helix3_dist(active) / scale, rmsd_npxxy(active, inactive)])
	traj_coords = [helix6_helix3_dist(traj, residues_map) / scale, rmsd_npxxy(traj, inactive, residues_map), rmsd_npxxy(traj, active, residues_map), rmsd_connector(traj, inactive, residues_map), rmsd_connector(traj, active, residues_map)]
	traj_coords = np.transpose(np.vstack(traj_coords))
	active_vectors = traj_coords[:,[0,1]] - np.transpose(active_tuple)
	inactive_vectors = traj_coords[:,[0,1]] - np.transpose(inactive_tuple)

	inactive_distances = np.linalg.norm(inactive_vectors, axis = 1)
	active_distances = np.linalg.norm(active_vectors, axis = 1)
	distances = [inactive_distances, active_distances]
	#print distances[1]
	return [traj_coords, distances]

def featurize_pnas_distance_traj(traj_dir, ianctive, active, features_dir):
	#pnas_distances = 
	return

def convert_np_to_map(data):
	data_map = {}
	for i in range(0, len(data)):
		traj_data = data[i]
		for j in range(0, np.shape(traj_data)[0]):
			try:
				data_map["traj%d_frame%d" %(i,j)] = traj_data[j,:]
			except:
				data_map["traj%d_frame%d" %(i,j)] = [traj_data[j]]
	return data_map

def featurize_pnas_distance(traj_dir, features_dir, ext, inactive_dir, active_dir, inactive_distances_dir, active_distances_dir, coords_dir, inactive_distances_csv, active_distances_csv, coords_csv, scale = 7.14, residues_map = None):
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	inactive = md.load(inactive_dir)
	active = md.load(active_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	trajs = get_trajectory_files(traj_dir, ext = ext)
	#trajs = [t for t in trajs if "clone0.lh5" in t]
	#traj_objs = md.load(trajs)
	featurize_partial = partial(compute_pnas_coords_and_distance, inactive = inactive, active = active, scale = scale, residues_map = residues_map)
	pool = mp.Pool(16)
	features = pool.map(featurize_partial, trajs)
	#for traj in trajs:
	#cxx	featurize_partial(traj)
	pool.terminate()
	

	coords = [f[0] for f in features]
	inactive_distances = [f[1][0] for f in features]
	active_distances = [f[1][1] for f in features]

	verbosedump(coords, coords_dir)
	verbosedump(inactive_distances, inactive_distances_dir)
	verbosedump(active_distances, active_distances_dir)

	write_map_to_csv(coords_csv, convert_np_to_map(coords), ["frame", "tm3_tm6_dist", "rmsd_npxxy_inactive", "rmsd_npxxy_active", "rmsd_connector_inactive", "rmsd_connector_active"])
	write_map_to_csv(active_distances_csv, convert_np_to_map(active_distances), ["frame", "pnas_distance_active"])
	print("Completed featurizing")

def load_pdb_traj(pdb_file):
	print pdb_file
	return md.load_frame(pdb_file, index = 0)

def featurize_pnas_distance_pdbs(traj_dir, new_filename, features_dir, inactive_dir, active_dir, inactive_distances_dir, active_distances_dir, coords_dir, scale = 7.14):
	#if not os.path.exists(features_dir): os.makedirs(features_dir)

	inactive = md.load(inactive_dir)
	active = md.load(active_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	samples = get_trajectory_files(traj_dir, ext = ".pdb")[0:1]
	pool = mp.Pool(mp.cpu_count())
	trajs = pool.map(load_pdb_traj, samples)
	trajs_joined = trajs[0].join(trajs[1:])

	trajs_joined.save_hdf5(new_filename)

	features = compute_pnas_coords_and_distance(new_filename, inactive, active)

	coords = [f[0] for f in features]
	inactive_distances = [f[1][0] for f in features]
	active_distances = [f[1][1] for f in features]

	verbosedump(coords, coords_dir)
	verbosedump(inactive_distances, inactive_distances_dir)
	verbosedump(active_distances, active_distances_dir)

	print("Completed featurizing")
	

def load_features(filename):
	return np.transpose(verboseload(filename))

def featurize_sasa_individual(traj_file, bp_residues, anton, stride):
	print(traj_file)
	traj = md.load(traj_file, stride = stride)
	top = traj.topology
	protein_atoms = [a.index for a in top.atoms if a.residue.is_protein]
	traj = traj.atom_slice(protein_atoms)
	
	print(traj)

	#try:
	traj = fix_traj(traj)
	top = traj.topology
	sasa = md.shrake_rupley(traj, mode = 'atom')
	#except:
	#	return []

	bp_atoms = [a.index for a in top.atoms if a.residue.resSeq in bp_residues]
	print(np.shape(sasa))
	sasa = sasa[:, bp_atoms]
	total_sasa = sasa.sum(axis=1)
	return(total_sasa)

def featurize_sasa(traj_dir, traj_ext, bp_residues, sasa_file, residues_map = None, anton = False, skip = 1, stride = 0):
	trajs = get_trajectory_files(traj_dir, traj_ext)

	
	print("therer are initially %d bp residues" %len(bp_residues))
	if residues_map is not None:
		bp_residues = map_residues(residues_map, bp_residues)
	print("therer are now %d bp residues" %len(bp_residues))

	featurize_partial = partial(featurize_sasa_individual, bp_residues = bp_residues, anton = anton, stride = stride)
	pool = mp.Pool(mp.cpu_count())
	sasa = pool.map(featurize_partial, [trajs[i] for i in range(0,len(trajs),skip)])
	pool.terminate()
	#sasa = []
	#for traj in trajs[0:5]:
	#	sasa.append(featurize_partial(traj))
	sasa = np.concatenate(sasa)
	#print("SASA is:")
	#print(sasa)
	sample_names = [sample.split(".")[0] for sample in trajs]
	sasa_map = {}
	for i in range(0, len(sasa)):
		sasa_map[sample_names[i]] = [sasa[i]]
	write_map_to_csv(sasa_file, sasa_map, ["sample", "sasa"])
	#print(sasa)
	np.savetxt(sasa_file, sasa, delimiter=",")
	#plt.hist(sasa, bins=50)
	

	print("Completed featurizing")


#Compute standardized values for each feature 
def standardize_features(features_dir, features_ext, standardized_features_dir):
	if not os.path.exists(standardized_features_dir): os.makedirs(standardized_features_dir)
	feature_files = get_trajectory_files(features_dir, features_ext)
	features = load_file_list(feature_files)
	concatenated_features = np.concatenate(features)
	means = np.mean(concatenated_features, axis = 0)
	stdevs = np.std(concatenated_features, axis = 0)
	standardized_features = []
	for X in features: 
		X -= means
		X /= stdevs 
		standardized_features.append(X)

	print("Finished standardizing features")
	for i in range(0, len(feature_files)):
		filename = feature_files[i].split("/")[len(feature_files[i].split("/"))-1]
		new_filename = "%s/%s" %(standardized_features_dir, filename)
		verbosedump(standardized_features[i], new_filename)
	print("Finished saving all standardized features")
	return 

def sample_tIC_extremes(tica_projected_coords_dir, features_dir, standardized_features_dir, tica_extremes_dir, ext, percentile):
	if not os.path.exists(tica_extremes_dir): os.makedirs(tica_extremes_dir)

	tica_coords = verboseload(tica_projected_coords_dir)
	tica_concatenated = np.concatenate(tica_coords)
	print("np shape of tica_concatenated is")
	print(np.shape(tica_concatenated))
	del(tica_coords)

	standardized_features = load_file_list(get_trajectory_files(standardized_features_dir, ext))
	standardized_features_concatenated = np.concatenate(standardized_features)
	print("np shape of standardized_features_concatenated is ")
	print(np.shape(standardized_features_concatenated))
	del(standardized_features)

	features = load_file_list(get_trajectory_files(features_dir, ext))
	features_concatenated = np.concatenate(features)
	del(features)
	print("np shape of features_concatenated is ")
	print(np.shape(features_concatenated))

	print("Loaded standardized features, finding extrema now")

	for j in range(0,np.shape(tica_concatenated)[1]):
		print("Analyzing tIC%d" %j)
		tIC = tica_concatenated[:,j]
		print(np.shape(tIC))
		low_percentile = np.percentile(tIC, percentile)
		high_percentile = np.percentile(tIC, (100. - percentile))
		low_indices = np.where(tIC < low_percentile)[0]
		high_indices = np.where(tIC > high_percentile)[0]
		low_features = standardized_features_concatenated[low_indices,:]
		print(low_features)
		print(np.shape(low_features))
		high_features = standardized_features_concatenated[high_indices,:]
		np.savetxt("%s/tIC.%d_standardized_low_values.csv" %(tica_extremes_dir, (j+1)), low_features,  delimiter=",")
		np.savetxt("%s/tIC.%d_standardized_high_values.csv" %(tica_extremes_dir, (j+1)),high_features,  delimiter=",")
		
		low_features = features_concatenated[low_indices,:]
		high_features = features_concatenated[high_indices,:]
		np.savetxt("%s/tIC.%d_low_values.csv" %(tica_extremes_dir, (j+1)),low_features,  delimiter=",")
		np.savetxt("%s/tIC.%d_high_values.csv" %(tica_extremes_dir, (j+1)) ,high_features,  delimiter=",")

	print("Finished finding extreme values")



#Saves coefficients in each tIC for each feature (i.e., pair of residues)
def find_most_important_residues_in_tIC(traj_file, tica_object, tic_features_csv, contact_residues,tic_residue_csv, feature_coefs_csv, duplicated_feature_coefs_csv, cutoff):
	try:
		tica = verboseload(tica_object)
	except:
		tica = load_dataset(tica_object)
	print traj_file
	traj = md.load_frame(traj_file, 0)
	#traj = fix_traj(traj)
	top = traj.topology 
	#residue_pairs = compute_contacts_below_cutoff([traj_file, [0]], cutoff = cutoff, contact_residues = contact_residues, anton = True)
	residue_pairs = generate_features(tic_features_csv)
	new_residue_pairs = []
	for pair in residue_pairs:
		new_residue_pairs.append(("%s%d.%d" %(pair[0][2], pair[0][1], pair[0][0])), ("%s%d.%d" %(pair[1][2], pair[1][1], pair[1][0])))
	residue_pairs = new_residue_pairs
	#print traj_file

	
	top_indices_per_tIC = {}
	feature_coefs_per_tIC = {}
	duplicated_feature_coefs_per_tIC = {}


	#for each tIC:
		#for each feature, get the absolute component value
		#add to feature_coefs_per_tIC dictionary the absolute coefficient for that tIC
		#duplicate them for the analysis where we look at residues individually
		#sort by absolute coefficient value

	#for each tIC:
		#

	for i in range(0, np.shape(tica.components_)[0]):
		print i
		index_components = [(j,abs(tica.components_[i][j])) for j in range(0,np.shape(tica.components_)[1])]
		feature_coefs_per_tIC[i] = [component[1] for component in index_components]
		duplicated_feature_coefs_per_tIC[i] = [j for k in feature_coefs_per_tIC[i] for j in (k, k)] 
		index_components = sorted(index_components, key= lambda x: x[1],reverse=True)
		print(index_components[0:10])
		list_i = [index_components[j][0] for j in range(0,len(index_components))]
		top_indices_per_tIC[i] = list_i
	
	top_residues_per_tIC = {}
	for i in range(0, np.shape(tica.components_)[0]):
		top_residues_per_tIC[i] = []
		for index in top_indices_per_tIC[i]:
			residues = residue_pairs[index]
			top_residues_per_tIC[i].append(residues)
		top_residues_per_tIC[i] = [item for sublist in top_residues_per_tIC[i] for item in sublist]

	residue_list = residue_pairs

	feature_coefs_per_tIC["residues_0"] = [pair[0] for pair in residue_list]
	feature_coefs_per_tIC["residues_1"] = [pair[1] for pair in residue_list]
	duplicated_feature_coefs_per_tIC["residues"] = [residue for residue_pair in residue_list for residue in residue_pair]

	write_map_to_csv(tic_residue_csv, top_residues_per_tIC, [])
	write_map_to_csv(feature_coefs_csv, feature_coefs_per_tIC, [])
	write_map_to_csv(duplicated_feature_coefs_csv, duplicated_feature_coefs_per_tIC, [])
	return
	#print(top_residues_per_tIC)

def save_features_to_residues_map(traj_file, contact_residues, feature_residues_csv, cutoff, residues_map = None, exacycle = False):
	if residues_map is not None:
		contact_residues = [r for r in contact_residues if r in residues_map.keys()]
		if exacycle: contact_residues = [residues_map[key] for key in contact_residues]

	traj = md.load_frame(traj_file, 0)
	#traj = fix_traj(traj)
	top = traj.topology 
	residue_pairs, residue_infos = compute_contacts_below_cutoff([traj_file, [0]], cutoff = cutoff, contact_residues = contact_residues, anton = False)
	if exacycle:
		reverse_residues_map = {v: k for k, v in residues_map.items()}
		new_residue_pairs = []
		for residue_pair in residue_pairs:
			new_residue_pair = [reverse_residues_map[residue_pair[0]], reverse_residues_map[residue_pair[1]]]
			new_residue_pairs.append(new_residue_pair)
		residue_pairs = new_residue_pairs

		new_reisdue_infos = []
		for residue_info in residue_infos:
			new_residue_info = [(reverse_residues_map[residue_info[0][0]], residue_info[0][1], residue_info[0][2]), (reverse_residues_map[residue_info[1][0]], residue_info[1][1], residue_info[1][2])]
			new_residue_infos.append(new_residue_info)
		residue_infos = new_reisdue_infos

	print("There are: %d residue pairs" %len(residue_pairs))
	f = open(feature_residues_csv, "wb")
	f.write("feature, residue.1.resSeq, residue.1.res, residue.1.chain, residue.2.resSeq, residue.2.res, residue.2.chain,\n")
	for i in range(0, len(residue_infos)):
		f.write("%d, %d, %d, %d, %d, %d, %d,\n" %(i, residue_infos[i][0][0], residue_infos[i][0][1], residue_infos[i][0][2], residue_infos[i][1][0], residue_infos[i][1][1], residue_infos[i][1][2]))
	f.close()
	return 





'''
def compute_contacts_below_cutoff(traj_file_frames, cutoff = 100000.0, contact_residues = [], anton = True):
	traj_file = traj_file_frames[0]
	print(traj_file)
	frames = sorted(traj_file_frames[1])
	anton = True
	traj = md.load_frame(traj_file, index = 0)
	#print(traj)
	#traj = traj.slice(frames)
	#print(traj)
	#traj = md.open(traj_file)
	#print(frames)
	#for frame in frames:
	#	traj.append(traj.read(frame))
	#print(traj)
	#return
	#traj.seek(offset = 0, whence=2)
	#n_frames = traj.tell()
	#print(n_frames)
	#random_frame = random.randint(0,n_frames)

	#frame = md.load_frame(traj_file, index = random_frame)
	#print(frame)
	
	if anton:
		print("Fixing traj")
		fixed_traj = fix_traj(traj)
		fixed_top = fixed_traj.topology
	else: 
		fixed_traj = traj
		fixed_top = traj.topology
	distance_residues = []
	res_objects = [r for r in fixed_top.residues]
	for r in contact_residues:
		for res in res_objects:
			if res.resSeq == r and len(res._atoms) > 5:
				#print res._atoms
				distance_residues.append(res.index)
			elif res.resSeq == r and len(res._atoms) > 0 and len(res._atoms) <= 5:
				print("bad residue %d has %d atoms " %(r, len(res._atoms)))

	if len(contact_residues) != len(distance_residues):
		print "Residues are missing"
		print len(contact_residues)
		print len(distance_residues)
		#sys.exit()
		#return None
	
	combinations = itertools.combinations(distance_residues, 2)
	pairs = [c for c in combinations]
	print("Starting num pairs = %d" %len(pairs))
	#print pairs

	final_pairs = []

	distances = md.compute_contacts(fixed_traj, contacts = pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
	#print(distances)
	print(np.shape(distances))
	for j in range(0, np.shape(distances)[1]):
		distance_column = distances[:,j]
		#print(distance_column)
		if min(distance_column) < cutoff:
			final_pairs.append(pairs[j])
		#if (j % 1000 == 0): print j

	#print(final_pairs)
	return(final_pairs)
'''

'''
def featurize_custom_anton(traj_dir, features_dir, traj_ext, dihedral_residues = None, dihedral_types = None, contact_residues = None, agonist_bound = False, residues_map = None, contact_cutoff = None, n_random_samples = 0, clusters_map_file = None, n_clusters=1000, contact_csv=""):
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	all_trajs = get_trajectory_files(traj_dir, traj_ext)[0:1]
	trajs = []
	for fulltraj in all_trajs:
		#if "clone0.lh5" not in fulltraj: continue
		traj = fulltraj.split("/")
		filename = traj[len(traj)-1]
		#if agonist_bound is not False and filename[0] not in agonist_bound: continue
		filename_noext = filename.split(".")[0]
		if os.path.exists("%s/%s.h5.h5" %(features_dir, filename_noext)):
			print("already featurized")	
		else:
			trajs.append(fulltraj)

	if residues_map is not None:
		contact_residues = [r for r in contact_residues if r in residues_map.keys()]
		dihedral_residues = [d for d in dihedral_residues if d in residues_map.keys()]		

	print contact_residues	
	rand_trajs = []

	
	if os.path.exists(contact_csv):
		with open(contact_csv, 'rb') as f:
			reader = csv.reader(f)
			contact_list = map(tuple, reader)
	else:
		if n_random_samples > 0:
			for i in range(0,n_random_samples):
				traj_id = random.randint(0,len(all_trajs))
				rand_trajs.append(traj_id)

			compute_contacts_partial = partial(compute_contacts_below_cutoff, cutoff = contact_cutoff, contact_residues = contact_residues, anton = True)

			trajs = [all_trajs[i] for i in rand_trajs]
			final_contact_residue_pairs = []
			#for traj in trajs:
			#	final_contact_residue_pairs.append(compute_contacts_partial(traj))
			pool = mp.Pool(mp.cpu_count())
			print("Calculating all potential contact pairs")
			final_contact_residue_pairs = pool.map(compute_contacts_partial, trajs)
			pool.terminate()
			print("Finished finding all potential contact pairs")

		elif os.path.exists(clusters_map_file):
			f  = open(clusters_map_file, 'r')
			clusters_map = json.load(f)	
			clusters_map = {int(k):v for k,v in clusters_map.items()}
			traj_frame_map = {}
			for i in range(0,n_clusters):
				traj_frame_map[i] = []

			for key, cluster_list in clusters_map.iteritems():
				for traj_frame_tuple in cluster_list:
					traj = traj_frame_tuple[0]
					frame = traj_frame_tuple[1]
					traj_frame_map[traj].append(frame)
			traj_file_frames = [[all_trajs[i],traj_frame_map[i]] for i in range(0,len(all_trajs))]

			compute_contacts_partial = partial(compute_contacts_below_cutoff, cutoff = contact_cutoff, contact_residues = contact_residues, anton = True)
			pool = mp.Pool(mp.cpu_count())
			print("Calculating all potential contact pairs")
			final_contact_residue_pairs = []
			#for traj_file_frame in traj_file_frames:
			#	final_contact_residue_pairs.append(compute_contacts_partial(traj_file_frame))
			final_contact_residue_pairs = pool.map(compute_contacts_partial, traj_file_frames)
			pool.terminate()
			print("Finished finding all potential contact pairs")
		else:
			print("featurizing based on first frame of first trajectory")
			final_contact_residue_pairs = compute_contacts_below_cutoff([all_trajs[0], [0]], cutoff = contact_cutoff, contact_residues = contact_residues, anton = True)
			print("finished finding residue pairs")
			print(len(final_contact_residue_pairs))

		contact_set = set()
		for pair_list in final_contact_residue_pairs:
			for pair in pair_list:
				contact_set.add(pair)
		contact_list = list(contact_set)
		print("Final contact residue pairs for featurization:")
		#print(contact_list)
		print(len(contact_list))
		f = open(contact_csv, 'w')
		f.write('\n'.join('%s %s' % x for x in contact_list))
		f.close()
		print("written list to file, now featurizing with those contacts")
		

	#contact_residue_pairs = compute_contacts_below_cutoff(all_trajs[0], cutoff = contact_cutoff, contact_residues = contact_residues, anton = False)
	
	featurize_partial = partial(read_and_featurize_iter, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = dihedral_types, contact_residue_pairs = contact_list)
	if len(contact_list) > 4000:
		for traj in trajs:
			featurize_partial(traj)
	else:
		pool = mp.Pool(mp.cpu_count()/4)
		featurize_partial = partial(read_and_featurize_iter, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = dihedral_types, contact_residue_pairs = contact_list)
		pool.map(featurize_partial, trajs)
		pool.terminate()
	#for traj in trajs:
	#	print traj
	#	featurize_partial(traj)


	print("Completed featurizing")
'''


'''
def read_and_featurize_custom(traj_file, features_dir = None, condition=None, dihedral_types = ["phi", "psi", "chi1", "chi2"], dihedral_residues = None, contact_residue_pairs = []):
	#if "clone4.lh5" not in traj_file: return
	#top = md.load_frame(traj_file,index = 0).topology
	#atom_indices = [a.index for a in top.atoms if a.residue.resSeq != 130]
	#atom_indices = [a.index for a in top.atoms]
	print traj_file
	try:
		traj = md.load(traj_file)
	except:
		print "Could not load %s" %traj_file
		sys.exit()
		return
	top = traj.topology
'''
'''
	a = time.time()
	featurizer = DihedralFeaturizer(types = ['phi', 'psi', 'chi1', 'chi2'])
	features = featurizer.transform(traj)
	b = time.time()
	#print(b-a)
	print("original features has dim")
	print(np.shape(features))
'''
'''		
	a = time.time()

	if len(dihedral_residues) > 0:
		print "Analyzing dihedrals"
		dihedral_indices = []

		for dihedral_type in dihedral_types:
			if dihedral_type == "phi": dihedral_indices.append(phi_indices(top, dihedral_residues))
			if dihedral_type == "psi": dihedral_indices.append(psi_indices(top, dihedral_residues))
			if dihedral_type == "chi1": dihedral_indices.append(chi1_indices(top, dihedral_residues))
			if dihedral_type == "chi2": dihedral_indices.append(chi2_indices(top, dihedral_residues))

		#print("new features has dim %d" %(2*len(phi_tuples) + 2*len(psi_tuples) + 2*len(chi2_tuples)))

		#print("feauturizing manually:")
		dihedral_angles = []

		for dihedral_type in dihedral_indices:
			angles = np.transpose(ManualDihedral.compute_dihedrals(traj=traj,indices=dihedral_type))
			dihedral_angles.append(np.sin(angles))
			dihedral_angles.append(np.cos(angles))

		manual_features = np.transpose(np.concatenate(dihedral_angles))
		print("Dihedral features for %s has shape: " %traj_file)
		print(np.shape(manual_features)) 

	if len(contact_residue_pairs) > 0:
		contact_features = md.compute_contacts(traj, contacts = contact_residue_pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
		#print contact_features
		#print(np.shape(contact_features))
		if len(dihedral_residues) > 0:
			manual_features = np.column_stack((manual_features, contact_features))
		else:
			manual_features = contact_features

	b = time.time()

	print("new features %s has shape: " %traj_file)
	print(np.shape(manual_features))

	traj_lastname = traj_file.split("/")[len(traj_file.split("/"))-1]
	traj_noext = traj_lastname.split(".")[0]
	features_file =  "%s/%s.h5" %(features_dir, traj_noext)

	verbosedump(manual_features, features_file)

def read_and_featurize_custom_anton(traj_file, features_dir = None, condition=None, dihedral_types = ["phi", "psi", "chi1", "chi2"], dihedral_residues = None, contact_residues = None):
	#if "23" not in traj_file and "24" not in traj_file: return
	top = md.load_frame(traj_file,index = 0).topology
	#atom_indices = [a.index for a in top.atoms if a.residue.resSeq != 130]
	atom_indices = [a.index for a in top.atoms]
	traj = md.load(traj_file, atom_indices=atom_indices)
	print traj_file
	#print traj
	#print("loaded trajectory")

'''
'''
	a = time.time()
	featurizer = DihedralFeaturizer(types = ['phi', 'psi', 'chi2'])
	features = featurizer.transform(traj)
	b = time.time()
	#print(b-a)
	print("original features has dim")
	print(np.shape(features))
'''
'''
	a = time.time()
	dihedral_indices = []
	residue_order = []
	if len(dihedral_residues) > 0:
		for dihedral_type in dihedral_types:
			if dihedral_type == "phi": dihedral_indices.append(phi_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "psi": dihedral_indices.append(psi_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "chi1": dihedral_indices.append(chi1_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "chi2": dihedral_indices.append(chi2_indices(fix_topology(top), dihedral_residues))

		#print("new features has dim %d" %(2*len(phi_tuples) + 2*len(psi_tuples) + 2*len(chi2_tuples)))

		#print("feauturizing manually:")
		dihedral_angles = []

		for dihedral_type in dihedral_indices:
			angles = np.transpose(ManualDihedral.compute_dihedrals(traj=traj,indices=dihedral_type))
			dihedral_angles.append(np.sin(angles))
			dihedral_angles.append(np.cos(angles))

		manual_features = np.transpose(np.concatenate(dihedral_angles))

	if len(contact_residue_pairs) > 0:
		fixed_traj = fix_traj(traj)
		fixed_top = fixed_traj.topology
		contact_features = md.compute_contacts(fixed_traj, contacts = contact_residues, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
		if len(dihedral_residues) > 0: 
			manual_features = np.column_stack((manual_features, contact_features))
		else:
			manual_features = contact_features


	b = time.time()

	print("new features %s has shape: " %traj_file)
	print(np.shape(manual_features))

	if condition is None:
		condition = get_condition(traj_file)

	verbosedump(manual_features, "%s/%s.h5" %(features_dir, condition))
'''

'''
def featurize_custom(traj_dir, features_dir, traj_ext, dihedral_residues = None, dihedral_types = None, contact_residues = None, agonist_bound = False, residues_map = None, contact_cutoff = None, n_random_samples = 1, contact_csv=""):
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	all_trajs = get_trajectory_files(traj_dir, traj_ext)[0:1]
	trajs = []
	for fulltraj in all_trajs:
		#if "clone0.lh5" not in fulltraj: continue
		traj = fulltraj.split("/")
		filename = traj[len(traj)-1]
		#if agonist_bound is not False and filename[0] not in agonist_bound: continue
		filename_noext = filename.split(".")[0]
		if os.path.exists("%s/%s.h5" %(features_dir, filename_noext)):
			print("already featurized")	
		else:
			trajs.append(fulltraj)

	pool = mp.Pool(mp.cpu_count())

	if residues_map is not None:
		dihedral_residues = map_residues(residues_map, dihedral_residues)
		contact_residues = map_residues(residues_map, contact_residues)

	print contact_residues	
	rand_trajs = []

	for i in range(0,n_random_samples):
		traj_id = random.randint(0,len(all_trajs))
		rand_trajs.append(traj_id)

	compute_contacts_partial = partial(compute_contacts_below_cutoff, cutoff = contact_cutoff, contact_residues = contact_residues, anton = False)
	pool = mp.Pool(mp.cpu_count())
	trajs = [all_trajs[i] for i in rand_trajs]
	final_contact_residue_pairs = []
	for traj in trajs:
		final_contact_residue_pairs.append(compute_contacts_partial(traj))
	#final_contact_residue_pairs = pool.map(compute_contacts_partial, )
	pool.terminate()

	contact_set = set()
	for pair in final_contact_residue_pairs:
		contact_set.add(pair)
	contact_list = list(contact_set)
	#print("Final contact residue pairs for featurization:")
	#print(contact_list)
	f = open(contact_csv, 'w')
	f.write("\n".join(contact_list))
	f.close()
	

	#contact_residue_pairs = compute_contacts_below_cutoff(all_trajs[0], cutoff = contact_cutoff, contact_residues = contact_residues, anton = False)

	featurize_partial = partial(read_and_featurize_custom, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = dihedral_types, contact_residue_pairs = contact_list)
	pool.map(featurize_partial, trajs)
	pool.terminate()
	#for traj in trajs:
	#	featurize_partial(traj)


	print("Completed featurizing")
'''

