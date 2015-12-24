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
import pickle

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

def convert_residue_pairs_to_mdtraj_indices(top, residue_pairs):
	resIndices = []
	for pair in residue_pairs:
		resSeq0 = pair[0][1]
		chain0 = pair[0][0]
		resSeq1 = pair[1][1]
		chain1 = pair[1][0]
		indices = [r.index for r in top.residues if pair.residue_i.is_mdtraj_res_equivalent(r) and not r.is_water]
		if len(indices) == 0:
			print("FATAL: No residues in trajectory for residue %d" %resSeq0)
			return None
		else:
			ind_i = indices[0]
			for j in indices:
				if j != ind_i: 
					#print("Warning: multiple res objects for residue %d " %resSeq0)
					if "CB" in [str(a) for a in r.atoms for r in top.residues if r.index == ind_i]:
						ind_i = j

		indices = [r.index for r in top.residues if pair.residue_j.is_mdtraj_res_equivalent(r) and not r.is_water]
		if len(indices) == 0:
			print("FATAL: No residues in trajectory for residue %d" %resSeq0)
			return None
		else:
			ind_j = indices[0]
			for j in indices:
				if j != ind_j: 
					#print("Warning: multiple res objects for residue %d " %resSeq0)
					if "CB" in [str(a) for a in r.atoms for r in top.residues if r.index == ind_j]:
						ind_j = j

		resIndices.append((ind_i, ind_j))
	print("looking at %d pairs for trajectory" %len(resIndices))
	return(resIndices)



def read_and_featurize(traj_file, features_dir = None, condition=None, 
											 dihedral_types = ["phi", "psi", "chi1", "chi2"], 
											 dihedral_residues = None, residue_pairs = [], 
											 iterative = True):

	dihedral_indices = []
	residue_order = []
	if len(dihedral_residues) > 0:
		for dihedral_type in dihedral_types:
			if dihedral_type == "phi": dihedral_indices.append(phi_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "psi": dihedral_indices.append(psi_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "chi1": dihedral_indices.append(chi1_indices(fix_topology(top), dihedral_residues))
			if dihedral_type == "chi2": dihedral_indices.append(chi2_indices(fix_topology(top), dihedral_residues))

		dihedral_angles = []

		for dihedral_type in dihedral_indices:
			angles = np.transpose(ManualDihedral.compute_dihedrals(traj=traj,indices=dihedral_type))
			dihedral_angles.append(np.sin(angles))
			dihedral_angles.append(np.cos(angles))

		manual_features = np.transpose(np.concatenate(dihedral_angles))

	if len(residue_pairs) > 0:
		top = md.load_frame(traj_file, index=0).topology
		residue_pairs = convert_residue_pairs_to_mdtraj_indices(top, resSeq_pairs)
		contact_features = []
		if iterative:
			try:
				for chunk in md.iterload(traj_file, chunk = 1000):
					chunk_features = md.compute_contacts(chunk, contacts = residue_pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
					contact_features.append(chunk_features)
				contact_features = np.concatenate(contact_features)
			except Exception,e:
				print str(e)
				print("Failed")
				return
		else:
			try:
				traj = md.load(traj_file)
				contact_features =  md.compute_contacts(traj, contacts = residue_pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
			except Exception,e:
				print str(e)
				print("Failed for traj")
				return
		if len(dihedral_residues) > 0: 
			manual_features = np.column_stack((manual_features, contact_features))
		else:
			manual_features = contact_features

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
	residue_to_mdtraj_index = {}
	residue_full_infos = []

	for i in range(0, len(contact_residues)):
		residue = contact_residues[i]
		indices = [r.index for r in top.residues if residue.is_mdtraj_res_equivalent(r)]
		if len(indices) == 0:
			print("No residues in trajectory for residue %d chain %d" %(residue.resSeq, residue.chain_id)
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
			residue_to_mdtraj_index[residue] = ind
	
	residue_combinations = itertools.combinations(distance_residues, 2)
	residue_pairs = [c for c in residue_combinations]
	residue_index_combinations = []

	for combination in residue_pairs:
		res0 = combination[0]
		res1 = combination[1]
		mdtraj_index0 = residue_to_mdtraj_index[res0]
		mdtraj_index1 = residue_to_mdtraj_index[res1]
		mdtraj_index_combinations.append((mdtraj_index0, mdtraj_index1))


	final_residue_pairs = []
	final_mdtraj_index_pairs = []

	distances = md.compute_contacts(frame, contacts = mdtraj_index_combinations, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
	#print(distances)
	print(np.shape(distances))
	for i in range(0, len(distances[0])):
		distance = distances[0][i]
		#print(distance)
		if distance < cutoff:
			final_mdtraj_index_pairs.append(mdtraj_index_combinations[i])
			final_residue_pairs.append(sorted(residue_pairs[i]))

	print("There are %d residue-residue contacts below cutoff in structure." %len(final_residue_pairs))
	
	return final_residue_pairs


def which_trajs_to_featurize(traj_dir, traj_ext, features_dir):
	all_trajs = get_trajectory_files(traj_dir, traj_ext)
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

def featurize_contacts_custom(traj_dir, features_dir, traj_ext, structures, 
															contact_residue_pairs_file = None, 
															dihedral_residues = None, 
															dihedral_types = None, 
															contact_residues = None, 
															agonist_bound = False, 
															residues_map = None, 
															contact_cutoff = None, 
															parallel = False, exacycle = False):
	'''
	Nb: The input to this function, either contact_residues or contact_residue_pairs_file, must contain instances 
	of object Residue(). The .resSeq attribute of each such instance must refer to residue numbering in your reference
	structure/PDB. This is to standardize it across multiple simulation conditions. The residues_map must be given to map 
	the ref structure/PDB residue ID numbers to the "resSeq" attributes that mdtraj requires.

	you can also input a residue_map, a dictionary that maps residue_object --> residue_object. The reason for this is that there
		is not a consensus residue numbering for the same protein.
	'''
	trajs = which_trajs_to_featurize(traj_dir, traj_ext, features_dir)

	try:
		contact_residue_pairs = generate_features(contact_residue_pairs_file)
		print("Features already computed")
		if exacycle: contact_residue_pairs = [residues_map[key] for key in contact_residue_pairs]
	if True: #except:
		contact_residue_pairs = []
		for structure in structures: 
			structure_contact_residue_pairs = compute_contacts_below_cutoff([structure,0], cutoff = contact_cutoff, contact_residues = contact_residues, anton = False)
			for pair in structure_contact_residue_pairs:
				if pair not in contact_residue_pairs: contact_residue_pairs.append(pair)
		contact_residue_pairs.sort(key=operator.itemgetter(1))
		print("There are %d pairs of residues to be used in contact featurization." % len(contact_residue_pairs))
		print("Saving contact feature residue pairs to disk.")
		with open(contact_residue_pairs_file, "wb") as f:
			pickle.dump(contact_residue_pairs, f)

	print("About to featurize trajectories based on the chosen featurization scheme.")
	featurize_partial = partial(read_and_featurize, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = dihedral_types, resSeq_pairs = contact_residue_pairs, iterative = True)
	if parallel:
		pool = mp.Pool(mp.cpu_count())
		pool.map(featurize_partial, trajs)
		pool.terminate()
	else:
		for traj in trajs:
			featurize_partial(traj)

	print("Completed featurizing")

def save_feature_residues_pkl(traj_dir, features_dir, traj_ext, structures, contact_residue_pairs_file = None, dihedral_residues = None, dihedral_types = None, contact_residues = None, agonist_bound = False, residues_map = None, contact_cutoff = None, parallel = False, exacycle = False):
	if residues_map is not None:
		contact_residues = [r for r in contact_residues if r in residues_map.keys()]
		if exacycle: contact_residues = [residues_map[key] for key in contact_residues]	

	if contact_residue_pairs_file == "" or (not os.path.exists(contact_residue_pairs_file)):
		contact_residue_pairs = []
		for structure in structures: 
			contact_residue_pairs.append(compute_contacts_below_cutoff([structure,0], cutoff = contact_cutoff, contact_residues = contact_residues, anton = False)[0])
		contact_residue_pairs = sorted(list(set(contact_residue_pairs)))
		with open(contact_residue_pairs_file, "wb") as f:
			pickle.dump(contact_residue_pairs, f)
	else:
		print("Features already computed")
		contact_residue_pairs = generate_features(contact_residue_pairs_file)
		if exacycle: contact_residue_pairs = [residues_map[key] for key in contact_residue_pairs]
	print(contact_residue_pairs)
	print("Number of contact pairs = %d" %len(contact_residue_pairs))

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

def compute_pnas_coords_and_distance(traj_file, inactive, active, scale = 7.14, residues_map = None, structure=None, connector_residues=[], npxxy_residues=[], tm6_tm3_residues=[]):
	print "featurizing %s" %traj_file
	if structure is not None:
		traj = md.load(traj_file, top=structure)
	else:
		traj = md.load(traj_file)
	inactive_tuple = np.array([helix6_helix3_dist(inactive, residues=tm6_tm3_residues) / scale, rmsd_npxxy(inactive, inactive, residues=npxxy_residues)])
	active_tuple = np.array([helix6_helix3_dist(active, residues=tm6_tm3_residues) / scale, rmsd_npxxy(active, inactive, residues=npxxy_residues)])
	traj_coords = [helix6_helix3_dist(traj, residues_map, residues=tm6_tm3_residues) / scale, rmsd_npxxy(traj, inactive, residues_map, residues=npxxy_residues), rmsd_npxxy(traj, active, residues_map, resiudes=npxxy_residues), rmsd_connector(traj, inactive, residues_map, residues=connector_residues), rmsd_connector(traj, active, residues_map, residues=connector_residues)]
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

def featurize_pnas_distance(traj_dir, features_dir, ext, inactive_dir, 
														active_dir, inactive_distances_dir, 
														active_distances_dir, coords_dir, 
														inactive_distances_csv, active_distances_csv, 
														coords_csv, scale = 7.14, residues_map = None, 
														structure=None, connector_residues=[], 
														npxxy_residues=[], tm6_tm3_residues=[]):
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	inactive = md.load(inactive_dir)
	active = md.load(active_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	trajs = get_trajectory_files(traj_dir, ext = ext)

	featurize_partial = partial(compute_pnas_coords_and_distance, inactive = inactive, 
															active = active, scale = scale, residues_map = residues_map, structure=structure, 
															connector_residues=[], npxxy_residues=[], tm6_tm3_residues=[])
	pool = mp.Pool(mp.cpu_count())
	features = pool.map(featurize_partial, trajs)
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
		try:
			new_residue_pairs.append(("%s%d.%d" %(pair[0][2], pair[0][1], pair[0][0])), ("%s%d.%d" %(pair[1][2], pair[1][1], pair[1][0])))
		except:
			new_residue_pairs.append(pair)
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

def save_features_to_residues_map(traj_file, residue_pairs, feature_residues_csv, cutoff, residues_map = None, exacycle = False):
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
	f.write("feature, residue.1.resSeq, residue.2.resSeq, residue.1.res, residue.1.chain, residue.2.res, residue.2.chain,\n")
	for i in range(0, len(residue_infos)):
		f.write("%d, %d, %d, %d, %d, %d, %d,\n" %(i, residue_infos[i][0][0], residue_infos[i][1][0], residue_infos[i][0][1], residue_infos[i][0][2], residue_infos[i][1][1], residue_infos[i][1][2]))
	f.close()
	return 

