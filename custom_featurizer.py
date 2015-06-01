import mdtraj as md
import numpy as np
from analysis import *
from msmbuilder.utils import verbosedump, verboseload

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
			print "No phs found for %s " %c.name

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

	top = fix_topology(top)
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

	top = fix_topology(top)
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
		if None not in dihedral:
			dihedral = tuple(dihedral)
			chi2_tuples.append(dihedral)
			#print residue.resSeq
		elif dihedral != [None, None, None, None] and str(residue.name)[0:3] in chi2_residues:
			print "no chi2 found for %s" %str(residue)		


	return chi2_tuples


def read_and_featurize_custom(traj_file, features_dir = None, condition=None, dihedral_types = ["phi", "psi", "chi1", "chi2"], dihedral_residues = None, contact_residues = None):
	#if "23" not in traj_file and "24" not in traj_file: return
	top = md.load_frame(traj_file,index = 0).topology
	#atom_indices = [a.index for a in top.atoms if a.residue.resSeq != 130]
	atom_indices = [a.index for a in top.atoms]
	traj = md.load(traj_file, atom_indices=atom_indices)
	print traj_file
	#print traj
	#print("loaded trajectory")

	'''
	a = time.time()
	featurizer = DihedralFeaturizer(types = ['phi', 'psi', 'chi2'])
	features = featurizer.transform(traj)
	b = time.time()
	#print(b-a)
	print("original features has dim")
	print(np.shape(features))
	'''
	a = time.time()
	dihedral_indices = []

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

	if contact_residues is not None:
		fixed_traj = fix_traj(traj)
		fixed_top = fixed_traj.topology
		distance_residues = []
		res_objects = [r for r in fixed_top.residues]
		for r in contact_residues:
			for res in res_objects:
				if res.resSeq == r and len(res._atoms) > 5:
					#print res._atoms
					distance_residues.append(res.index)
		if len(contact_residues) != len(distance_residues):
			print "SOMETHING WENT WRONG"
			print len(contact_residues)
			print len(distance_residues)
			sys.exit()
			return None
		
		combinations = itertools.combinations(distance_residues, 2)
		pairs = [c for c in combinations]
		#print pairs
		contact_features = md.compute_contacts(traj, contacts = pairs, scheme = 'closest-heavy', ignore_nonprotein=False)[0]
		#print contact_features
		#print(np.shape(contact_features))
		manual_features = np.column_stack((manual_features, contact_features))

	b = time.time()

	print("new features %s has shape: " %traj_file)
	print(np.shape(manual_features))

	if condition is None:
		condition = get_condition(traj_file)

	verbosedump(manual_features, "%s/%s.h5" %(features_dir, condition))


def featurize_custom(traj_dir, features_dir, dihedral_residues = None, dihedral_types = None, contact_residues = None):
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	#agonist_bound = ['D', 'J']
	all_trajs = get_trajectory_files(traj_dir)
	trajs = []
	for fulltraj in all_trajs:
		traj = fulltraj.split("/")
		filename = traj[len(traj)-1]
		if filename[0] in agonist_bound:
			condition = get_condition(fulltraj)
			if os.path.exists("%s/%s.h5" %(features_dir, condition)):
				print("already featurized")
			else:
				trajs.append(fulltraj)
	
	featurize_partial = partial(read_and_featurize_custom, features_dir = features_dir, dihedral_residues = dihedral_residues, dihedral_types = dihedral_types, contact_residues = contact_residues)
	pool = mp.Pool(mp.cpu_count())
	pool.map(featurize_partial, trajs)
	pool.terminate()
	#print trajs[0]
	#featurize_partial(trajs[0])
	#for traj in trajs:
	#	featurize_partial(traj)

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
			condition = get_condition(fulltraj)
			if os.path.exists("%s/%s.h5" %(features_dir, condition)):
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

def compute_pnas_coords_and_distance(traj_file, inactive, active, scale = 1.0):
	print "featuring %s" %traj_file
	traj = md.load(traj_file)
	inactive_tuple = np.array([helix6_helix3_dist(inactive) / scale, rmsd_npxxy(inactive, inactive)])
	active_tuple = np.array([helix6_helix3_dist(active) / scale, rmsd_npxxy(active, inactive)])
	traj_coords = [helix6_helix3_dist(traj) / scale, rmsd_npxxy(traj, inactive)]
	traj_coords = np.transpose(np.vstack(traj_coords))
	active_vectors = traj_coords - np.transpose(active_tuple)
	inactive_vectors = traj_coords - np.transpose(inactive_tuple)

	inactive_distances = np.linalg.norm(inactive_vectors, axis = 1)
	active_distances = np.linalg.norm(active_vectors, axis = 1)
	distances = [inactive_distances, active_distances]

	return [traj_coords, distances]

def featurize_pnas_distance_traj(traj_dir, ianctive, active, features_dir):
	#pnas_distances = 
	return

def featurize_pnas_distance(traj_dir, features_dir, ext, inactive_dir, active_dir, inactive_distances_dir, active_distances_dir, coords_dir, scale):
	if not os.path.exists(features_dir): os.makedirs(features_dir)

	inactive = md.load(inactive_dir)
	active = md.load(active_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	trajs = get_trajectory_files(traj_dir, ext = ext)
	featurize_partial = partial(compute_pnas_coords_and_distance, inactive = inactive, active = active, scale = scale)
	pool = mp.Pool(mp.cpu_count())
	features = pool.map(featurize_partial, trajs)
	pool.terminate()
	

	coords = [f[0] for f in features]
	inactive_distances = [f[1][0] for f in features]
	active_distances = [f[1][1] for f in features]

	verbosedump(coords, coords_dir)
	verbosedump(inactive_distances, inactive_distances_dir)
	verbosedump(active_distances, active_distances_dir)

	print("Completed featurizing")

def load_pdb_traj(pdb_file):
	print pdb_file
	return md.load_frame(pdb_file, index = 0)

def featurize_pnas_distance_pdbs(traj_dir, new_filename, features_dir, inactive_dir, active_dir, inactive_distances_dir, active_distances_dir, coords_dir):
	#if not os.path.exists(features_dir): os.makedirs(features_dir)

	inactive = md.load(inactive_dir)
	active = md.load(active_dir)

	agonist_bound = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
	samples = get_trajectory_files(traj_dir, ext = ".pdb")
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