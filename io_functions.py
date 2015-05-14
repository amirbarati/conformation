import fileinput

def write_map_to_csv(filename, data_map, firstline):
	csvfile = open(filename, "wb")
	csvfile.write(firstline)
	for key in sorted(data_map.keys()):
		csvfile.write("%s, " %key)
		for value in data_map[key]:
			csvfile.write("%f, " %value)
		csvfile.write("\n")
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

def convert_csv_to_map_nocombine(filename):
	print filename
	rmsds = open(filename, "rb")
	lines = rmsds.readlines()
	if "docking" in filename:
		print lines
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

		for i in range(1,len(line)):
			rmsd = float(line[i].split('\\')[0])
			#print rmsd
			#if rmsd > 3.0: print "%s %f" %(cluster, rmsd)

			if cluster in rmsd_map.keys():
				rmsd_map[cluster].append(rmsd)
			else:
				rmsd_map[cluster] = [rmsd]
	return rmsd_map

def convert_csv_to_map(filename):
	rmsds = open(filename, "rb")
	lines = rmsds.readlines()
	rmsd_map = {}
	i = 0
	for line in lines:
		print line
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
	return rmsd_map


def get_condition(filename):
	filename = filename.split('/')[len(filename.split('/'))-1]
	pieces = filename.split('-')
	condition = "%s-%s" %(pieces[0], pieces[1])
	return condition 

def read_trajectory(directory, filename, stride=10):
	print("reading %s" %(filename))
	traj = md.load(filename, stride=stride, top="/home/harrigan/compute/wetmsm/gpcr/des/system_mae_to_pdb/des_trajs/DESRES-Trajectory_pnas2011b-H-05-all/system.pdb")
	return traj