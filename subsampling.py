import mdtraj as md
from io_functions import *
import subprocess
import multiprocessing as mp

def reimage(traj_file):
	print("Examining %s" %traj_file)
	traj = md.load(traj_file)
	traj = fix_traj(traj)
	top = traj.topology
	traj_name = traj_file.spit(".h5")[0]
	directory = traj_file.split("/")[len(traj_file.split("/"))-2]
	traj.save_dcd("%s.dcd" %traj_name)
	top.save_pdb("%s.pdb" %traj_name)
	del(traj)
	del(top)
	f = open("%s.in" %traj_name, 'w')
	f.write("parm %s" %("%s.pdb" %traj_name))
	f.write("trajin %s" %("%s.dcd" %traj_name))
	f.write("autoimage")
	f.write("trajout %s" %("%s_reimaged.dcd" %traj_name))
	f.write("trajout %s onlyframes 1" %("%s_reimaged.pdb" %traj_name))
	f.write("go")
	f.close()
	subprocess.call("ml load amber/14-intel")
	os.chdir(directory)
	subprocess.call("cpptraj %s.in" %traj_name)
	traj = md.load("%s_reimaged.dcd" %traj_name, top = "%s_reimaged.pdb" %traj_name)
	traj.save("%s_reimaged.h5" %traj_name)


def reimage_amber(traj_dir):
	trajs = get_trajectory_files(traj_dir, ext = ".h5")
	pool = mp.Pool(mp.cpu_count())
	pool.map(reimage, trajs)
	pool.close()
	#map(reimage, trajs)

def subsample_traj(traj, stride=5, top=None):
	directory = traj.split("/")
	simulation = directory[len(directory)-2]
	dcd_file = directory[len(directory)-1]
	condition = "%s-%s" %(simulation.split('-')[1], simulation.split('-')[2])
	print("analyzing simulation %s file %s" %(simulation, dcd_file))
	top_file = top

	top = md.load_frame(traj, 0, top=top_file).topology
	atom_indices = [a.index for a in top.atoms if str(a.residue)[0:3] != "POP" and not a.residue.is_water and str(a.residue)[0:2] != "NA" and str(a.residue)[0:2] != "CL"]

	traj = md.load(traj, stride=stride, top=top_file, atom_indices=atom_indices)
	print("traj loaded")

	new_file = "%s_stride%d.h5" %(dcd_file.rsplit( ".", 1 )[ 0 ] , stride)
	new_root_dir = "/scratch/users/enf/b2ar_analysis/subsampled_allprot"

	

	new_condition_dir = "%s/%s" %(new_root_dir, condition)

	new_file_full = "%s/%s/%s" %(new_root_dir, condition, new_file)
	print("saving trajectory as %s" %new_file_full)
	traj.save(new_file_full)


def subsample(directory, stride=5):
	conditions = sorted(os.listdir(directory))
	for condition in conditions: 
		subdirectory = "%s/%s" %(directory, condition)
		traj_dir = condition.split("_")[1]
		agonist_residues = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

		if traj_dir.split('-')[1] not in agonist_residues:
			print("skipping simulation")
			continue
		#	if traj_dir.split('-')[2] != "00": continue

		new_root_dir = "/scratch/users/enf/b2ar_analysis/subsampled_allprot"
		sim_dir = "%s-%s" %(traj_dir.split('-')[1], traj_dir.split('-')[2])
		new_condition_dir = "%s/%s" %(new_root_dir, sim_dir)
		if os.path.exists(new_condition_dir):
			print("We have already subsampled_allprot this simulation")
			continue
		else:
			os.makedirs(new_condition_dir)

		top_file = "%s/system.pdb" %subdirectory
		top_file = alias_pdb(top_file,sim_dir)
		

		traj_dir = "%s/%s" %(subdirectory, traj_dir)
		traj_files = get_trajectory_files(traj_dir)

		subsample_parallel = partial(subsample_traj, top = top_file, stride=stride)

		num_workers = mp.cpu_count()
		pool = mp.Pool(num_workers)
		pool.map(subsample_parallel, traj_files)
		pool.terminate()
		gc.collect()
		#subsample_parallel(traj_files[0])

		subsampled_allprot_trajs = get_trajectory_files("/scratch/users/enf/b2ar_analysis/subsampled_allprot/%s" %sim_dir)

		combined_traj = md.load(subsampled_allprot_trajs)
		combined_traj.save("/scratch/users/enf/b2ar_analysis/subsampled_allprot_combined/%s.h5" %sim_dir)
		first_frame = md.load_frame("/scratch/users/enf/b2ar_analysis/subsampled_allprot_combined/%s.h5" %sim_dir, index=1)
		first_frame.save_pdb("/scratch/users/enf/b2ar_analysis/%s_firstframe.pdb" %condition)
		print("trajectory has been subsampled_allprot, combined, and saved")