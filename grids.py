from PDB_Order_Fixer import PDB_Order_Fixer
import mdtraj as md
import os
import numpy as np
import h5py
import datetime
import multiprocessing as mp
import copy
import gc
from functools import partial 
import time
import fileinput
import subprocess
from subprocess import Popen
import sys
from io_functions import *
#import pytraj.io as mdio
#from pytraj import adict

'''
If simulation was run under periodic boundary conditions, this will reimage the trajectory.
It takes as input the trajectory file (traj_file), the directory that that trajectory is in (traj_dir),
the directory to which you would like to save the new trajectory, and the extension (ext) of the file. 

In the docking pipeline, I do this to all PDB files containing receptors to which I would like to dock. 
You can skip this step if it already has been reimaged. 

It requires that Pytraj be installed. This can be very annoying to do. It was easy to install on 
Sherlock but not on Biox3. 

'''

'''
def reimage_traj(traj_file, traj_dir, save_dir, ext):
	if ext == ".pdb":
		file_lastname = traj_file.split("/")[len(traj_file.split("/"))-1]
		filename = file_lastname.split(".")[0]
		h5_filename = file_lastname
		new_h5_file = "%s/%s" %(save_dir, h5_filename)
		if os.path.exists(new_h5_file):
			print "already reimaged"
			return 

		traj_pytraj = mdio.load(traj_file, top = traj_file)[:]
		#traj_pytraj.fixatomorder()
		traj_pytraj.autoimage()

		
		traj_pytraj.save(new_h5_file)
		print "saving %s" %h5_filename

	else:
		traj_file_lastname = traj_file.split("/")[len(traj_file.split("/"))-1]
		traj_filename = traj_file_lastname.split(".")[0]
		traj_dcd = "%s/%s.dcd" %(traj_dir, traj_filename)
		traj_pdb = "%s/%s.pdb" %(traj_dir, traj_filename)
		traj = md.load(traj_file)
		traj_frame = md.load_frame(traj_file, index=0)
		traj.save_dcd(traj_dcd)
		traj_frame.save_pdb(traj_pdb)

		traj_pytraj = mdio.load(traj_dcd, top = traj_pdb)[:]
		traj_pytraj.autoimage()

		file_lastname = traj_file.split("/")[len(traj_file.split("/"))-1]
		filename = file_lastname.split(".")[0]
		dcd_filename = "%s_temp.dcd" %filename
		top_filename = "%s_temp.pdb" %filename
		h5_filename = file_lastname
		new_dcd_file = "%s/%s" %(save_dir, dcd_filename)
		new_top_file = "%s/%s" %(save_dir, top_filename)
		new_h5_file = "%s/%s" %(save_dir, h5_filename)
		print new_dcd_file
		print new_top_file
		traj_pytraj.save(new_dcd_file)
		traj_pytraj.save(new_top_file)

		new_traj = md.load(new_dcd_file, top = traj_pdb)
		new_traj.save(new_h5_file)
		os.remove(traj_dcd)
		os.remove(traj_pdb)
		os.remove(new_dcd_file)
		os.remove(new_top_file)
	return
'''

'''
If sim was run under periodic boundary conditions, this will reimage all trajectories in directory traj_dir
'''

def reimage_trajs(traj_dir, ext = ".pdb"):
	print "traj dir = %s" %traj_dir
	new_dir = "%s_reimaged" %traj_dir
	print "new dir = %s" %new_dir

	if not os.path.exists(new_dir): os.makedirs(new_dir)

	trajs = get_trajectory_files(traj_dir, ext = ext)

	reimage = partial(reimage_traj, save_dir = new_dir, traj_dir = traj_dir, ext = ext)

	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(reimage, trajs)
	pool.terminate()
	#reimage(trajs[0])
	#for traj in trajs:
	#	reimage(traj)
	return

'''
The following two functions take as input a directory containing PDB files containing structures to which you would like to dock.
It will prepare the protein with Schrodinger's tools (add hydrogens, SS bonds (no, not that SS!), bond orders, etc.) and then save
an .mae file, which is required for docking.
'''

def pprep_prot(pdb, ref):
	pdb_name = pdb.split("/")[len(pdb.split("/"))-1]
	new_pdb = pdb_name.rsplit( ".", 1 )[ 0 ]
	new_pdb = "%s.mae" %(new_pdb)
	if os.path.exists(new_pdb): 
		print "already prepped and mae'd protein"
		return
	command = "$SCHRODINGER/utilities/prepwizard -WAIT -disulfides -fix -noepik -noimpref -noprotassign -reference_st_file %s -NOLOCAL %s %s" %(ref, pdb_name, new_pdb)
	print command
	os.system(command)
	return

def pprep(pdb_dir, ref):
	pdbs = get_trajectory_files(pdb_dir, ext = ".pdb")
	os.chdir(pdb_dir)
	
	pprep_partial = partial(pprep_prot, ref = ref)

	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(pprep_partial, pdbs)
	pool.terminate()


'''
The f ollowing two functions take as input the directory contianing the ligands you would like to dock to your receptors, 
and prepares them with Schrodinger LigPrep, and then saves them in .maegz format, required for the actual docking. 
You can change the settings listed in "ligfile.write" lines. Perhaps we should add this instead as optional inputs
in the function definition. 
'''


def prepare_ligand(lig, lig_dir):
	os.chdir(lig_dir)
	lig_last_name = lig.split("/")[len(lig.split("/"))-1]
	lig_no_ext = lig_last_name.split(".")[0]

	lig_input = "%s/%s.inp" %(lig_dir, lig_no_ext)
	lig_output = "%s-out.maegz" %lig_no_ext

	if os.path.exists("%s/%s" %(lig_dir,lig_output)):
		print "already prepared ligand"
		return

	ligfile = open(lig_input, "wb")
	ligfile.write("INPUT_FILE_NAME   %s \n" %lig_last_name)
	ligfile.write("OUT_MAE   %s \n" %lig_output)
	ligfile.write("FORCE_FIELD   14 \n")
	ligfile.write("EPIK   yes \n")
	ligfile.write("DETERMINE_CHIRALITIES   no \n")
	ligfile.write("IGNORE_CHIRALITIES   no \n")
	ligfile.write("NUM_STEREOISOMERS   32 \n")
	ligfile.write("NUM_RING_CONF   6 \n")
	ligfile.close()

	cmd = "$SCHRODINGER/ligprep -WAIT -inp %s" %lig_input
	print cmd
	subprocess.call(cmd, shell=True)


def prepare_ligands(lig_dir, ext = ".mae"):
	ligs = get_trajectory_files(lig_dir, ext)
	print ligs
	lig_partial = partial(prepare_ligand, lig_dir = lig_dir)

	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(lig_partial, ligs)
	pool.terminate()
	print "finished preparing ligands"

'''
To dock, Schrodinger has to generate grid files (in .zip format) for each receptor. This needs as input the (x,y,z) coordinates 
for the center of the grid, and parameters for the size of the box surrounding that point in which Glide will try to dock your ligand(s).
ALl you need to do is pass to "generate_grids() the following: 
mae_dir, a directory containing mae files of the receptors to which you will dock 
grid_center, a *string* containing the x,y,z coords of the center of the grid, e.g: grid_center = "64.4, 16.9, 11.99"
grid_dir: the directory where you want Schrodinger to save the .zip grid files
remove_lig: if there is a co-crystallized or docked ligand already in your .mae files, you will need to remove it first. to automatically
do this, set remove_lig to the 3-letter upper case string residue name denoting that ligand. for B2AR PDB ID: 3P0G, I would pass: remove_lig = "BIA"

"
'''

def generate_grid_input(mae, grid_center, grid_dir, remove_lig = None):
	mae_name = mae.rsplit( ".", 1)[0]
	mae_last_name = mae_name.split("/")[len(mae_name.split("/"))-1]

	output_dir = grid_dir
	
	new_mae = "%s/%s.mae" %(output_dir, mae_last_name)

	grid_job = "%s/%s.in" %(output_dir, mae_last_name)
	grid_file = "%s/%s.zip" %(output_dir, mae_last_name)

	if (os.path.exists(grid_job) and os.path.exists(new_mae)) or (os.path.exists(grid_file)):
		print "Already created that grid job, skipping"
		return

	if remove_lig == None:
		cmd = "cp %s %s" %(mae, new_mae)
		print cmd
		subprocess.call(cmd, shell=True)
	else:
		cmd = "$SCHRODINGER/run $SCHRODINGER/mmshare-v29013/python/common/delete_atoms.py -asl \"res.pt %s \" %s %s" %(remove_lig, mae, new_mae)
		print cmd
		subprocess.call(cmd, shell=True)

	gridfile = open(grid_job, "wb")
	gridfile.write("GRIDFILE   %s.zip \n" %mae_last_name)
	gridfile.write("OUTPUTDIR   %s \n" %output_dir)
	gridfile.write("GRID_CENTER   %s \n" %grid_center)
	gridfile.write("INNERBOX   10, \n")
	gridfile.write("OUTERBOX   25.0, \n")
	gridfile.write("RECEP_FILE   %s \n" %new_mae)
	#gridfile.write("RECEP_VSCALE   1.0 \n")
	#gridfile.write("LIGAND_MOLECULE  1 \n")
	#gridfile.write("WRITEZIP   TRUE \n")
	gridfile.close()

def generate_grid(grid_file, grid_dir):
	grid_zip = grid_file.rsplit( ".", 1)[0]
	grid_zip = "%s.zip" %grid_zip
	if os.path.exists(grid_zip):
		print "already generated grid; skipping"
		return

	os.chdir(grid_dir)
	grid_command = "$SCHRODINGER/glide %s -OVERWRITE -WAIT" %grid_file

	subprocess.call(grid_command, shell = True)
	print "completed grid generation job"
	return 

def unzip(zip_file):
	output_folder = "/".join(zip_file.split("/")[0:len(zip_file.split("/"))-1])
	os.chdir(output_folder)
	zip_file_last_name = zip_file.split("/")[len(zip_file.split("/"))-1].split(".")[0]
	if os.path.exists("%s/%s.grd" %(output_folder, zip_file_last_name)): return

	cmd = "unzip %s" %zip_file
	subprocess.call(cmd, shell = True)
	return

def generate_grids(mae_dir, grid_center, grid_dir, remove_lig = None):
	print grid_dir
	if not os.path.exists(grid_dir): os.makedirs(grid_dir)

	maes = get_trajectory_files(mae_dir, ".mae")

	generate_grid_input_partial = partial(generate_grid_input, grid_dir = grid_dir, grid_center = grid_center, remove_lig = remove_lig)
	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(generate_grid_input_partial, maes)
	pool.terminate()

	grid_files = get_trajectory_files(grid_dir, ".in")

	grid_partial = partial(generate_grid, grid_dir = grid_dir)

	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(grid_partial, grid_files)
	pool.terminate()

	zips = get_trajectory_files(grid_dir, ".zip")
	pool = mp.Pool(num_workers)
	pool.map(unzip, zips)
	pool.terminate()

	return grid_dir

'''
the function, dock_conformations() is to dock a single ligand to many conformations. you will probably  be using the function,
dock_ligands_and_receptors(), however.
'''

def dock(dock_job):
	cmd = "$SCHRODINGER/glide %s -OVERWRITE -WAIT" %dock_job
	subprocess.call(cmd, shell = True)
	print "Completed docking job %s" %dock_job
	return


def dock_conformations(grid_dir, docking_dir, ligand_dir, precision = "SP", chosen_jobs = False, parallel = False):
	if not os.path.exists(docking_dir): os.makedirs(docking_dir)
	os.chdir(docking_dir)

	#grid_subdirs = [x[0] for x in os.walk(grid_dir)]
	#grid_subdirs = grid_subdirs[1:]
	grid_files = get_trajectory_files(grid_dir, ".grd")
	dock_jobs = []
	for grid_file in grid_files:
		grid_filename = grid_file.split("/")[len(grid_file.split("/"))-1]
		grid_file_no_ext = grid_filename.split(".")[0]

		if chosen_jobs is not False:
			if grid_file_no_ext not in chosen_jobs:
				#print "%s not in chosen jobs " %grid_file_no_ext
				continue
		#print grid_file_no_ext
		maegz_name = "%s/%s_pv.maegz" %(docking_dir, grid_file_no_ext)
		log_name = "%s/%s.log" %(docking_dir, grid_file_no_ext)
		log_size = 0
		if os.path.exists(log_name): log_size = os.stat(log_name).st_size
		if os.path.exists(maegz_name):# and log_size > 3000:
			print "already docked %s" %grid_file_no_ext
			continue
		dock_job_name = "%s/%s.in" %(docking_dir, grid_file_no_ext)
		dock_jobs.append(dock_job_name)

		dock_job_input = open(dock_job_name, "wb")
		dock_job_input.write("GRIDFILE  %s \n" %grid_file)
		dock_job_input.write("LIGANDFILE   %s \n" %ligand_dir)
		if precision == "XP":
			dock_job_input.write("POSTDOCK_XP_DELE   0.5 \n")
		dock_job_input.write("PRECISION   %s \n" %precision)
		if precision == "XP":
			dock_job_input.write("WRITE_XP_DESC   False \n")
		dock_job_input.write("OUTPUTDIR   %s \n" %docking_dir)
		dock_job_input.close()

	print("Written all docking job input files")
	#print dock_jobs

	if parallel:
		num_workers = mp.cpu_count()
		pool = mp.Pool(num_workers)
		pool.map(dock, dock_jobs)
		pool.terminate()
	else:
		for job in dock_jobs:
			dock(job)

	print("Done docking.")

def failed(log_file):
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
	if score == 0: return False 
	return True
def failed_docking_jobs(docking_dir, ligand, precision):
	logs = get_trajectory_files(docking_dir, ext = ".log")
	failed_jobs = []

	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	job_results = pool.map(failed, logs)
	pool.terminate()

	for i in range(0,len(logs)):
		if job_results[i] == False:
			failed_jobs.append(logs[i])

	failed_jobs = [j.split("/")[len(j.split("/"))-1].split(".")[0] for j in failed_jobs]
	return failed_jobs

def dock_helper(args):
	dock_conformations(*args)

class NoDaemonProcess(mp.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(mp.pool.Pool):
    Process = NoDaemonProcess

'''
This is the function for docking multiple ligands to multiple receptors. 

grid_dir: the directory where the .zip files for each receptor to whicih you would like to dock is located. 
docking_dir = the directory to which you would like Glide to save all the results of docking 
ligands_dir = the directory containing the .maegz files containing the LigPrep prepared ligands for docking 
precision --> each Glide docking job can be done in SP or XP level of precision. XP is more accurate but takes about 7 times as long as each SP calculations
	you can change this to precision = "XP" if you would like to try that. Literature shows that it is in fact a little more accurate.
chosen_ligands --> if, in your ligands_dir directory you only want to dock particular ligands, pass here a list of strings of the ligand names,
	and it will only dock those ligands. in the example folder provided, for example, if you pass ["procaterol", "ta-2005"], it will only dock
	those two ligands
chosen_receptors --> same as chosen_ligands. if you pass ["cluster301_sample0", "cluster451_sample5"] it will only use those two receptors for docking
parallel --> if you set it to "both" it will run in parallel over both ligands and receptors. I don't recommend this generally. 
	if you pass "ligand": it will parallelize over all ligands. Recommened if n_liagnds > n_receptors
	if you pass "receptor": it will parallelize over receptors. Recommedned if n_receptors > n_ligands
'''

def dock_ligands_and_receptors(grid_dir, docking_dir, ligands_dir, precision = "SP", ext = "-out.maegz", chosen_ligands = False, chosen_receptors = False, parallel = False):
	ligands = get_trajectory_files(ligands_dir, ext = ext)

	if parallel == "both":
		lig_dirs = []
		docking_dirs = []
		args = []
		for ligand in ligands:
			lig_last_name = ligand.split("/")[len(ligand.split("/"))-1]
			lig_no_ext = lig_last_name.split("-out.")[0]
			if chosen_ligands is not False:
				if lig_no_ext not in chosen_ligands: continue
			lig_dir = "%s/%s" %(docking_dir, lig_no_ext)
			if not os.path.exists(lig_dir): os.makedirs(lig_dir)
			docking_dirs.append(lig_dir)
			lig_dirs.append(ligand)
			args.append((grid_dir, lig_dir, ligand, precision, chosen_receptors, True))
		
		num_workers = 5
		pool = MyPool(num_workers)
		pool.map(dock_helper, args)
		pool.terminate()

	elif parallel == "ligand":
		lig_dirs = []
		docking_dirs = []
		args = []
		for ligand in ligands:
			lig_last_name = ligand.split("/")[len(ligand.split("/"))-1]
			lig_no_ext = lig_last_name.split("-out.")[0]
			if chosen_ligands is not False:
				if lig_no_ext not in chosen_ligands: continue
			lig_dir = "%s/%s" %(docking_dir, lig_no_ext)
			if not os.path.exists(lig_dir): os.makedirs(lig_dir)
			docking_dirs.append(lig_dir)
			lig_dirs.append(ligand)
			args.append((grid_dir, lig_dir, ligand, precision, chosen_receptors, False))
		
		num_workers = mp.cpu_count()
		pool = mp.Pool(num_workers)
		pool.map(dock_helper, args)
		pool.terminate()

	elif parallel == "receptor":
		lig_dirs = []
		docking_dirs = []
		args = []
		for ligand in ligands:
			lig_last_name = ligand.split("/")[len(ligand.split("/"))-1]
			lig_no_ext = lig_last_name.split("-out.")[0]
			if chosen_ligands is not False:
				if lig_no_ext not in chosen_ligands: continue
			lig_dir = "%s/%s" %(docking_dir, lig_no_ext)
			if not os.path.exists(lig_dir): os.makedirs(lig_dir)
			docking_dirs.append(lig_dir)
			lig_dirs.append(ligand)
			args.append((grid_dir, lig_dir, ligand, precision, chosen_receptors, True))
		
		for arg in args:
			dock_helper(arg)

	else:
		for ligand in ligands:
			lig_last_name = ligand.split("/")[len(ligand.split("/"))-1]
			lig_no_ext = lig_last_name.split("-out.")[0]
			if chosen_ligands is not False:
				if lig_no_ext not in chosen_ligands: continue
			lig_dir = "%s/%s" %(docking_dir, lig_no_ext)
			if not os.path.exists(lig_dir): os.makedirs(lig_dir)
			dock_conformations(grid_dir, lig_dir, ligand, precision = precision, chosen_jobs = chosen_receptors)

'''

Identical as above functions for docking, but for MM-GBSA calculations 
'''

def mmgbsa_individual(job):
	cmd = "$SCHRODINGER/prime_mmgbsa -WAIT %s" %job
	subprocess.call(cmd, shell = True)
	print "Completed mmgbsa job %s" %job
	return

def mmgbsa(docking_dir, mmgbsa_dir, chosen_jobs = False):
	if not os.path.exists(mmgbsa_dir): os.makedirs(mmgbsa_dir)
	os.chdir(mmgbsa_dir)

	dock_files = get_trajectory_files(docking_dir, "pv.maegz")
	mmgbsa_jobs = []
	for dock_file in dock_files:
		dock_filename = dock_file.split("/")[len(dock_file.split("/"))-1]
		dock_file_no_ext = dock_filename.rsplit(".", 1)[0]
		dock_file_no_pv = dock_file_no_ext.split("_pv")[0]
		if chosen_jobs is not False:
			if dock_file_no_pv not in chosen_jobs:
				continue
		mmgbsa_out_name = "%s/%s-out.maegz" %(mmgbsa_dir, dock_file_no_pv)
		log_name = "%s/%s.log" %(mmgbsa_dir, dock_file_no_pv)
		log_size = 0
		if os.path.exists(log_name): log_size = os.stat(log_name).st_size

		if os.path.exists(mmgbsa_out_name) and log_size > 2599:
			print "Already ran mmgbsa with %s" %dock_file_no_pv
			continue
		cmd = "cp %s %s" %(dock_file, mmgbsa_dir)
		os.system(cmd)
		job_name = "%s/%s.inp" %(mmgbsa_dir, dock_file_no_ext)
		mmgbsa_jobs.append(job_name)
		job_input = open(job_name, "wb")
		job_input.write("STRUCT_FILE %s \n" %dock_filename)
		job_input.write("OUT_TYPE COMPLEX \n")
		job_input.write("FLEXDIST 5.0 \n")
		job_input.write("OVERWRITE \n")
		job_input.close()

	print("Written all mmgbsa input files")

	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(mmgbsa_individual, mmgbsa_jobs)
	pool.terminate()
	print "Done with MM GBSA calculations"

def mmgbsa_ligands_and_receptors(docking_dir, mmgbsa_dir, ligands, chosen_receptors = False):
	for ligand in ligands:
		lig_dir = "%s/%s" %(docking_dir, ligand)
		lig_mmgbsa_dir = "%s/%s" %(mmgbsa_dir, ligand)
		mmgbsa(lig_dir, lig_mmgbsa_dir, chosen_jobs = chosen_receptors)




