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
import pytraj.io as mdio
from pytraj import adict

'''
If sim was run under periodic boundary conditions, this will reimage the trajectory
'''

def reimage_traj(traj_file, traj_dir, save_dir, ext):
	if ext == ".pdb":
		traj_pytraj = mdio.load(traj_file, top = traj_file)[:]
		#traj_pytraj.fixatomorder()
		traj_pytraj.autoimage()

		file_lastname = traj_file.split("/")[len(traj_file.split("/"))-1]
		filename = file_lastname.split(".")[0]
		h5_filename = file_lastname
		new_h5_file = "%s/%s" %(save_dir, h5_filename)
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
'''

def pprep_prot(pdb):
	pdb_name = pdb.split("/")[len(pdb.split("/"))-1]
	new_pdb = pdb_name.rsplit( ".", 1 )[ 0 ]
	new_pdb = "%s.mae" %(new_pdb)
	if os.path.exists(new_pdb): 
		print "already prepped and mae'd protein"
		return
	ref = "/scratch/users/enf/b2ar_analysis/3P0G_pymol_prepped.pdb"

	command = "$SCHRODINGER/utilities/prepwizard -WAIT -disulfides -fix -noepik -noimpref -noprotassign -reference_st_file %s -NOLOCAL %s %s" %(ref, pdb_name, new_pdb)
	print command
	os.system(command)
	return

def pprep(pdb_dir):
	pdbs = get_trajectory_files(pdb_dir, ext = ".pdb")
	os.chdir(pdb_dir)
	
	num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(pprep_prot, pdbs)
	pool.terminate()

'''
INPUT_FILE_NAME   ligprep_5.sdf
OUT_MAE   ligprep_5-out.maegz
FORCE_FIELD   14
EPIK   yes
DETERMINE_CHIRALITIES   no
IGNORE_CHIRALITIES   no
NUM_STEREOISOMERS   32
NUM_RING_CONF   6


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

def generate_grid_input(mae, grid_center, tica_dir, n_clusters, n_samples, grid_dir):
	mae_name = mae.rsplit( ".", 1)[0]
	mae_last_name = mae_name.split("/")[len(mae_name.split("/"))-1]

	output_dir = grid_dir
	new_mae = "%s/%s.mae" %(output_dir, mae_last_name)

	grid_job = "%s/%s.in" %(output_dir, mae_last_name)
	grid_file = "%s/%s.zip" %(output_dir, mae_last_name)

	if (os.path.exists(grid_job) and os.path.exists(new_mae)) or (os.path.exists(grid_file)):
		print "Already created that grid job, skipping"
		return

	cmd = "$SCHRODINGER/run $SCHRODINGER/mmshare-v29013/python/common/delete_atoms.py -asl \"res.pt BIA \" %s %s" %(mae, new_mae)
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

def generate_grids(mae_dir, grid_center, tica_dir, n_clusters, n_samples, grid_dir):
	if not os.path.exists(grid_dir): os.makedirs(grid_dir)

	maes = get_trajectory_files(mae_dir, ".mae")

	generate_grid_input_partial = partial(generate_grid_input, grid_dir = grid_dir, grid_center = grid_center, tica_dir = tica_dir, n_clusters = n_clusters, n_samples = n_samples)
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

	return grid_dir

def dock(dock_job):
	cmd = "$SCHRODINGER/glide %s -OVERWRITE -WAIT" %dock_job
	subprocess.call(cmd, shell = True)
	print "Completed docking job %s" %dock_job
	return


def dock_conformations(grid_dir, docking_dir, ligand_dir, precision = "SP", chosen_jobs = False, num_workers = False):
	if not os.path.exists(docking_dir): os.makedirs(docking_dir)
	os.chdir(docking_dir)

	grid_files = get_trajectory_files(grid_dir, ".zip")
	dock_jobs = []
	for grid_file in grid_files: 
		grid_filename = grid_file.split("/")[len(grid_file.split("/"))-1]
		grid_file_no_ext = grid_filename.rsplit(".", 1)[0]
		if chosen_jobs is not False:
			if grid_file_no_ext not in chosen_jobs:
				print "%s not in chosen jobs " %grid_file_no_ext
				continue
		maegz_name = "%s/%s_pv.maegz" %(docking_dir, grid_file_no_ext)
		log_name = "%s/%s.log" %(docking_dir, grid_file_no_ext)
		log_size = 0
		if os.path.exists(log_name): log_size = os.stat(log_name).st_size
		if os.path.exists(maegz_name) and log_size > 3000:
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

	if num_workers is False: num_workers = mp.cpu_count()
	pool = mp.Pool(num_workers)
	pool.map(dock, dock_jobs)
	pool.terminate()

	print("Done docking.")

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

def dock_ligands_and_receptors(grid_dir, docking_dir, ligands_dir, precision = "SP", ext = "-out.maegz", chosen_ligands = False, chosen_receptors = False, parallel = False):
	ligands = get_trajectory_files(ligands_dir, ext = ext)

	if parallel:
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
			args.append((grid_dir, lig_dir, ligand, precision, chosen_receptors, 2))
		
		num_workers = 5
		pool = MyPool(num_workers)
		pool.map(dock_helper, args)
		pool.terminate()

	else:

		for ligand in ligands:
			lig_last_name = ligand.split("/")[len(ligand.split("/"))-1]
			lig_no_ext = lig_last_name.split("-out.")[0]
			if chosen_ligands is not False:
				if lig_no_ext not in chosen_ligands: continue
			lig_dir = "%s/%s" %(docking_dir, lig_no_ext)
			if not os.path.exists(lig_dir): os.makedirs(lig_dir)
			dock_conformations(grid_dir, lig_dir, ligand, precision = precision, chosen_jobs = chosen_receptors)

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




