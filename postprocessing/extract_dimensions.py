# ------------------------------------------------------------------
#
#	Extract box dimensions from LAMMPS	
#
# ------------------------------------------------------------------

import sys
py_path = './io_operations/'
sys.path.insert(0, py_path)

import math
import select_d_section as dsec
import post_utils as utils

def extract_box_dims_and_lims(fname):
	''' 
		fname is the name of the dump file
		Returns a 4 component nested list: 
			dim[0] - lower x limits (all times)
			dim[1] - upper x limits (all times)
			dim[2] - box dimension in y (all times)
			dim[3] - box dimension in z (all times)
	
		All units are the simulation defaults
	'''

	# Extract box dimensions
	box_dims = dsec.extract_all_sections(fname, 'ITEM: BOX')
	
	# Collect the data
	x0 = []
	xF = []
	dy = []
	dz = []
	for box in box_dims:
		flt_box = utils.conv_to_float(box)
		x0.append(flt_box[0][0])
		xF.append(flt_box[0][1]) 
		dy.append(flt_box[1][1] - flt_box[1][0]) 
		dz.append(flt_box[2][1] - flt_box[2][0])	
	
	data = []
	data.append(x0)
	data.append(xF)
	data.append(dy)
	data.append(dz)

	return data

def extract_box_lims_last(fname):
	''' 
		fname is the name of the dump file
		Returns a 6 component list: 
			dim[0] - lower x limits 
			dim[1] - upper x limits 
			dim[2] - lower y limits 
			dim[3] - upper y limits 
			dim[4] - lower z limits 
			dim[5] - upper z limits 
		
		All units are the simulation defaults
	'''

	# Extract the last box dimensions
	box = dsec.extract_all_sections(fname, 'ITEM: BOX')[-1]
	flt_box = utils.conv_to_float(box)
	data = []
	data.append(flt_box[0][0])	
	data.append(flt_box[0][1])
	data.append(flt_box[1][0])
	data.append(flt_box[1][1])
	data.append(flt_box[2][0])
	data.append(flt_box[2][1])

	return data

def compute_volume_change_ratio(pre_dfile, dfile, nbins, types):
	''' 
		Compute volume change ratio, Vt/V0 for all the 
		time steps recorded in dfile dump file and nbins

		V0 is from the last recorded time in pre_file
		dump file - assumed equilibrated, pre-deformed
		configuration 
		
		Accounts only for atom types in types list

		Currently only for bins in x direction

		Returns nested list where each sublist has 
		the ratios for each bin at that time step
	
	'''

	# --- Initial coordinates
	# Collect the last step in pre_dfile
	# Extract atom data
	atoms_0 = dsec.extract_all_sections(pre_dfile, 'ITEM: ATOMS')[-1]	
	# Extract box dimensions
	box_0 = dsec.extract_all_sections(pre_dfile, 'ITEM: BOX')[-1]
	# x direction
	b_ind = 0 
	d_ind = 2
	dims = box_0[b_ind].strip().split()
	Lx_0 = float(dims[0])
	Lx_tot0 = float(dims[1]) - Lx_0
	# y and z for adjusting dimensions
#	Ly_0 = float(box_0[1].strip().split()[0])
#	Lz_0 = float(box_0[2].strip().split()[0])
	# Map atom IDs and coordinates for species listed in types
	all_atoms_0 = {}
	for at in atoms_0:
		at = at.strip().split()
		if not(at[1] in types):
			continue
#		all_atoms_0[at[0]] = [float(at[2]) - Lx_0, float(at[3]) - Ly_0, float(at[4]) - Lz_0]
		all_atoms_0[at[0]] = [float(at[2]), float(at[3]), float(at[4])]

	# --- Volume ratio for each time step
	# Load the dfile - all simulation steps 
	# Extract atom data
	atoms_all_t = dsec.extract_all_sections(dfile, 'ITEM: ATOMS')
	# Extract box dimensions
	box_t = dsec.extract_all_sections(dfile, 'ITEM: BOX')

	# Each sublist will be the J factor for each bin at that step 
	volume_ratios = []
	for atoms, box in zip(atoms_all_t, box_t):
		# Create bins, compute their volumes
		dims = box[b_ind].strip().split()
		Lx = float(dims[0])
		Lx_tot = float(dims[1]) - Lx
		bin_width = Lx_tot/nbins
		# y and z for adjusting dimensions
		dims = box_0[1].strip().split()
		Ly = float(dims[0])
		Ly_tot = float(dims[1]) - float(dims[0])
		dims = box_0[2].strip().split()
		Lz = float(dims[0])
		Lz_tot = float(dims[1]) - float(dims[0])
		# Simple for now
		V_bin = bin_width*Ly_tot*Lz_tot 

		# Find atoms in each bin
		bin_atoms = {}
		for at in atoms:
			at = at.strip().split()
			if not(at[1] in types):
				continue
			# So that the binning is simple 
			# Shifted the interval to 0->Lx_tot
			norm_pos = float(at[d_ind]) - Lx
			bin_num = str(max(min(math.floor(norm_pos/bin_width),nbins-1), 0))
			# Store ID as a string
			if bin_num in bin_atoms: 
				bin_atoms[bin_num].append((at[0], float(at[2]), float(at[3]), float(at[4])))
			else:
				bin_atoms[bin_num] = [(at[0], float(at[2]), float(at[3]), float(at[4]))]

		# Map the atoms in each bin to old bins, determine the old bin 
		# coordinates and volume
		F_bin = []
		for bn, atom_list in bin_atoms.items():
			# x_min, x_max, y_min, .., z_max for each bin
			bin_dims = [1000, -1000, 1000, -1000, 1000, -1000]
			for aID, xc, yc, zc in atom_list:
				# Original coordinates of this atom
				xc0, yc0, zc0 = all_atoms_0[aID]
				# Correct for periodicity in the bin direction
				# ------ is this right?
				if (abs(xc0-xc) > abs(xc0 - Lx_tot - xc)):
					# If the atom crossed the PBC
					dx = xc0 - Lx_tot - xc
					xc0 = xc + dx + Lx_tot
				# Find if on the edge
				# x limits
				if xc0 < bin_dims[0]:
					bin_dims[0] = xc0		
				if xc0 > bin_dims[1]:
					bin_dims[1] = xc0	
				# y limits
				if yc0 < bin_dims[2]:
					bin_dims[2] = yc0		
				if yc0 > bin_dims[3]:
					bin_dims[3] = yc0
				# z limits
				if zc0 < bin_dims[4]:
					bin_dims[4] = zc0		
				if zc0 > bin_dims[5]:
					bin_dims[5] = zc0					
	
			# Compute and save volume of this bin 
			V0_bin = (bin_dims[1]-bin_dims[0])*(bin_dims[3]-bin_dims[2])*(bin_dims[5]-bin_dims[4])
			F = V_bin/V0_bin
			F_bin.append(F)
		
		# Append the ratio for all bins in this time step
		volume_ratios.append(F_bin)		
	
	# Return nested list	
	return volume_ratios

if __name__ == '__main__':
	
	pre_fname = '/home/user/Research/nafion/MDNafion/simulations/stress_driven_sensing/seed_1/pre_nafion.d'
	fname = '/home/user/Research/nafion/MDNafion/simulations/stress_driven_sensing/seed_1/compression_nafion.d'
	nbins = 10
	types = ['1', '2', '3', '4', '5']
	
	volume_ratios = compute_volume_change_ratio(pre_fname, fname, nbins, types)
	print(volume_ratios[0])
	print(volume_ratios[-1])
	print(len(volume_ratios))










