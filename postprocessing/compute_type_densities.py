# ------------------------------------------------------------------
#
#	Spatial distribution of density compuations
#
# ------------------------------------------------------------------

import math 

import sys
py_path = './io_operations/'
sys.path.insert(0, py_path)

import select_d_section as dsec
import extract_spatial_density as esd 
import post_utils as utils

def compute_spatial_density(fname, typeID, time, nbins, idir, wpos=[], all_steps=False, den=False):
	''' Compute spatial density distribution of given atom type ''' 

	#
	# fname - name of the d file 
	# typeID - atom type as defined in .data file
	# time - array with simulation times to plot
	# nbins - number of bins in the target direction
	# idir - string representing direction
	# wpos - positions of walls in a non-periodic system 
	#

	# Extract time steps
	time_steps = dsec.extract_all_sections(fname, 'ITEM: TIMESTEP')

	# Extract atom data
	atoms_all_t = dsec.extract_all_sections(fname, 'ITEM: ATOMS')	

	# Extract box dimensions
	box_all_t = dsec.extract_all_sections(fname, 'ITEM: BOX')

	# Direction settings
	# b_ind - index of coordinates in the BOX data
	# d_ind - column with atom positions in that direction in the .d file
	if idir == 'x':
		b_ind = 0
		d_ind = 2 
	elif idir == 'y':
		b_ind = 1
		d_ind = 3
	else:
		b_ind = 2
		d_ind = 4

	densities = []
	for step, atoms, box in zip(time_steps, atoms_all_t, box_all_t):
	
		if (not (int(step[0]) in time)) and (all_steps == False):
			continue
	
		# Spatial bins
		dims = box[b_ind].strip().split()
		if wpos:
			L_0 = wpos[0]
			L_tot = wpos[1] - wpos[0]
		else:
			L_0 = float(dims[0])
			L_tot = float(dims[1]) - L_0
		bin_width = L_tot/nbins
		number_density = [0]*nbins		

		# Calculate how many atoms in each bin	
		atom_data = [] 
		for at in atoms:
			at = at.strip().split()
			if at[1] == str(typeID):
				# So that the binning is simple 
				# Shifted the interval to 0->L_tot
				norm_pos = float(at[d_ind]) - L_0
				number_density[min(max(math.floor(norm_pos/bin_width), 0),nbins-1)] += 1
		
		# Compute volume in cm^3 then convert to A^3
		cm32A3 = 1.0e24 
		flt_box = utils.conv_to_float(box)
	
		vol = utils.compute_volume(flt_box)*cm32A3		
		bin_vol = vol/nbins
		
		if den:
			# Divide by bin volumes (all volumes equal)
			densities.append([x/bin_vol for x in number_density])
		else:
			densities.append([x for x in number_density])
	return densities

				
def compute_type_numbers(dfile_name, bfile_name, typeID, time, den_lim, all_steps = False):
	''' Compute number of atoms of given atom type in a membrane region defined by densities ''' 

	#
	# dfile_name - name of the d file 
	# bfile_name - name of the file with binned densities from chunk/ave
	# typeID - atom type as defined in .data file
	# time - list with simulation times to return 
	# den_lim - minimum density that defines the membrane region
	# all_steps - True if process all the time steps
	# 
	# Returns a list of values for each requested time 
	#

	# Extract time steps
	time_steps = dsec.extract_all_sections(dfile_name, 'ITEM: TIMESTEP')

	# Extract atom data
	atoms_all_t = dsec.extract_all_sections(dfile_name, 'ITEM: ATOMS')	

	# Extract box dimensions
	box_all_t = dsec.extract_all_sections(dfile_name, 'ITEM: BOX')

	# Collect the binned density data
	den_bins = esd.collect_density_one_dir(bfile_name)

	atom_count = []
	for step, atoms, box in zip(time_steps, atoms_all_t, box_all_t):
	
		if (not (int(step[0]) in time)) and (all_steps == False):
			continue

		# Box lenghth in x 
		dims = box[0].strip().split()
		Lx = float(dims[1]) - float(dims[0])
		
		# Find the region spanned by the memmbrane based
		# on the density bounds at this or closest time step
		step_i = int(step[0].strip())
		dt_min = 1000000.0
		dt_key = 'wrong key'
		for key in den_bins.keys():
			dt = abs(step_i-int(key))
			if dt < dt_min:
				dt_key = key
		bin_data = den_bins[dt_key]
		dbin = bin_data[1][0] - bin_data[0][0]

		# For the region - take into account possible crossing
		# of the periodic boundary; recognizes crossing by 
		# checking bin distances for density values that 
		# meet the criterion >= den_lim
		mem_reg = [x[0] for x in bin_data if x[1] >= den_lim]
		# Flag for later use
		is_crossing = False
		# Not crossing if sequential
		is_sequence = [i for i in range(len(mem_reg)-1) if math.isclose(mem_reg[i+1]-mem_reg[i], dbin)]
		if (len(is_sequence) == len(mem_reg)):
			# Not crossing
			x_min = mem_reg[0]*Lx
			x_max = mem_reg[-1]*Lx
		else:
			# Crossing
			water_reg = [x[0] for x in bin_data if x[1] < den_lim] 
			x_min = (water_reg[0]-dbin)*Lx
			x_max = (water_reg[1]+dbin)*Lx
			is_crossing = True

		# Calculate how many atoms of typeID are in the region
		count_i = 0
		for at in atoms:
			at = at.strip().split()
			if at[1] == str(typeID):
				x_pos_at = float(at[2])
				if not is_crossing:
					if  x_pos_at >= x_min and x_pos_at <= x_max:
						count_i += 1
				else:
					if x_pos_at <= x_min or x_pos_at >= x_max:
						count_i += 1	

		atom_count.append(count_i)

	return atom_count 

				



