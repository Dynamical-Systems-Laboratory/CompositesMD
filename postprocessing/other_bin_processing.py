# ------------------------------------------------------------------
#
#	Computations related to various binned quantities 
#
# ------------------------------------------------------------------

import math 

import sys
py_path = './io_operations/'
sys.path.insert(0, py_path)

import select_d_section as dsec
import post_utils as utils

def compute_spatial_vx(fname, time, mass_dict, nbins, wpos=[], all_steps=False):
	''' Compute spatial distribution of average m(vx)^2 '''

	#
	# fname - name of the d file 
	# time - array with simulation times to plot
	# mass_dict - string of atom type : atom mass pairs
	# nbins - number of bins in the target direction
	# wpos - positions of walls in a non-periodic system
	# all_steps - include all time steps
	#

	# Extract time steps
	time_steps = dsec.extract_all_sections(fname, 'ITEM: TIMESTEP')
	# Extract atom data
	atoms_all_t = dsec.extract_all_sections(fname, 'ITEM: ATOMS')	
	# Extract box dimensions
	box_all_t = dsec.extract_all_sections(fname, 'ITEM: BOX')

	x_vel_all = []
	for step, atoms, box in zip(time_steps, atoms_all_t, box_all_t):
	
		if (not (int(step[0]) in time)) and (all_steps == False):
			continue
	
		# Spatial bins
		dims = box[0].strip().split()
		if wpos:
			L_0 = wpos[0]
			L_tot = wpos[1] - wpos[0]
		else:
			L_0 = float(dims[0])
			L_tot = float(dims[1]) - L_0
		bin_width = L_tot/nbins
		atoms_in_bins = [0]*nbins
		temp_vx = [0]*nbins

		# Sum x velocities in each bin	
		for at in atoms:
			at = at.strip().split()
			# So that the binning is simple 
			# Shifted the interval to 0->L_tot
			norm_pos = float(at[2]) - L_0
			ind = max(min(math.floor(norm_pos/bin_width),nbins-1), 0)
			atoms_in_bins[ind] += 1
			# x velocity component
			temp_vx[ind] += mass_dict[at[1]]*float(at[5])*float(at[5])

		# Divide by number of atoms in each bin 
		ave_vx = [x/max(y,1) for x,y in zip(temp_vx, atoms_in_bins)]
		x_vel_all.append(ave_vx)
		
	return x_vel_all

