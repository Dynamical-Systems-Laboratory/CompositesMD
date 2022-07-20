# ------------------------------------------------------------------
#
#	Computations for assymetry factors Sf and Qf 
#
# ------------------------------------------------------------------

import sys
py_path = './io_operations/'
sys.path.insert(0, py_path)

import math, warnings

def compute_sf(sdist_atom):
	''' Compute Sf factor - measure of spatial asymmetry 
			sdist_atom is a nested list containing some form 
			of binned data on spatial distribution of this atom 
			type or types e.g. number of atoms in each bin. 
			The nested lists are the bins at different time steps.
			Returns Sf (sfactor) for every time step i.e. nested list. '''

	sfactor = []
	for sdtime in sdist_atom:
		# Even number of bins is better
		if len(sdtime)%2:
			warnings.warn('Number of bins is not even. Central bin will not be considered in the computation.')
		
		# For one time step
		lhs = 0
		rhs = len(sdtime)-1
		# Difference between mirror bins
		Sij = []
		while lhs < rhs:
			Sij.append(abs(sdtime[rhs]-sdtime[lhs]))
			lhs += 1
			rhs -= 1

		# Sf factor for that time step
		sfactor.append(math.sqrt(sum(x**2 for x in Sij)))

	return sfactor

def compute_qf(sdist_neg, sdist_pos):
	''' Compute Qf factor - measure of net charge asymmetry.

			sdist_neg and sdist_pos are binned spatial distributions
			of negative and positive ions.
			Each is a nested list containing some form 
			of binned data on spatial distribution of this atom 
			type or types e.g. number of atoms in each bin. 
			The nested lists are the bins at different time steps.
	
			Returns Qf (qfactor) for every time step in a list. '''

	qfactor = []
	for sneg, spos in zip(sdist_neg, sdist_pos):
		# Both datasets need to have equal lenghts
		if len(sneg) != len(spos):
			raise ValueError('Binned data for positive and negative ions should have the same length')
		# Even number of bins is better
		if len(sneg)%2:
			warnings.warn('Number of bins is not even. Central bin will not be considered in the computation.')
		
		# For one time step
		lhs = 0
		rhs = len(sneg)-1
		
		# Net charge in each bin 
		Qi = [x - y for x,y in zip(spos, sneg)]

		# Difference of net charges between mirror bins
		Qij = []
		while lhs < rhs:
			Qij.append(abs(Qi[rhs]-Qi[lhs]))
			lhs += 1
			rhs -= 1

		# Qf factor for that time step
		qfactor.append(math.sqrt(sum(x**2 for x in Qij)))
		
	return qfactor

def phi_double_sum(sdist_neg, sdist_pos, len_box):
	''' Compute the difference between approximate potential in LHS and RHS
			of the matrix - measure of net charge asymmetry.
			
			sdist_neg, sdist_pos, and len_box are binned spatial distributions
			of negative and positive ions, and box dimension in the 
			deformation direction.

			Spatial distributions are nested lists containing some form 
			of binned data on spatial distribution of this atom 
			type or types e.g. number of atoms in each bin.
			The nested lists are the bins at different time steps.

			len_box is a list with lenght at a given time step.

			Returns potential difference for every time step in a list. '''

	diff_fi = []
	for sneg, spos, len_x in zip(sdist_neg, sdist_pos, len_box):
		# Both datasets need to have equal lenghts
		if len(sneg) != len(spos):
			raise ValueError('Binned data for positive and negative ions should have the same length')
		# Even number of bins is better
		if len(sneg)%2:
			raise ValueError('Use even number of bins')
		
		# Net charge in each bin 
		Qi = [x - y for x,y in zip(spos, sneg)]

		# Sum of net charges in each bin in each half of the membran
		beta_L = sum(Qi[0:int(len(sneg)/2)])  
		beta_R = sum(Qi[int(len(sneg)/2):])

		# Difference between net charges
		diff_fi.append(len(sneg)/8.0*len_x*len_x*abs(beta_L-beta_R))

	return diff_fi

def phi_spatial(sdist_neg, sdist_pos, len_x, len_y, len_z):
	''' Compute the approximate potential in each spatial bin 			
			sdist_neg and sdist_pos are binned spatial distributions
			of negative and positive ions in the deformation direction.
			
			len_x, len_y, and len_z are the bin x, y, and z dimensions 
				respectively

			Spatial distributions are nested lists containing some form 
			of binned data on spatial distribution of this atom 
			type or types e.g. number of atoms in each bin.
			The nested lists are the bins at different time steps.

			Returns potential in each bin for every time step in a nested list. '''

	time_fi = []
	for sneg, spos, dxi, ly, lz in zip(sdist_neg, sdist_pos, len_x, len_y, len_z):
		# Both datasets need to have equal lenghts
		if len(sneg) != len(spos):
			raise ValueError('Binned data for positive and negative ions should have the same length')
		
		# Net charge in each bin 
		Qi = [(x - y) for x,y in zip(spos, sneg)]
		
		bin_fi = [0]*len(sneg)
		# Compute the double sum for each bin
		for ib in range(len(sneg)):
			bin_fi[ib] = (ib+1)*dxi*sum(Qi[0:ib+1])/(ly*lz)
		
		# Append the complete solution 
		time_fi.append(bin_fi)

	return time_fi



