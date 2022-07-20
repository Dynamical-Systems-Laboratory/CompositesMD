# ------------------------------------------------------------------
#
#	Computations for binned densities from LAMMPS
#
# ------------------------------------------------------------------

import sys
py_path = './io_operations/'
sys.path.insert(0, py_path)

import math
import load_bin_densities as lbd

def compute_ratios(file_num, file_den, times=[], all_t=True): 
	''' Collects number or mass densities from density*.out files
			and computes their ratios
   
    	file_num -  name of the file with data for the numerator of the ratio 
    	file_den -  name of the file with data for the denominator of the ratio
		times - optional list of select times for which (only) to collect
					the data; need to be integers 
		all_t - collect all time steps provided, has to be set to False
					to collect the data from time steps listed in times array only

		Returns a nested list of ratios of the mass or number densities in each bin at every
		time step i.e. [[bin_1(t0), ..., bin_N(t0)], ..., [bin_1(tK), ..., bin_N(tK)]]

		Return values are floating point numbers; 0 in the denominator appears as a 0 entry 
		in the result.
	'''

	# Number of header lines at the top of the file
	nhd = 3

	# Collect the data
	if all_t == True:
		data_num = lbd.load_bin_density(file_num, nhd)
		data_den = lbd.load_bin_density(file_den, nhd)
	else:
		data_num = lbd.load_bin_density(file_num, nhd, times, False)
		data_den = lbd.load_bin_density(file_den, nhd, times, False)

	# Compute the ratios
	ratios = []
	for row_num, row_den in zip(data_num,  data_den):
		temp_ratios = []
		for bin_num, bin_den in zip(row_num, row_den):
			if math.isclose(bin_den, 0.0, rel_tol=1e-16):
				temp_ratios.append(0.0)
			else:
				temp_ratios.append(bin_num/bin_den)

		ratios.append(temp_ratios)

	return ratios
	

