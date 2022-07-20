# ------------------------------------------------------------------
#
#   Tools for extraction of density data from density*.out files
#
# ------------------------------------------------------------------

def load_bin_density(fname, nhd, times=[], all_t=True):
	''' Collects number or mass densities from density*.out files
   
    	fname -  name of the file with the density data 
		nhd - number of header lines
		times - optional list of select times for which (only) to collect
					the data; need to be integers 
		all_t - collect all time steps provided, has to be set to False
					to collect the data from time steps listed in times array only

		Returns a nested list of mass or number densities in each bin at every
		time step i.e. [[bin_1(t0), ..., bin_N(t0)], ..., [bin_1(tK), ..., bin_N(tK)]]

		Return values are floating point numbers
	'''

	data_out = []
	with open (fname, 'r') as fin:
		# Skip the headers
		for ih in range(nhd):
			next(fin)

		for line in fin:
			cur_step = []

			# Should always be the line with time and bin value at this point
			temp_line = line.strip().split()

			if all_t == False:
				# Skip if not the desired step
				if not (int(temp_line[0]) in times):
					continue

			# Number of recorded bins 
			nbins = int(temp_line[1])
			
			# Collect data in all nbins
			for ib in range(nbins):
				line = next(fin)
				temp_line = line.strip().split()
				# Density is the last value in a row 
				cur_step.append(float(temp_line[-1]))

			data_out.append(cur_step)	

	return data_out
