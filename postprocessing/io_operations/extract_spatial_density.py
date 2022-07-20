# ------------------------------------------------------------------
#
#   Tools for extraction of data from file with binned mass 
#		density of all species as a function of time
#
# ------------------------------------------------------------------

def collect_density_one_dir(fname):
	''' Extracts density data from file generated with chunk/ave in x direction '''

	# fname - name of the file with the input data
	# Returns a dictionary of nested lists with time as a key (converted
	#	to string) and value a sublist per bin [x coordinate of bin, density]
	#	both values converted to floats

	data = {}
	with open(fname, 'r') as fin:
		# Skip the header
		for i in range(3):
			next(fin)
		# Extract
		for line in fin:
			temp = line.strip().split()
			# First element of first line is time
			time = temp[0]
			data[time] = []
			# Second element of first line is bin number
			nbins = int(temp[1])
			for ibin in range(nbins):
				# Extract all the bin information
				line = next(fin)
				sub_temp = line.strip().split()
				data[time].append([float(sub_temp[1]), float(sub_temp[-1])])
	return data

if __name__ == '__main__':
	
	fname = '../../equilibration/DPD_based/water_models/seed_1/density.out' 
	
	data = collect_density_one_dir(fname)

	print(data['10151500'])


