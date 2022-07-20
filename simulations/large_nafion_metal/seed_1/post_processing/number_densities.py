import sys
py_path = '../../../../postprocessing/io_operations/'
sys.path.insert(0, py_path)

import sys
py_path = '../../../../postprocessing/'
sys.path.insert(0, py_path)

import glob, os
import compute_type_densities as ctd
import io_module as io

#
# Input
#

# Input files
dfile = '../nafion.d'
# Output files
ofile = 'number_density_'

# Remove all old files because this code appends
for filename in glob.glob(ofile + '*.txt'):
	os.remove(filename)

# Atom types (6 - Na+, 4 - S, 9 - O in H2O)
#       1   12.01115  # C
#       2   18.99840  # F
#       3   15.99940  # GFOc
#       4   32.06400  # GGS
#       5   15.99940  # GOs      
#	    6    22.990 	 # Na
#       7   195.09  # P
#       8   1.00797  # WH
#       9   15.99940 # WO

typeID = range(1,10)
type_name = ['c', 'f', 'oc', 's', 'os', 'ion', 'pt', 'wh', 'wo']

out_file = {}
for i, name in zip(typeID, type_name):
	out_file[str(i)] = ofile + type_name[i-1] + '.txt'

# Number of bins
nbins = 50
# Direction 
idir = 'x'

# Times at which to compute
times = list(range(5153000,8152000,1000))
times.append(8152000)

#
# Compute number density as a function of x direction
#

for ID in typeID:
	densities = ctd.compute_spatial_density(dfile, ID, times, nbins, idir)
	io.write_data(out_file[str(ID)], densities, len(times))

