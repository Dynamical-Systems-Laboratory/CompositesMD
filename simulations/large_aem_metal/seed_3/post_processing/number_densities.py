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

# Input file
dfile = '../aem-metal.d'
# Output files
ofile = 'number_density_'

# Remove all old files because this code appends
for filename in glob.glob(ofile + '*.txt'):
	os.remove(filename)

# Atom types (6 - Na+, 4 - S, 9 - O in H2O)
#       1   12.01115  # C
#       2   35.453  # Cl
#       3   18.99840  # F
#       4   1.00797  # H
#       5   14.0067  # N
#       6   195.09  # P
#       7   1.00797  # WH
#       8   15.99940  # WO
#       9   12.01115  # c
#      10   1.00797  # ch

typeID = range(1,11)
type_name = ['c', 'ion', 'f', 'h', 'n', 'pt', 'wh', 'wo', 'cb', 'ch']

out_file = {}
for i, name in zip(typeID, type_name):
	out_file[str(i)] = ofile + type_name[i-1] + '.txt'

# Number of bins
nbins = 50
# Direction 
idir = 'x'

# Times at which to compute
times = list(range(5155000,8150000,10000))
times.append(8150000)

#
# Compute number density as a function of x direction
#

for ID in typeID:
	densities = ctd.compute_spatial_density(dfile, ID, times, nbins, idir)
	io.write_data(out_file[str(ID)], densities, len(times))

