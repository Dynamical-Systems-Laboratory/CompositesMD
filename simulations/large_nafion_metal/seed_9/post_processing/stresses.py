import sys
py_path = '../../../../postprocessing/io_operations/'
sys.path.insert(0, py_path)

import sys
py_path = '../../../../postprocessing/'
sys.path.insert(0, py_path)

import glob, os
import stress_processing as sp
import io_module as io

#
# Input
#

# Input file
dfile = '../nafion.d'
# Output file
ofile = 'stresses_and_eng_out_'

# Remove all old files because this code appends
for filename in glob.glob(ofile + '*.txt'):
	os.remove(filename)

# Stress directions
out_file = {'xx' : ofile + 'xx.txt', 'yy' : ofile + 'yy.txt', 'zz' : ofile + 'zz.txt',
				'xy': ofile + 'xy.txt', 'xz': ofile + 'xz.txt', 'yz': ofile + 'yz.txt',
 				'ke_xx': ofile + 'ke_xx.txt', 'ke_yy': ofile + 'ke_yy.txt', 
				'ke_zz': ofile + 'ke_zz.txt','virial_xx': ofile + 'virial_xx.txt',
				'virial_yy': ofile + 'virial_yy.txt', 'virial_zz': ofile + 'virial_zz.txt',
				'kinetic': ofile + 'kinetic_energy.txt', 'potential': ofile + 'potential_energy.txt'}


# Number of bins
nbins = 50
# Direction 
idir = 'x'

# Times at which to compute
times = list(range(5154000,6153000,1000))
times.append(6153000)

#
# Compute stresses as a function of x direction
#

stresses = sp.compute_spatial_stress_and_energy(dfile, times, nbins, idir)
for str_time in stresses:
	for str_dir, str_value in str_time.items():
		with open(out_file[str_dir], 'a') as fout:
			for val in str_value:
				fout.write(str(val))
				fout.write(' ')
			fout.write('\n')


