# ------------------------------------------------------------------
#
#	Tests for bin density file IO module 
#
# ------------------------------------------------------------------

import sys
py_path = '../../postprocessing/io_operations/'
sys.path.insert(0, py_path)

py_path = '../../postprocessing/'
sys.path.insert(0, py_path)

py_path = '../common/'
sys.path.insert(0, py_path)

from colors import *
import utils as utils
import load_bin_densities as lbd
import io_module as io
import compare_files as comp

#
# Test
#

# Inpute file
dfile_in = 'test_data/density_s.out'
# Files with expected data
exp_all_times_file = 'test_data/exp_density_s.txt'
exp_some_times_file = 'test_data/exp_select_density_s.txt'
# Output files
file_all_t = 'test_data/out_density_s.txt'
file_some_t = 'test_data/out_select_density_s.txt'

# Number of header lines at the top of the file
nhd = 3
# Number of bins
nbins = 20
# Number of time steps
n_all_t = 4
# Times to consider
times = [5000, 10145000, 10150000] 

# All time steps
data = lbd.load_bin_density(dfile_in, nhd)
io.write_rows(file_all_t, data, n_all_t)
utils.test_pass(comp.equal_files(exp_all_times_file, file_all_t), 'Binned density, all time steps')

# Select time steps
data = lbd.load_bin_density(dfile_in, nhd, times, False)
io.write_rows(file_some_t, data, len(times))
utils.test_pass(comp.equal_files(exp_some_times_file, file_some_t), 'Binned density, select time steps')
