# ------------------------------------------------------------------
#
#	Tests for bin density computations module 
#
# ------------------------------------------------------------------

import sys
py_path = '../../postprocessing/io_operations/'
sys.path.insert(0, py_path)

py_path = '../../postprocessing/'
sys.path.insert(0, py_path)

py_path = '../common/'
sys.path.insert(0, py_path)

import math
from colors import *
import utils as utils
import asymmetries as asym
import io_module as io
import compare_files as comp

#
# Supporting functions
#

def check_sf(data, exp_sf):
	''' Tests if Sf factor is computed correctly '''
	
	sfactor = asym.compute_sf(data)

	if len(sfactor) != len(exp_sf):
		return False

	for exp_sft, sft in zip(exp_sf, sfactor):
		if not math.isclose(exp_sft, sft, rel_tol=1e-4):
			return False
	return True

def check_qf(data_neg, data_pos, exp_qf):
	''' Tests if Qf factor is computed correctly '''
	
	qfactor = asym.compute_qf(data_neg, data_pos)

	if len(qfactor) != len(exp_qf):
		return False

	for exp_qft, qft in zip(exp_qf, qfactor):
		if not math.isclose(exp_qft, qft, rel_tol=1e-4):
			return False
	return True

#
# Even number of bins 
# 

neg_bins = [[4, 4, 7, 7, 1, 1], [2, 6, 3, 8, 3, 2], [8, 2, 7, 5, 10, 3]]
pos_bins = [[9, 2, 2, 6, 7, 4], [5, 3, 9, 2, 4, 5], [10, 2, 6, 9, 6, 6]]

sf_neg_exp = [4.2426, 5.8310, 9.6437]
qf_exp = [9.1652, 12.6491, 6.4807]

utils.test_pass(check_sf(neg_bins, sf_neg_exp), 'Structural symmetry factor - even number of bins')
utils.test_pass(check_qf(neg_bins, pos_bins, qf_exp), 'Net charge symmetry factor - even number of bins')

#
# Odd number of bins 
# 

neg_bins = [[4, 4, 7, 7, 1], [2, 6, 3, 8, 3], [8, 2, 7, 5, 10]]
pos_bins = [[9, 2, 2, 6, 7], [5, 3, 9, 2, 4], [10, 2, 6, 9, 6]]

sf_neg_exp = [4.2426, 2.2361, 3.6056]
qf_exp = [1.4142, 3.6056, 7.2111]

utils.test_pass(check_sf(neg_bins, sf_neg_exp), 'Structural symmetry factor - odd number of bins')
utils.test_pass(check_qf(neg_bins, pos_bins, qf_exp), 'Net charge symmetry factor - odd number of bins')


