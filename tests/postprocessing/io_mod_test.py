#!/usr/local/bin/python3.6

# ------------------------------------------------------------------
#
#	Tests for IO module 
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
import io_module as io
import utils as utils

#
# Supporting functions
#

def correct_array_io(test_arrays, ncol, file_out):
	''' Check if written array equal to expected '''

	for i, data in enumerate(test_arrays):
		
		array_file = file_out + str(i) + '.txt'
		io.write_data(array_file, data, ncol[i])
		
		with open(array_file, 'r') as fin:
			wrote_array = fin.readlines();
	
		# Convert to int and do a direct comparison

		# Number of rows
		if not (len(wrote_array) == len(data[0])):
			return False

		# Content
		for iwa, line in enumerate(wrote_array):
			line = line.strip().split()
			# Check number of columns
			if not(len(line) ==	ncol[i]):
				return False
			# If empty line
			if not line:
				return False
			for jwa, entry in enumerate(line):
				if not (int(entry) == data[jwa][iwa]):
					return False
	return True

def correct_list_io(test_list, test_x, file_out):
	''' Check if written list equal to expected 
			test_list - 1D list to write 
			test_x - 1D list of independent variable 
						values to write 
			file_out - file name for writing '''

	# Write the list and independent variable
	io.write_list(file_out, test_x, test_list) 
	line_no = 0
	with open(file_out, 'r') as fin:
		for line, var, value in zip(fin, test_x, test_list):
			temp = line.strip().split()
			if not math.isclose(float(temp[0]), var):
				return False
			if not math.isclose(float(temp[1]), value):
				return False
			line_no += 1
	if line_no != len(test_x):
		return False
	return True


#
# Tests
#

#
# Saving ND arrays
#

# File to write to
file_out = './test_data/io_mod_test_out_'

# Contents (each sub-array is a specific test)
# Uses ints and a direct comparison
test_arrays = [[[1, 2, 3], [2, 4, 6]], [[4, 5, 6], [7, 8, 9], [10, 11, 12]]]
test_arrays.append([[1,2], [3,4], [5,6], [7,8]])
ncol = [2,3,4]

# Write, load, and compare to expected
utils.test_pass(correct_array_io(test_arrays, ncol, file_out), 'Writing N-dimensional arrays')

#
# Saving 1D lists
#

fout_1D = './test_data/write_1d_list.txt'
x_var = [0.1, 0.2, 0.9]
values = [1.3, 100.5, 1.1e-6]

utils.test_pass(correct_list_io(values, x_var, fout_1D), 'Writing 1D arrays with independent variable')


