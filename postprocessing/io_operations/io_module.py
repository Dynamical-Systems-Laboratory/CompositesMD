# ------------------------------------------------------------------
#
#	Module for general I/O operations 				  
#
# ------------------------------------------------------------------ 

def write_data(fname, data, num_col):
	''' Save data with a fixed number of columns to file '''

	#
	# fname - file name
	# data - data set in a form of nested lists, one list per column
	# num_col - number of columns in the data
	# 

	with open(fname, 'w') as fout:
		for i in range(len(data[0])):
			for col in range(num_col):
				fout.write(str(data[col][i]))
				fout.write(' ')
			fout.write('\n')

def write_rows(fname, data, num_rows):
	''' Save data with a fixed number of rows to file '''

	#
	# fname - file name
	# data - data set in a form of nested lists, one list per row 
	#			row can be arbitrary length
	# num_rows - number of rows in the data
	# 

	with open(fname, 'w') as fout:
		for i in range(num_rows):
			temp = [str(x) for x in data[i]]
			fout.write((' ').join(temp))
			fout.write('\n')

def write_list(fname, x, data):
	''' Save a list data and corresponding independent 
			variable to file fname '''

	with open(fname, 'w') as fout:
		for var, d in zip(x, data):
			fout.write((' ').join([str(var), str(d)]))
			fout.write('\n')	
