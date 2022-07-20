import sys, os
py_path = '../postprocessing/io_operations/'
sys.path.insert(0, py_path)

import random, copy
import extract_data_file as ed

def new_box_save(data_file, new_data_file, data, box_dims):
	''' Save new data after adjusting box dimensions'''

	# Save into modified file
	found_atoms = False			
	with open(data_file, 'r') as fin:
		with open(new_data_file, 'w') as fout:
			for iline in fin:

				# Box dimensions
				if 'xlo' in iline:
					new_box = (' ').join([box_dims[0][0], box_dims[0][1], 'xlo', 'xhi'])
					fout.write(new_box + '\n')
					continue
				if 'ylo' in iline:
					new_box = (' ').join([box_dims[1][0], box_dims[1][1], 'ylo', 'yhi'])
					fout.write(new_box + '\n')
					continue				
				if 'zlo' in iline:
					new_box = (' ').join([box_dims[2][0], box_dims[2][1], 'zlo', 'zhi'])
					fout.write(new_box + '\n')
					continue				

				# Pre-Atoms portion - just write
				if ('Atoms' not in iline) and (found_atoms == False):
					fout.write(iline)
					continue

				# Atoms - write new data
				if 'Atoms' in iline:
					found_atoms = True
					fout.write('Atoms\n\n')
					iline = next(fin)
					for oline in data:
						fout.write(oline + '\n')
						iline = next(fin)
					continue

				# The remaining of data file - no change
				if found_atoms == True:
					fout.write(iline)  			
	
def shrink_cell(data_file, new_data_file):
	''' Make the cell slightly larger than the current chain '''
	
	with open(data_file, 'r') as fin:
		data = ed.extract_data_section(fin, 'Atoms')

	# Find current chain extent
	x_bounds = [1e5, -1e5]
	y_bounds = [1e5, -1e5]
	z_bounds = [1e5, -1e5]
	for line in data:
		entry = line.strip().split()
		if float(entry[2]) < x_bounds[0]:
			x_bounds[0] = float(entry[2])
		if float(entry[2]) > x_bounds[1]:
			x_bounds[1] = float(entry[2])	
		if float(entry[3]) < y_bounds[0]:
			y_bounds[0] = float(entry[3])
		if float(entry[3]) > y_bounds[1]:
			y_bounds[1] = float(entry[3])
		if float(entry[4]) < z_bounds[0]:
			z_bounds[0] = float(entry[4])
		if float(entry[4]) > z_bounds[1]:
			z_bounds[1] = float(entry[4])
	# Box dimensions increased by box_eps to accomodate
	# for atom radius
	box_eps = 1.0
	Lx =  (x_bounds[1]-x_bounds[0]) + 2.0*box_eps
	Ly =  (y_bounds[1]-y_bounds[0]) + 2.0*box_eps
	Lz =  (z_bounds[1]-z_bounds[0]) + 2.0*box_eps

	# Translate all atoms into the box
	for i, line in enumerate(data):
		entry = line.strip().split()
		entry[2] = str(float(entry[2])-x_bounds[0]+box_eps)
		entry[3] = str(float(entry[3])-y_bounds[0]+box_eps)
		entry[4] = str(float(entry[4])-z_bounds[0]+box_eps)
		data[i] = (' ').join(entry)

	# Save while adjusting box dimensions
	box_dims = [['0', str(Lx)], ['0', str(Ly)], ['0', str(Lz)]]
	new_box_save(data_file, new_data_file, data, box_dims)

def find_hydrogens(mol_ID, data):
	''' Return IDs of hydrogens in a water molecule identified
			by oxygen atom '''
	Hs = []
	for line in data:
		line = line.strip().split()
		if (line[1] == '7' and line[5] == mol_ID):
			Hs.append(line)
	if not (len(Hs) == 2):
		raise RuntimeError('There should be two hydrogens in a water molecule')
	return Hs

def shrink_cell_w_water(data_file, new_data_file):
	''' Make the cell larger than the current chain, bringing the water
			molecules into it '''
	
	with open(data_file, 'r') as fin:
		data = ed.extract_data_section(fin, 'Atoms')
		bond = ed.extract_data_section(fin, 'Bonds')

	# Find current chain extent
	x_bounds = [1e5, -1e5]
	y_bounds = [1e5, -1e5]
	z_bounds = [1e5, -1e5]
	for line in data:
		entry = line.strip().split()
		# Skip water
		if entry[1] == '7' or entry[1] == '8':
			continue
		if float(entry[2]) < x_bounds[0]:
			x_bounds[0] = float(entry[2])
		if float(entry[2]) > x_bounds[1]:
			x_bounds[1] = float(entry[2])	
		if float(entry[3]) < y_bounds[0]:
			y_bounds[0] = float(entry[3])
		if float(entry[3]) > y_bounds[1]:
			y_bounds[1] = float(entry[3])
		if float(entry[4]) < z_bounds[0]:
			z_bounds[0] = float(entry[4])
		if float(entry[4]) > z_bounds[1]:
			z_bounds[1] = float(entry[4])
	# Box dimensions increased by box_eps to accomodate
	# for atom radius
	box_eps = 10.0
	# Water section
	#dx_w = 20
	#x_water = [x_bounds[1], x_bounds[1]+dx_w]
	dy_w = 20
	y_water = [y_bounds[1], y_bounds[1]+dy_w]
	
	# Translate all Nafion/ion atoms into the polymer
	# box and all the water atoms to the side box 
	orig_data = copy.deepcopy(data)
	for i, line in enumerate(data):
		entry = line.strip().split()
		
		# Translate hydrogen(7) based on oxygen(8)
		if entry[1] == '7' or entry[1] == '8':
			if entry[1] == '8':
				Hs = find_hydrogens(entry[5], orig_data)
				dH1 = [float(entry[2]) - float(Hs[0][2]), 
								float(entry[3]) - float(Hs[0][3]), 
								float(entry[4]) - float(Hs[0][4])] 
				dH2 =  [float(entry[2]) - float(Hs[1][2]), 
								float(entry[3]) - float(Hs[1][3]), 
								float(entry[4]) - float(Hs[1][4])]
				# Translate oxygen, position hydrogens
				# based on the original distance
#				xi = random.uniform(x_water[0], x_water[1])
#				yi = random.uniform(y_bounds[0], y_bounds[1])
				xi = random.uniform(x_bounds[0], x_bounds[1])
				yi = random.uniform(y_water[0], y_water[1])				
				zi = random.uniform(z_bounds[0], z_bounds[1])
				entry[2], entry[3], entry[4] = str(xi), str(yi), str(zi)
				orig_data[i] = (' ').join(entry)
				Hs[0][2], Hs[0][3], Hs[0][4] = str(xi+dH1[0]), str(yi+dH1[1]), str(zi+dH1[2])
				orig_data[int(Hs[0][0])-1] = (' ').join(Hs[0])
				Hs[1][2], Hs[1][3], Hs[1][4] = str(xi+dH2[0]), str(yi+dH2[1]), str(zi+dH2[2])
				orig_data[int(Hs[1][0])-1] = (' ').join(Hs[1])
			else:
				continue
		else:
			entry[2] = str(float(entry[2])-x_bounds[0]+box_eps)
			entry[3] = str(float(entry[3])-y_bounds[0]+box_eps)
			entry[4] = str(float(entry[4])-z_bounds[0]+box_eps)
			orig_data[i] = (' ').join(entry)
	
	data = copy.deepcopy(orig_data)

	# Compute the box size
	#x_bounds[1] = x_water[1]
	y_bounds[1] = y_water[1]
	Lx =  (x_bounds[1]-x_bounds[0]) + 2.0*box_eps
	Ly =  (y_bounds[1]-y_bounds[0]) + 2.0*box_eps
	Lz =  (z_bounds[1]-z_bounds[0]) + 2.0*box_eps

	# Save while adjusting box dimensions
	box_dims = [['0', str(Lx)], ['0', str(Ly)], ['0', str(Lz)]]
	new_box_save(data_file, new_data_file, data, box_dims)

in_file = 'dpd.data'
out_file = 'single_chain.data'

shrink_cell_w_water(in_file, out_file)

