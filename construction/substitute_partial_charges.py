#
# Module for implementing correct partial charges in Nafion molecule
#

import math

class SubPCharges:
	''' Class with substitution process for every atom '''

	def __init__(self, atom_data, all_types, chfile):
		''' Sub charges in atom_data while depleting all_types'''

		self.data = atom_data
		self.nafion_types = all_types
		self.pcharge_conv = {}
		self.sum_of_charges = 0.0

		# Read partial charges file
		self.read_chfile(chfile)

	def read_chfile(self, chfile):
		with open(chfile, 'r') as fin:
			for line in fin:
				line = line.strip().split(',')
				self.pcharge_conv[line[0]] = float(line[1])
			
	def double_oxygen(self):
		pcharge = self.pcharge_conv['O=']
		IDs = []
		for key, value in self.nafion_types['5'].items():
			# Assumes IDs are continues from 1 to N_atoms
			self.data[int(key)-1] = (' ').join(value[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(value[-1]) + '\n'
			self.sum_of_charges += pcharge
			IDs.append(key)
		self.remove_IDs(self.nafion_types['5'], IDs)

		IDs = []
		# Atoms labelled cl are also o=-s
		for key, value in self.nafion_types['6'].items():
			# Assumes IDs are continues from 1 to N_atoms
			self.data[int(key)-1] = (' ').join(value[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(value[-1]) + '\n'
			self.sum_of_charges += pcharge
			IDs.append(key)
		self.remove_IDs(self.nafion_types['6'], IDs)
	
	def rest_of_the_chain(self):
		pcharge = self.pcharge_conv['S']

		IDs = []
		# For every sulfur
		for key, value in self.nafion_types['4'].items():
			
			# Substitute
			self.data[int(key)-1] = (' ').join(value[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(value[-1]) + '\n'
			self.sum_of_charges += pcharge
			# Store ID od S
			IDs.append(key)

			# Find and sub C5 and its Fs
			C5 = self.sub_closest_carbon(value, self.pcharge_conv['C5'])
			self.sub_closest_Fs(C5[0], self.pcharge_conv['F(C5)'], 2)
	
			# Find C4 + Fs then delete C5
			C4 = self.sub_closest_carbon(self.nafion_types['1'][C5[0]], self.pcharge_conv['C4'])
			self.sub_closest_Fs(C4[0], self.pcharge_conv['F(C4)'], 2)
			self.remove_IDs(self.nafion_types['1'], [C5[0]])
			
			# Find O2 then delete C4
			O2 = self.sub_closest_ether_oxygen(self.nafion_types['1'][C4[0]], self.pcharge_conv['O2']) 
			self.remove_IDs(self.nafion_types['1'], [C4[0]])

			# Find C2 + Fs then delete O2
			C2 = self.sub_closest_carbon(self.nafion_types['3'][O2[0]], self.pcharge_conv['C2'])
			self.sub_closest_Fs(C2[0], self.pcharge_conv['F(C2)'], 1)
			self.remove_IDs(self.nafion_types['3'], [O2[0]])
	
			# Subsitute the branch C3 and C1 by distinguishing one closest to O1
			O1 = self.sub_c1_c3_and_o1(self.nafion_types['1'][C2[0]])	

			# Substitute C0
			C0 = self.sub_closest_carbon(self.nafion_types['3'][O1[0]], self.pcharge_conv['C0'])
			self.sub_closest_Fs(C0[0], self.pcharge_conv['F(C0)'], 2)
			self.remove_IDs(self.nafion_types['3'], [O1[0]])
			self.remove_IDs(self.nafion_types['1'], [C0[0]])		

		self.remove_IDs(self.nafion_types['4'], IDs)
		
	def sub_c1_c3_and_o1(self, c2_values):
		# The branch - C1 and C3
		Cs = self.find_closest_atoms(c2_values[0][0], c2_values[-1], self.nafion_types['1'], 2)

		# Identify C1 - it is closer to O
		dCO1 = self.find_closest_atoms(Cs[0][0], self.nafion_types['1'][Cs[0][0]][-1], self.nafion_types['3'], 1) 
		dCO2 = self.find_closest_atoms(Cs[1][0], self.nafion_types['1'][Cs[1][0]][-1], self.nafion_types['3'], 1)
		if dCO1[0][0] < dCO2[0][0]:
			C1 = Cs[0][0]
			C3 = Cs[1][0]
		else:
			C1 = Cs[0][0]
			C3 = Cs[1][0]
		
		# Subsitute C3, its Fs and delete
		c3_val = self.nafion_types['1'][C3] 
		pcharge = self.pcharge_conv['C3']
		self.sum_of_charges += pcharge
		self.data[int(C3)-1] = (' ').join(c3_val[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(c3_val[-1]) + '\n'
		self.sub_closest_Fs(C3, self.pcharge_conv['F(C3)'], 3)			
		self.remove_IDs(self.nafion_types['1'], [C3])

		# Substitute C1 and its Fs
		c1_val = self.nafion_types['1'][C1] 
		pcharge = self.pcharge_conv['C1']
		self.sum_of_charges += pcharge
		self.data[int(C1)-1] = (' ').join(c1_val[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(c1_val[-1]) + '\n'
		self.sub_closest_Fs(C1, self.pcharge_conv['F(C1)'], 2)			

		# Substitute O1 and then remove C1
		O1 = self.sub_closest_ether_oxygen(self.nafion_types['1'][C1], self.pcharge_conv['O1'])
		self.remove_IDs(self.nafion_types['1'], [C1])
		return O1

	def sub_closest_carbon(self, value, pcharge):
			# Find closest C and subsitute
			C = self.find_closest_atoms(value[0][0], value[-1], self.nafion_types['1'], 1)[0]
			self.sum_of_charges += pcharge
			c_val = self.nafion_types['1'][C[0]] 
			self.data[int(C[0])-1] = (' ').join(c_val[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(c_val[-1]) + '\n'
			# Return for reference
			return C
	
	def sub_closest_Fs(self, C_ID, pcharge, n_tot):
			# Substitute F attached to C
			val = self.nafion_types['1'][C_ID]
			CFs = self.find_closest_atoms(C_ID, val[-1], self.nafion_types['2'], n_tot)
			for f in CFs:
				self.sum_of_charges += pcharge
				f_val = self.nafion_types['2'][f[0]] 
				self.data[int(f[0])-1] = (' ').join(f_val[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(f_val[-1]) + '\n'
				self.remove_IDs(self.nafion_types['2'], [f[0]])	

	def sub_closest_ether_oxygen(self, value, pcharge):
			# Find closest O and subsitute
			O = self.find_closest_atoms(value[0][0], value[-1], self.nafion_types['3'], 1)[0]
			self.sum_of_charges += pcharge
			o_val = self.nafion_types['3'][O[0]] 
			self.data[int(O[0])-1] = (' ').join(o_val[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(o_val[-1]) + '\n'
			# Return for reference
			return O

	def rest_of_nafion(self):
		''' Substitute partial charges in the remaining of Nafion '''
		# It should only be the backbone 
		
		# Cs
		IDs = []
		for key, value in self.nafion_types['1'].items():
			# Substitute
			pcharge = self.pcharge_conv['C(skel)']
			self.data[int(key)-1] = (' ').join(value[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(value[-1]) + '\n'
			self.sum_of_charges += pcharge
			# Store ID od S
			IDs.append(key)
		self.remove_IDs(self.nafion_types['1'], IDs)
		
		# Fs
		IDs = []
		for key, value in self.nafion_types['2'].items():
			# Substitute
			pcharge = self.pcharge_conv['F(skel)']
			self.data[int(key)-1] = (' ').join(value[0]) + ' ' + str(pcharge) + ' ' +  (' ').join(value[-1]) + '\n'
			self.sum_of_charges += pcharge
			# Store ID od S
			IDs.append(key)
		self.remove_IDs(self.nafion_types['2'], IDs)
		print(self.sum_of_charges)

	def find_closest_atoms(self, atom_ID, coords, atom_dict, n_tot):
		''' Return n_tot closes atoms to coords '''
		
		distances = []
		# Compute and store distances and IDs
		for key, atom in atom_dict.items():
			# Exclude self
			if key != atom_ID: 
				dx = float(coords[0]) - float(atom[-1][0])
				dy = float(coords[1]) - float(atom[-1][1])
				dz = float(coords[2]) - float(atom[-1][2])
			
				dist = math.sqrt(dx*dx + dy*dy + dz*dz)

				distances.append((key, dist))

		# Then sort and return n_tot smallest
		temp = sorted(distances, key=lambda x: x[1])
		return temp[:n_tot]

	def remove_IDs(self, del_dict, IDs):
		''' Remove all entries from del_dict with keys in IDS '''

		for del_id in IDs:
			del del_dict[del_id]

	def write_data_file(self, fname, begining, ending):
		with open(fname, 'w') as fout:
			for line in begining:
				fout.write(line)
			for line in self.data:
				fout.write(line)
			fout.write('\n')
			for line in ending:
				fout.write(line)

def load_data_file(dfile):
	''' Loads the data file with separated Atoms section '''

	begining = []
	ending = []
	all_atoms = []

	with open(dfile, 'r') as fin:

		# Read everything befor Atoms data
		for line in fin:
			begining.append(line)
			if 'Atoms' in line:
				line = next(fin)
				begining.append(line)
				break

		# Read atoms data
		for line in fin:
			all_atoms.append(line)
			if 'Bonds' in line:
				# Remove last two lines and append to 
				# trailing part
				ending.append(all_atoms.pop())
				all_atoms.pop()
				break
	
		# Read the remainder
		for line in fin:
			ending.append(line)

	return begining, all_atoms, ending

def split_atoms(data):
	''' Distribute all Nafion atoms into type dicts with IDs as keys '''
	# data - as loaded from data file
	# Hardcoded to Nafion

	# All the atom types in a Nafion molecule
	# and cl that is there as an O= substitute
	all_types = { '1' : {}, '2' : {},
				  '3' : {}, '4' : {},
				  '5' : {}, '6' : {}}
	
	for atom in data:
		temp = atom.strip().split()
		if int(temp[2]) < 7:
			all_types[temp[2]][temp[0]] = [temp[0:3], float(temp[3]), temp[4:]]
		
	return all_types

if __name__ == '__main__':

	#
	# Preprocessing
	#

	# Retrieve data file in pieces for easier processing
	begining, all_atoms, ending = load_data_file('nafion_pre_charge.data')
	# Select and group Nafion atoms 
	all_types = split_atoms(all_atoms)

	# 
	# Substitute partial charges
	#

	sub_charges = SubPCharges(all_atoms, all_types, 'nafion_charges.txt')
	sub_charges.double_oxygen()
	sub_charges.rest_of_the_chain()
	sub_charges.rest_of_nafion()
	
	#
	# Write new data file
	#

	sub_charges.write_data_file('nafion.data', begining, ending)



