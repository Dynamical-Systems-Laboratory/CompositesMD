
#
# Module for converting COMPASS LAMMPS input into OPLS-AA LAMMPS input
#

import subprocess
import copy
import math

class DataFile:
	''' Creating and handling of LAMMPS .data file '''
		
	def __init__(self, orig_name, conv_name):
		# These need to be full paths
		self.compass = orig_name
		self.oplsaa = conv_name
		
	def write_all(self, atoms):
		''' Copy the original compass .data file while substituting Atoms data '''
		with open(self.compass, 'r') as fin:
			with open(self.oplsaa, 'w') as fout:
				for line in fin:
					if 'Atoms' in line:
						# This moves the fin and fout beyond Atoms
						fin, fout = self.sub_atoms(line, fin, fout, atoms)
					else:
						fout.write(line)
	def sub_Cl(self):
		''' Substitute Cl mass with O mass '''
		# Read the entire opls file
		# This has to be done after charge sub
		print('*'*10, ' DataFile.sub_Cl(): Have you substituted the charges first? ', '*'*10)

		with open(self.oplsaa, 'r') as fin:
			temp = fin.readlines()
		# Now write it back but substitute Cl mass
		# while doing so

		# Guard against bugs - only one sub right in 
		# the Masses section
		sub_cl = False
		with open(self.oplsaa, 'w') as fout:
			for line in temp:
				if 'Masses' in line:
					sub_cl = True
				if sub_cl and '35.45300' in line:
					line = line.strip().split()
					line[1] = '15.99940'
					line = ' '*7 + (' ').join(line) + '\n'
					sub_cl = False
				fout.write(line)


	def sub_atoms(self, line, fin, fout, atoms):
		''' Substitute atom data into fout while 
				incrementing the fin '''
		# Write Atoms and newline
		fout.write(line)
		next(fin)
		fout.write('\n')
		# Then write substituted values
		for atom in atoms:
			fout.write(atom)
			next(fin)
		return fin, fout
  
class ParamFile:
	''' Creating and handling of LAMMPS .param file '''

	# Right now it only just copies the potential file
	# But will be adding checking in the future

	def __init__(self, orig_name, new_name):
		self.orig = orig_name
		self.new = new_name

	def copy_orig(self):
		subprocess.run([(' ').join(['cp', self.orig, self.new])], shell=True, check=True)

class AddIon:
	''' Class for handling ion additions to the modeled setup '''
	
	# Note that speciec spec does not need to be defined in nafion OPLS-AA 
	# potentials files

	def __init__(self, file_name, species, charge, mass, radius):
		self.fname = file_name
		self.ion = species
		self.charge = charge
		self.mass = mass
		self.R = radius
		# This is something in between Cl- and O ion
		self.R_O = 1.40 
		self.buff = []

	def add_cation(self, nspec, natoms):
		''' Positions a cation of species next to -SO3- group '''

		# Load the whole file 
		with open(self.fname, 'r') as fin:
			orig = fin.readlines()
 
		# Find Mass section 
		for nln, line in enumerate(orig):
			if 'Masses' in line:
				m_ind = nln
				break
		
		# Collect until Masses end and add ion
		self.buff = copy.deepcopy(orig[0:m_ind + nspec + 2])
		self.buff.append(' '*7 + (' ').join([str(nspec+1), self.mass, '#', self.ion]) + '\n')
		self.buff.append('\n')
  
		# Find Atoms labelled with "# cl", record their information
		atom_flag = False
		cl_info = []
		s_info = []
		for nln, line in enumerate(orig):
			if 'Atoms' in line:
				atm_ind = nln
				atom_flag = True
			if atom_flag and '# cl' in line and nln > m_ind + nspec + 2:
				cl_info.append(line.strip().split())
			if atom_flag and '# s' in line and nln > m_ind + nspec + 2:
				s_info.append(line.strip().split())
			if 'Bonds' in line:
				bnd_ind = nln
				break

		# Simple version for now, sort based on 
		# x coordinate, assume closeness based on that
		cl_info = sorted(cl_info, key = lambda x: float(x[4]))	
		s_info = sorted(s_info, key = lambda x: float(x[4]))

		# Coppy Atoms information
		self.buff += orig[atm_ind:bnd_ind - 1]
		 
		# Now add cations
		cur_atoms = natoms
		for sulf_ind, atom in enumerate(cl_info):

			# Basic information
			cur_atoms += 1
			add_ion = ''
			add_ion = (' ').join([str(cur_atoms), atom[1], str(nspec+1), self.charge])
		  
			# For now move only in x direction
			# but determine the direction like in hydronium
			if float(s_info[sulf_ind][4]) > float(atom[4]):
				at_dir = -1.0
			else:
				at_dir = 1.0
		
			add_ion += (' ' + str(float(atom[4]) + at_dir*(self.R + self.R_O)) + ' ')
			add_ion += (' ').join(atom[5:7])
			add_ion += (' # ' + self.ion + '\n')

			# Add to buff
			self.buff.append(add_ion)

		# Copy the remaining entries
		self.buff += '\n'
		self.buff += orig[bnd_ind:]
		
		# Fix the header 
		# This assumes fixed position of header information
		if 'atoms' in self.buff[2]:
			temp = self.buff[2].strip().split()
			temp[0] = str(cur_atoms)
			self.buff[2] = '\t'+ ' ' + (' ').join(temp) + '\n'
		else:
			print('! !'*10 + 'atoms not on line 3')

		if 'atom types' in self.buff[8]:
			temp = self.buff[8].strip().split()
			temp[0] = str(int(temp[0]) + 1)
			self.buff[8] = '\t' + ' ' + (' ').join(temp) + '\n'
		else:
			print('! !'*10 + 'atom types not on line 9')

		return cur_atoms - natoms

	def write_all(self, fname):
		''' Write self.buff list to file fname '''
		with open(fname, 'w') as fout:
			for line in self.buff:
				fout.write(line)

class AddHydronium:
	''' Class for handling hydronium ion additions to the modeled setup '''
	
		# Based on current water model
		
	def __init__(self, file_name, types, ntypes):	
		# All properties have the following order
		# [H, O, H, H] if they are lists
		
		self.fname = file_name
		# List of type IDs
		self.types = types
		
		# Total number of types
		self.nspec = ntypes
		
		# List of radiuses
		self.R = [0.25, 0.6, 0.25, 0.25]
		# List of partial charges
		self.charges = [0.416, -0.248, 0.416, 0.416]
		
		self.OH_bond_len = 0.973
		# In radians
		self.beta_angle = 1.193805
		
		# This is the radius of atomic O 
		# I think having it too large for hydronium causes
		# some hydroniums to get detached and simulation explodes
		self.R_O = 0.8 
		self.buff = []

	def add_hydronium(self, nmolecules, natoms, nbonds, nangles):
		''' Positions a hydronium ion next to -SO3- group '''
	
		# nmolecules - number of molecules before introducing 
		#	hydronium
		# natoms - number of atoms before introducing hydronium
		
		# Load the whole file 
		with open(self.fname, 'r') as fin:
			orig = fin.readlines()
		
		# Find Mass section 
		for nln, line in enumerate(orig):
			if 'Masses' in line:
				m_ind = nln
				break
		# Find Atoms labelled with "# cl" and "# s", 
		# record their information
		cl_info = []
		s_info = []
		atom_flag = False
		for nln, line in enumerate(orig):
			if 'Atoms' in line:
				atm_ind = nln
				atom_flag = True
			if atom_flag and '# cl' in line and nln > m_ind + self.nspec + 2:
				cl_info.append(line.strip().split())
			if atom_flag and '# s' in line and nln > m_ind + self.nspec + 2:
				s_info.append(line.strip().split())
			if 'Bonds' in line:
				bnd_ind = nln
				atom_flag = False
			if 'Angles' in line:
				angle_ind = nln
			if 'Dihedrals' in line:
				dih_ind = nln
				break
		
		# Simple version for now, sort based on 
		# x coordinate, assume closeness based on that
		cl_info = sorted(cl_info, key = lambda x: float(x[4]))	
		s_info = sorted(s_info, key = lambda x: float(x[4]))
	
		# Copy Atoms information
		self.buff += orig[0:bnd_ind - 1]		   
		
		# Now add hydronium
		cur_atoms = natoms
		cur_molec = nmolecules

		# IDs of atoms in each hydronium molecule
		h3o_molecules = []

		for sulf_ind, atom in enumerate(cl_info):
			# Basic information
			# Each hydronium ion represents 4 atoms
			cur_molec += 1
			# Molecule building direction
			if float(s_info[sulf_ind][4]) > float(atom[4]):
				at_dir = -1.0
			else:
				at_dir = 1.0 
			# IDs of atoms in this molecule
			mol_str = []
			for i in range(4):
				cur_atoms += 1
				add_ion = ''
				add_ion = (' ').join([str(cur_atoms), str(cur_molec), str(self.types[i]), str(self.charges[i])])
				mol_str.append(str(cur_atoms))
				# For now move only in x direction
				# H
				if i == 0:
					add_ion += (' ' + str(float(atom[4]) + at_dir*(self.R[i] + self.R_O)) + ' ')
					i0 = 5
				# O
				elif i == 1:
					add_ion += (' ' + str(float(atom[4]) + at_dir*(self.R[i] + self.R_O + 2.0*self.R[0])) + ' ')
					i0 = 5
				# Side hydrogens
				else:
					base = float(atom[4]) + at_dir*(self.R_O + 2.0*(self.R[0] + self.R[1]))
					xpos = base + self.OH_bond_len*math.cos(self.beta_angle)
					if i == 3:
						ypos = float(atom[5]) + self.OH_bond_len*math.sin(self.beta_angle)
					else:
						ypos = float(atom[5]) - self.OH_bond_len*math.sin(self.beta_angle)	
					add_ion += (' ' + str(xpos) + ' ' + str(ypos) + ' ')
					i0 = 6
				
				add_ion += (' ').join(atom[i0:7]) 
				if i == 1:
					add_ion += ' # oh* \n'
				else:
					add_ion += ' # hhw \n'		
		
				# Add to buff
				self.buff.append(add_ion)
			h3o_molecules.append(mol_str)
	
		self.buff += '\n'
		# Copy bonds
		self.buff += orig[bnd_ind:angle_ind - 1]
		# And add new ones corresponding to H3O+
		cur_bonds = nbonds
		for mol in h3o_molecules:
			for i in [0,2,3]:
				cur_bonds += 1
				self.buff += (' ').join([str(cur_bonds), '8',  mol[i], mol[1], '# hhw,oh*']) + '\n'  	
		
		# Copy angles
		self.buff += '\n'
		self.buff += orig[angle_ind:dih_ind - 1]
		# And add new ones corresponding to H3O+
		cur_angles = nangles
		for mol in h3o_molecules:
			for ps in [(0,2),(0,3),(2,3)]:
				cur_angles += 1
				self.buff += (' ').join([str(cur_angles), '14', mol[ps[0]], mol[1], mol[ps[1]], '# hhw,oh*,hhw']) + '\n'

		# Copy the remaining entries
		self.buff += '\n'
		self.buff += orig[dih_ind:]
	
		# Fix the header 
		# This assumes fixed position of header information
		for tag, exp, num in zip(['atoms', 'bonds', 'angles'], [2, 3, 4], [cur_atoms, cur_bonds, cur_angles]):
			if tag in self.buff[exp]:
				temp = self.buff[exp].strip().split()
				temp[0] = str(num)
				self.buff[exp] = '\t'+ ' ' + (' ').join(temp) + '\n'
			else:
				print('! !'*10 + tag + ' not on line ' + str(exp + 1))		
		return cur_molec - nmolecules
	
	def write_all(self, fname):
		''' Write self.buff list to file fname '''
		with open(fname, 'w') as fout:
			for line in self.buff:
				fout.write(line)

class ChargeConv:
	''' Class for reading, storing, using, and manipulating 
			charge conversion rules '''

	def __init__(self, charge_conv):
		''' Initialize and collect conversion rules '''
		self.charges = {}
		with open(charge_conv, 'r') as fin:
			for line in fin:
				temp = line.strip().split()
				# Skip the comments
				if (not temp) or (temp[0] == '#'):
					continue
				# Add each conversion rule to the dictionary
				self.charges[temp[0]] = temp[1]

	def get_charge(self, atom_key):
		''' Retrieve partial charge that corresponds to 
				atom_key combination '''

		# Note: allows for constructed keys having different bond order
		# than keys loaded from the charge file; drawback is uniqueness
		# needs to be monitored
	   
		# It assumes that atoms in a bond are alphabetically ordered
		# i.e. they should be c,f not f,c

		# Split by hashes to obtain a list of bonded pairs
		atom_key_set = list(atom_key.split('#'))
		# This is due to specific format of the atom keys
		# as collected from a data file - removes the empty
		# string
		if '' in atom_key_set:
			atom_key_set.remove('')
		
		# Now make a list out of loaded tags
		# this is a tuple, second element stores the original key for reference
		tags_set_list = [(list(x.split('#')), x) for x in list(self.charges.keys())]

		# And for each see if it matches searched key
		for tag in tags_set_list:
			if sorted(atom_key_set) == sorted(tag[0]):
				# Match found
				F_key = tag[1] + '//c,f'
				if F_key in self.charges:
					return self.charges[tag[1]], True, self.charges[F_key]
				else:
					return self.charges[tag[1]], False, False

		# This will be returned either because of an error
		# or for the F atoms detected during C detection
		return False, False, False


class AtomsCharge:
	''' Class that reads, detects, and stores data on 
			each atoms partial charge '''

	def __init__(self, compass_data, conv_rules):
		''' Does all the detection and substitution '''
		# compass_data - lammps data file with compass potentials
		# conv_rules   - instance of ChargeConv class

		# Read the atoms and bonds info into a nested list, sublist per atom
		# All atom data are nested lists, one per line, split and stripped
		# Bonds are the whole line with blanks and newlines
		self.atoms_data = []
		self.bonds_data = []
		self.new_atoms_data = []
		self.fluorine_ids = {}

		with open(compass_data, 'r') as fin:
			lines = fin.readlines()
			# First find and store indices of lines
			# after which the data is located
			# This does assume the data is stored
			# as Atoms then Bonds then Angles
			idx = {}
			for i, line in enumerate(lines):
				if 'Atoms' in line:
					idx['Atoms'] = i 
				if 'Bonds' in line:
					idx['Bonds'] = i
				if 'Angles' in line:
					idx['Angles'] = i
			# Now collect each segment assuming
			# empty line between the data name and data
			self.atoms_data = lines[idx['Atoms'] + 2:idx['Bonds']]
			self.bonds_data = lines[idx['Bonds'] + 2:idx['Angles']]
			self.new_atoms_data = copy.deepcopy(self.atoms_data)

		# For each atom in the list find the bonds,
		# form a key, then use ChargeConv class to retrieve 
		# new partial charge
		for num, atom in enumerate(self.atoms_data):
			atom = atom.strip().split()
			self.atoms_data[num] = atom
			if not atom:
				continue
			atom_id = atom[0]
			atom_type = atom[-1]
			atom_key = ''
			for bond in self.bonds_data:
				bond = bond.strip().split()
				if not bond:
					continue
				if atom_id == bond[2] or atom_id == bond[3]:
					# Sort individual bond tags like c,f into
					# alphabetical order required by ChargeConv
					temp = sorted(bond[-1].split(','))
					atom_key += ( (',').join(temp) + '#')
					# If this is C and there is F atom
					# assumes that F only binds to C
					if 'c,f' in atom_key:
						if atom_id == 'f':
							continue
						if atom_id in self.fluorine_ids:
							 self.fluorine_ids[atom_id].append(int(bond[2]) if atom_id == bond[3] else int(bond[3]))
						else:
							self.fluorine_ids[atom_id] = [int(bond[2]) if atom_id == bond[3] else int(bond[3])]

			# Compare key and assign new charges
			# If there is a corresponding F atom, get_charge
			# will return it's charge too
			pcharge, fluorine, F_charge = conv_rules.get_charge(atom_key)
			
			# Substitute the atom
			# pcharge = False should signify c,f bond for type f
			if pcharge:
				atoms_split = self.new_atoms_data[num].strip().split()
				atoms_split[3] = pcharge
				self.new_atoms_data[num] = (' ').join(atoms_split) + '\n'
				if fluorine:
					f_ids = self.fluorine_ids[atom_id]
					# Substitute all F if any
					for fid in f_ids:
						fid -= 1
						atoms_split = self.new_atoms_data[fid].strip().split()
						if atoms_split:
							atoms_split[3] = F_charge
							self.new_atoms_data[fid] = (' ').join(atoms_split) + '\n'

	def write_all_new_atoms(self, fname):
		''' Save new atoms data to file '''
		with open(fname, 'w') as fout:
			for atom in self.new_atoms_data:
				fout.write(atom)

def charge_sum(atoms_data, split):
	''' Returns the sum of all atoms partial charges '''
	tot_charge = 0.0
	for atom in atoms_data:
		if not split:
			atom = atom.strip().split()
		if atom:
			tot_charge += float(atom[3])
	return tot_charge

# --------------------------------------------------

# Below is an example that demonstrates this module
# functionality for a complete run
# Paths are relative to the script as in github repository
# Run as ./compass_2_oplsaa.py 

if __name__ == '__main__':

	# Input EMC file 
	data_emc = 'constructed_setups/nafion_test_no_water/emc_files/emc_nafion.data'
	
	# Files with OPLS-AA potential information
	charges_conv = 'parameters/nafion_charges.txt'
	potential_conv = 'parameters/nafion_potential.txt'
	
	# Output files
	file_params = 'constructed_setups/nafion_test_no_water/nafion.params'
	file_data = 'constructed_setups/nafion_test_no_water/nafion.data'
	
	#
	# Convert charges
	#
	
	conv_rules = ChargeConv(charges_conv)
	atoms_pc = AtomsCharge(data_emc, conv_rules)
	
	#
	# Create .data file 
	#
	
	write_data = DataFile(data_emc, file_data)
	# Subsitute new charges
	write_data.write_all(atoms_pc.new_atoms_data)
	# Substitue Cl mass for O mass
	write_data.sub_Cl()

	#
	# Add Na+
	#
	
	ion_spec = 'na+'
	ion_charge = '+1'
	ion_mass = '22.990'
	# In Angstroms (current units)
	ion_radius = 1.02
	# Number of species without Na+ (based on Mass entries)
	nspec = 8
	# Number of atoms before adding ions
	natoms = 96

	add_na = AddIon(file_data, ion_spec, ion_charge, ion_mass, ion_radius)
	ncation = add_na.add_cation(nspec, natoms)
	add_na.write_all(file_data)

	#
	# Create .params file
	#
	
	write_params = ParamFile(potential_conv, file_params)
	write_params.copy_orig()

	#
	# Check charge
	#

	print('Total system charge: ', charge_sum(atoms_pc.new_atoms_data, False) + ncation*float(ion_charge))

