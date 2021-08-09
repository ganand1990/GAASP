#Deletes the previously created files
import os
os.system("rm NiCo*.data")
os.system("rm info.txt")
os.system("rm child*.data")

#data file creation

from random import sample
import os
import numpy as np
from ase.build import bulk
from ase.io import write
import random

for i in range(0, 100):
    atom = bulk("Ni", "fcc", a = 3.6, cubic = True)*(5,5,5)
    pos = atom.get_positions()
    sym = atom.get_chemical_symbols()
    ran_no = sample(range(0,len(atom)), 100)
    for j in range(0, len(ran_no)):
        sym[ran_no[j]] = "Co"
    atom.set_chemical_symbols(sym)
    '''
    This script generates the data file for disordered alloys,
    which can be used in conjugation with lammps .in file.
    '''
    import numpy as np
    from ase.build import bulk
    from ase.io import write
    import random

    def number_det(atom,elements):
        '''

        '''
        num_types = len(elements)
        num_elem = len(atom)/num_types


        #checking if total number of elements id divisible by type of elements
        check_int = isinstance(num_elem,int)
        num_each = []
        if(check_int == False):
            near_int = round(num_elem) 
            for a in range(num_types):
                if(a == (num_types-1)):
                    extra_num = near_int + (len(atom) - (num_types * near_int))
                    # By above, we make sure that whole number of atoms are for each species
                    # which leads to an extra atom for last index 
                    num_each.append(extra_num)
                else:
                    num_each.append(near_int)
        if(check_int == True):
            for a in range(num_types):
                num_each.append(num_elem)

        return num_each

    def element_assign(atom,elements,num_each):
        '''
        This module calculates the list of atoms. 
        Such list contains the the name of chemical
        species in a sequence. 
        Input:
        atom: Atoms object
        elements: list containing chemical species
        num_each: list containing number of each
        chemical species.
        Returns:
        elem_list : List of atoms in a sequence, such 
        that, same elements are together. 
        '''
        listof_symbols = []
        elem_list = []
        for a,b in enumerate(num_each):
            for c in range(b):
                elem_list.append(elements[a])

        return elem_list

    def random_gen(atom):
        '''

        '''
        num_list = []
        for a in range(len(atom)):
            num_list.append(a)

    #     random.shuffle(num_list)
        positions = atom.get_positions()
        rand_pos = []
        for a in range(len(atom)):
            rand_pos.append(positions[num_list[a]])

        return rand_pos

    def elem_dict(elements):
        '''

        '''
        elem_id = {}
        for a,b in enumerate(elements):
            elem_id[b] = a + 1
        return elem_id

    #variables
    outpt_file = 'NiCo'+str(i)
    # cryst_struc = 'fcc'
    # latparam = 3.6
    elements = ["Ni", "Co"]
    # supercell_x = 7
    # supercell_y = 7
    # supercell_z = 7


    #generating atoms object
    # at = at3
    atom_sym = sym
    # list_pos_inis = []
    # for a in range(len(at)):
    #     list_pos_ini = np.zeros(3)
    #     list_pos_ini[0] = at.get_positions()[a][0] + at.get_cell()[0][0]/2.
    #     list_pos_ini[1] = at.get_positions()[a][1] + at.get_cell()[1][1]/2.
    #     list_pos_ini[2] = at.get_positions()[a][2] + at.get_cell()[2][2]/2.
    #     list_pos_inis.append(list_pos_ini)

    # at.set_positions(list_pos_inis)
    #at.set_chemical_symbols('Cu')
    #print (at.get_chemical_symbols()[2])

    #determining number of each chemical species
    num_each = number_det(atom,elements)

    #list of atoms in a sequence
    elem_list = element_assign(atom,elements,num_each)

    rand_pos = random_gen(atom)
    #print (rand_pos[0][0],rand_pos[0][1],rand_pos[0][2])

    elem_id = elem_dict(elements)
    #Writting LAMMPS data file

    with open('{}.data'.format(outpt_file),'w') as fdata:
        fdata.write('Alloy datafile containing disordered positions\n\n')

        #Header
        #specify number of atoms and types of atoms
        fdata.write('{} atoms\n'.format(len(atom)))
        fdata.write('{} atom types\n'.format(len(elements)))

        #specify the box dimensions 

        fdata.write('{} {} xlo xhi\n'.format(atom.get_cell()[1][0],atom.get_cell()[0][0]))
        fdata.write('{} {} ylo yhi\n'.format(0.0,atom.get_cell()[1][1]))
        fdata.write('{} {} zlo zhi\n'.format(0.0,atom.get_cell()[2][2]))
        fdata.write('\n')

        #Atoms section
        fdata.write('Atoms\n\n')

        #Write each position
        for i in range(len(atom)):

            if atom_sym[i] == 'Ni':
                atom_type = 1
            else:
                atom_type = 2
            
            fdata.write('{0} '.format(i+1))
            fdata.write('{0} '.format(atom_type))
            fdata.write('{0:.4f} '.format(rand_pos[i][0]))  
            fdata.write('{0:.4f} '.format(rand_pos[i][1]))
            fdata.write('{0:.4f} \n'.format(rand_pos[i][2]))


#system command
import os
os.system("/mnt/d/'LAMMPS 64-bit 15Apr2020'/bin/lmp_serial.exe -in input0.in")
os.system("rm input*.in")
os.system("rm NiCo*.data")

#parent finder

from random import sample
filename = "info.txt"
file = open(filename, "r")
x = []
y = []
for lines in file:
    line = lines.split(" = ")
    x.append(float(line[1]))
#print(x)
file.close()
file = open(filename, "r")
for lines in file:
    line = lines.split(" = ")
    y.append(float(line[1]))
file.close()
x.sort()
z = []
par = [0, 0]
for i in range(0, len(x[:10])):
    z.append(y.index(x[:10][i]))
#print(z)
r = sample(range(0, 1000), 2)
for j in range(0, 2):
    if r[j] < 400:
        par[j] = z[sample(range(0, 4), 1)[0]]
    elif r[j] >= 400 and r[j] < 700:
        par[j] = z[sample(range(4, 7), 1)[0]]
    elif r[j] >= 700 and r[j] < 900:
        par[j] = z[sample(range(7, 9), 1)[0]]
    else:
        par[j] = z[9]

#child maker

from lammps_data_cs import read_lammps_data
f1 = open('NiCo'+str(par[0])+'.data','r')
at1 = read_lammps_data(f1, style='atomic')
f1.close()
f2 = open('NiCo'+str(par[1])+'.data','r')
at2 = read_lammps_data(f2, style='atomic')
f2.close()

pos1 = at1.get_positions()
pos2 = at2.get_positions()
import random
ran = random.sample(range(0, 500), 2)
# print(pos1[ran[0]])
# print(pos2[ran[1]])
for i in range(0, len(pos1)):
    if np.array_equal(pos1[i], pos2[ran[1]]):
        break
# print(pos1[i])
for j in range(0, len(pos2)):
    if np.array_equal(pos2[j], pos1[ran[0]]):
        break
# print(pos2[j])

pos1[[ran[0], i],:] = pos1[[i, ran[0]],:]
pos2[[ran[1], j],:] = pos2[[j, ran[1]],:]
# print(pos1[ran[0]])
# print(pos1[i])
# print(pos2[ran[1]])
# print(pos2[j])

ls = at1.get_chemical_symbols()
for i in range(0, len(ls)):
    if ls[i] == 'H':
        ls[i] = 'Ni'
    elif ls[i] == 'He':
        ls[i] = 'Co'
at1.set_chemical_symbols(ls)
ls = at2.get_chemical_symbols()
for i in range(0, len(ls)):
    if ls[i] == 'H':
        ls[i] = 'Ni'
    elif ls[i] == 'He':
        ls[i] = 'Co'
at2.set_chemical_symbols(ls)

fdata1 = open('child1.data','w')

fdata1.write('Alloy datafile containing disordered positions\n\n')

#Header
#specify number of atoms and types of atoms
fdata1.write('{} atoms\n'.format(len(at1)))
fdata1.write('{} atom types\n'.format(len(elements)))

#specify the box dimensions 

fdata1.write('{} {} xlo xhi\n'.format(at1.get_cell()[1][0],at1.get_cell()[0][0]))
fdata1.write('{} {} ylo yhi\n'.format(0.0,at1.get_cell()[1][1]))
fdata1.write('{} {} zlo zhi\n'.format(0.0,at1.get_cell()[2][2]))
fdata1.write('\n')

#Atoms section
fdata1.write('Atoms\n\n')
atom_sym1 = at1.get_chemical_symbols()
#Write each position
for i in range(len(at1)):

    if atom_sym1[i] == 'Ni':
        atom_type = 1
    else:
        atom_type = 2
            
    fdata1.write('{0} '.format(i+1))
    fdata1.write('{0} '.format(atom_type))
    fdata1.write('{0:.4f} '.format(pos1[i][0]))  
    fdata1.write('{0:.4f} '.format(pos1[i][1]))
    fdata1.write('{0:.4f} \n'.format(pos1[i][2]))
fdata1.close()
    
fdata2 = open('child2.data','w')

fdata2.write('Alloy datafile containing disordered positions\n\n')

#Header
#specify number of atoms and types of atoms
fdata2.write('{} atoms\n'.format(len(at2)))
fdata2.write('{} atom types\n'.format(len(elements)))

#specify the box dimensions 

fdata2.write('{} {} xlo xhi\n'.format(at2.get_cell()[1][0],at2.get_cell()[0][0]))
fdata2.write('{} {} ylo yhi\n'.format(0.0,at2.get_cell()[1][1]))
fdata2.write('{} {} zlo zhi\n'.format(0.0,at2.get_cell()[2][2]))
fdata2.write('\n')

#Atoms section
fdata2.write('Atoms\n\n')
atom_sym2 = at2.get_chemical_symbols()
#Write each position
for i in range(len(at2)):

    if atom_sym2[i] == 'Ni':
        atom_type = 1
    else:
        atom_type = 2
            
    fdata2.write('{0} '.format(i+1))
    fdata2.write('{0} '.format(atom_type))
    fdata2.write('{0:.4f} '.format(pos2[i][0]))  
    fdata2.write('{0:.4f} '.format(pos2[i][1]))
    fdata2.write('{0:.4f} \n'.format(pos2[i][2]))
fdata2.close()

#run the child
os.system("/mnt/d/'LAMMPS 64-bit 15Apr2020'/bin/lmp_serial.exe -in first.in")
