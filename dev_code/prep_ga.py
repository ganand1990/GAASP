'''
This scrips prepares the lammps data file
do that it may be read by ga.c code.
'''
from ase import io
import numpy as np
import os
import shutil
import glob
from pathlib import Path
#TODO modify this part 
tot_num=30 #total number of random population
to_choose=20 #number of parents to be chosen
pop_name='CoNi' #common string in all data files 
ene_file='info.txt' #name of the energy file
top_data=10 #number of lines with top part in lammps data file
max_coord_width=10
#Generating parent population
num_val = sum(1 for line in open('{}'.format(ene_file)))

#TEST
if(num_val != tot_num):
    print ("total number of files and their energy \
            values don't match!\n")
    exit(0)

energies={}
with open('{}'.format(ene_file),'r') as inp:
    lines=inp.readlines()
    for num,line in enumerate(lines):
        energies[num]=line.split('\n')[0]
#print (energies)

#sorting the dictionary values in the ascending order.
sort_energies = sorted(energies.items(), key=lambda x: x[1],reverse=True)
#finding the indentities of the population from the sorted dictionary
filename='parent' #common string in all parent files
for a in range(to_choose):
    src='{}{}.data'.format(pop_name,sort_energies[a][0])
    dst='{}.{}'.format(filename,a)
    shutil.copy(src,dst)

#since most of the files are being appended
#old files are being removed here.

for f in glob.glob("atomlist*"):
    os.remove(f)
for f in glob.glob("top*"):
    os.remove(f)
for f in glob.glob("coord*"):
    os.remove(f)

for a in range(to_choose): #TODO change to to_choose
    with open('parent.{}'.format(a),'r') as inp:
        lines=inp.readlines()
        for num,line in enumerate(lines):
            if(num < top_data):
                with open('top.{}'.format(a),'a') as inp1:
                    inp1.write(line)
            else:
                with open('atomlist.{}'.format(a),'a') as inp2:
                    inp2.write(line.split(' ')[1])
                    inp2.write('\n')
                b=[]
                c=[]
                d=[]
                xcoord=''
                ycoord=''
                zcoord=''
                with open('coord.{}'.format(a),'a') as inp3:
                    b.append(line.split(' ')[2])
                    c.append(line.split(' ')[3])
                    d.append(line.split(' ')[4])
                    while len(xcoord) < max_coord_width:
                        req_pad=max_coord_width-len(b[0]) #required padding
                        for e in range(req_pad):
                            b.append('0')
                        xcoord = ''.join(b)
                    while len(ycoord) < max_coord_width:
                        req_pad=max_coord_width-len(c[0])
                        for e in range(req_pad):
                            c.append('0')
                        ycoord = ''.join(c)
                    while len(zcoord) < max_coord_width:
                        req_pad=max_coord_width-len(d[0])
                        for e in range(req_pad):
                            d.append('0')
                        zcoord = ''.join(d)
                    inp3.write(xcoord)
                    inp3.write(' ')
                    inp3.write(ycoord)
                    inp3.write(' ')
                    inp3.write(zcoord)
                    inp3.write('\n')
                    

#start_time=time.time()
#at = io.read('CoNi0.data',format='lammps-data',style='atomic')
#print("--- %s seconds ---" % (time.time() - start_time))
#for item in at:
#    print (item)

#finding number of parent population files

