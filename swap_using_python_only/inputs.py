for z in range(0, 100):
    filename = 'input' + str(z) + '.in'
    file = open(filename, 'w')
    
    if z == 99:
        jump = ""
    else:
        jump = 'jump input' + str(z + 1) + '.in'
    file.write('# ---------- Initialize Simulation --------------------- \
    \nclear \
    \nunits metal \
    \ndimension 3 \
    \nboundary p p p \
    \natom_style atomic \
    \natom_modify map array\
    \n# ---------- Create Atoms --------------------- \
    \nread_data NiCo' + str(z) + '.data\
    \n# ---------- Define Interatomic Potential --------------------- \
    \npair_style eam/alloy \
    \npair_coeff * * NiCo-lammps-2014.alloy Ni Co\
    \nneighbor 2.0 bin  \
    \nneigh_modify delay 0 every 1 check yes page 100000 one 4000 \
    \n# ---------- Define Settings --------------------- \
    \ncompute eng all pe/atom \
    \n# EQUILIBRATION\
    \nreset_timestep	0\
    \ntimestep 0.005\
    \nvelocity all create 300 12345 mom yes rot no\
    \nfix 1 all npt temp 300 300 100 iso 0 0 1 \
    \n#set thermodynamic output\
    \nthermo 1000\
    \nthermo_style custom pe \
    \nrun 10000\
    \nunfix 1\
    \nprint "All done!" \
    \nprint "Energy = $(pe)" append info.txt\
    \n' + jump)
    file.close()

for z in range(0, 2):
    if z == 0:
        filename = 'first.in'
        jump = 'second.in'
        data = 'child1.data'
    if z == 1:
        filename = 'second.in'
        jump = ''
        data = 'child2.data'
    file = open(filename, 'w')
    file.write('# ---------- Initialize Simulation --------------------- \
    \nclear \
    \nunits metal \
    \ndimension 3 \
    \nboundary p p p \
    \natom_style atomic \
    \natom_modify map array\
    \n# ---------- Create Atoms --------------------- \
    \nread_data '+data+'\
    \n# ---------- Define Interatomic Potential --------------------- \
    \npair_style eam/alloy \
    \npair_coeff * * NiCo-lammps-2014.alloy Ni Co\
    \nneighbor 2.0 bin  \
    \nneigh_modify delay 0 every 1 check yes page 100000 one 4000 \
    \n# ---------- Define Settings --------------------- \
    \ncompute eng all pe/atom \
    \n# EQUILIBRATION\
    \nreset_timestep	0\
    \ntimestep 0.005\
    \nvelocity all create 300 12345 mom yes rot no\
    \nfix 1 all npt temp 300 300 100 iso 0 0 1 \
    \n#set thermodynamic output\
    \nthermo 1000\
    \nthermo_style custom pe \
    \nrun 10000\
    \nunfix 1\
    \nprint "All done!" \
    \nprint "Energy = $(pe)" append info.txt\
    \n' + jump)
    file.close()
