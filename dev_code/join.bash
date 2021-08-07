#!/bin/bash

#note that top files are assumed to be same. If they change
#then, this part needs to be modified.

rm energy-ga.txt

cp top.1 top

cat top conf.1 > first.data
cat top conf.2 > second.data


mpirun -np 2 lmp_mpi < first.in
