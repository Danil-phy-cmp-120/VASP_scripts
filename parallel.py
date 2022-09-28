# -*- coding: utf-8 -*-
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen.io.vasp.inputs as inputs
from pymatgen.core.structure import Structure
import os
from scipy.constants import pi


def get_valence_electrons(p="POTCAR"):
    f = open(p, "r")
    potcar = f.readlines()
    f.close()

    v_electrons = []
    for i in range(len(potcar)):
        inp = potcar[i].split()
        if len(inp) > 2 and inp[0] == 'PAW_PBE':
            v_electrons += [float(potcar[i+1])]
    return(np.array(v_electrons, dtype='int'))
    
def get_PREC_coefficient(p="INCAR"):
    PREC_coefficient = {'Normal':[3,2], 'Single':[3,1], 'Accurate':[4,2], 'Low':[3,2], 'Medium':[3,2], 'High':[4,2]}

    f = open(p, "r")
    incar = f.readlines()
    f.close()

    n = 0
    for i in incar:
        if i.find('PREC') != -1:
            n = i.split('=')[-1].replace(' ','').replace('\n','')
            
    if n != 0:   
        return(PREC_coefficient[n])
    else:
        return(PREC_coefficient['Normal'])    

def get_ENCUT(p="INCAR", l="POTCAR"):
    f = open(p, "r")
    incar = f.readlines()
    f.close()
    
    f = open(l, "r")
    potcar = f.readlines()
    f.close()

    ENCUT = 0
    for i in incar:
        if i.find('ENCUT') != -1:
            ENCUT = float(i.split('=')[-1])
    for i in potcar:
        if i.find('ENMAX') != -1:
            if float(i.split()[2].replace(';','')) > ENCUT:
                ENCUT = float(i.split('=')[1].replace(';',''))
  
    return(ENCUT)

    
def get_num_atoms(p="POSCAR"):
    f = open(p, "r")
    poscar = f.readlines()
    f.close()

    return(np.array(poscar[6].split(), dtype='int'))


structure = Structure.from_file('POSCAR')
poscar = inputs.Poscar(structure)
kpoints = inputs.Kpoints.from_file('KPOINTS').kpts

rytoev = 13.605826
aB = 0.529177249
cutof = np.array([np.sqrt(get_ENCUT()/rytoev) / (2*pi/(anorm/aB)) for anorm in poscar.structure.lattice.abc])
ng = get_PREC_coefficient()[0] * cutof + 0.5

nbands = int(sum(get_valence_electrons()*get_num_atoms())*0.6 + sum(get_num_atoms()))
nkdim = len(SpacegroupAnalyzer(structure).get_ir_reciprocal_mesh(kpoints))
        
memory = np.prod(ng) * nbands * nkdim * 16 + 4*((ng[0]*get_PREC_coefficient()[1])/2 + 1) * ng[1]*get_PREC_coefficient()[1] * ng[2]*get_PREC_coefficient()[1] *16
memory /= 1024**3


node = input('Enter name of node\n')
npar = int(input('Enter NPAR (default 2)\n'))

node_core = {'regular12':[12, 12, 11.5],'commom20':[20, 4, 63.0],'regular20':[20, 8, 63.0], 'low24':[24, 4, 62.9], 'big24':[24, 12, 126.0], 'new24':[24, 8, 188.0], 'common32':[32, 10, 188.0]} 

nnodes = int(np.ceil(nbands/npar/node_core[node][0]))

if nnodes < node_core[node][1]:
    print('NODES : {}\nNCORE = {}\nKPAR = {}'.format(nnodes, int(node_core[node][0]/npar), nnodes))
else:
    print('Not enough nodes to efficiently parallelize the task\n\t* select a larger NPAR\n\t* use other nodes\n\t* use non-paralleling settings (remove NPAR, NCORE, KPAR tags from INCAR)')
if memory < nnodes * node_core[node][2]:
    print('Required Memory = {:.3f} Gb'.format(memory))
else:
    print('Required Memory = {:.3f} Gb. Not enough memory on the nodes,\n increase the number of nodes to {} and use non-paralleling settings (remove NPAR, NCORE, KPAR tags from INCAR)'.format(memory, int(memory/node_core[node][2])))

