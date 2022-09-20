# -*- coding: utf-8 -*-
import numpy as np
import ase.io.vasp


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
    
def get_num_atoms(p="POSCAR"):
    f = open(p, "r")
    poscar = f.readlines()
    f.close()

    return(np.array(poscar[6].split(), dtype='int'))


node = input('Enter name of node\n')
npar = int(input('Enter NPAR (default 2)\n'))

node_core = {'regular12':12,'commom20':20,'regular20':20, 'regular24':24, 'low24':24, 'big24':24, 'new24':24, 'common32':32} 
nbands = int(sum(get_valence_electrons()*get_num_atoms())*0.6 + sum(get_num_atoms()))

nnodes = int(np.ceil(nbands/npar/node_core[node]))
print('NODES : {}\nNCORE = {}\nKPAR = {}'.format(nnodes, int(node_core[node]/npar), nnodes))
