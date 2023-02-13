import numpy as np
import itertools
import os
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from more_itertools import distinct_permutations

potentials = {'Mn':'Mn_pv', 'V':'V_sv', 'Sb':'Sb', 'Ni':'Ni','Cr':'Cr_pv','Sc':'Sc_sv','As':'As','Ti':'Ti_sv','Fe':'Fe','Ru':'Ru_pv','Co':'Co', 'Ge':'Ge_d', 'Rh':'Rh_pv', 'Si':'Si', 'Ni':'Ni','P':'P', 'Ga':'Ga_d','Sn':'Sn_d'}

structure_0 = Structure.from_file('POSCAR')
position_for_permutation = [8, 9, 10, 11, 12, 13, 14, 15]
magmom = [0,0,0,0,0,0,0,0,1, 2, 3, 4, 5, 6, 7, 8]

species_for_permutation = []
species_constant = []
for i in range(len(structure_0.species)):
    if i in position_for_permutation:
        species_for_permutation += [structure_0.species[i]]
    else:
        species_constant += [structure_0.species[i]]

sg = SpacegroupAnalyzer(structure_0)
operations = sg.get_space_group_operations()

unique_structures = []
for n, p in enumerate(distinct_permutations(species_for_permutation)):
        print(p)
            
    species_new = species_constant + list(p)
    structure = Structure(structure_0.lattice, species_new, structure_0.frac_coords)

    flag = False
    for us in unique_structures:
        if operations.are_symmetrically_equivalent(structure, us):
            flag = True
            break
    if flag == False:
        unique_structures += [structure]


for n, s in enumerate(unique_structures):
    if not os.path.isdir('unique_structures'):
        os.makedirs('unique_structures')
    poscar = inputs.Poscar(s, sort_structure = True)

    paw = [potentials[i.symbol] for i in s.types_of_species]
    potcar = inputs.Potcar(symbols=paw, functional="PBE_54")

    kpoints = inputs.Kpoints.automatic_density(s, 5000)

    incar = inputs.Incar({#'METAGGA': 'SCAN',
                          'LCHARG': False,
                          'LWAVE': False,
                          'ICHARG': 1,
                          'IBRION': 2,
                          'ISIF': 3,
                          'POTIM': 0.005,
                          'PREC': 'Accurate',
                          'EDIFF': 1e-7,
                          'EDIFFG': 1e-3,
                          'ISMEAR': 1,
                          'SIGMA': 0.1,
                          'ISPIN': 2,
                          'LMAXMIX': 6,
                          'ADDGRID': True,
                          'LORBIT': 11,
                          'LASPH': True,
                          'NSW ': 100,
                          'NEDOS': 1000,
                          'NELM': 101,
                          'ALGO': 'Fast',
                          'NCORE': 6,
                          'KPAR': 2,
                          'MAGMOM': })

    VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
    VaspInput.write_input('unique_structures/{}'.format(n))