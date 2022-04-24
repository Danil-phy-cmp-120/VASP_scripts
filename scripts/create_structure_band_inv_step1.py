import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = ['Al', 'Si', 'Ga_d', 'Ge_d', 'In_d', 'Sn']

for i in range(len(z_atoms)):
    
    if not os.path.isdir('band/inverse/step1/{}'.format(z_atoms[i])):
        os.makedirs('band/inverse/step1/{}'.format(z_atoms[i]))

    cell = Structure.from_file('{}/ION/inverse/{}/FIM1/CONTCAR'.format(os.getcwd(), z_atoms[i]))
    cell = SpacegroupAnalyzer(cell).find_primitive()

    poscar = inputs.Poscar(cell)
                        
    potcar = inputs.Potcar(symbols=['Rh_pv', 'Fe', z_potentials[i]], functional="PBE_54")

    kpoints = inputs.Kpoints(kpts = [[11, 11, 11]])

    incar = inputs.Incar({'LCHARG': True,
                          'LWAVE': True,
                          'ICHARG': 1,
                          'IBRION': 2,
                          'ISIF': 3,
                          'POTIM': 0.1,
                          'PREC': 'Accurate',
                          'ENCUT': 460,
                          'EDIFF': 1e-7,
                          'EDIFFG': 1e-3,
                          'ISMEAR': 1,
                          'SIGMA': 0.1,
                          'ISPIN': 2,
                          'LMAXMIX': 6,
                          'ADDGRID': True,
                          'LORBIT': 10,
                          'LASPH': True,
                          'NSW ': 0,
                          'NEDOS': 1000,
                          'NELM': 101,
                          'ALGO': 'Fast',
                          'NCORE': 12,
                          'KPAR': 1,
                          'MAGMOM': '2*0.5 3.0 -0.1'})

                       
                       
    VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
    VaspInput.write_input('band/inverse/step1/{}'.format(z_atoms[i]))