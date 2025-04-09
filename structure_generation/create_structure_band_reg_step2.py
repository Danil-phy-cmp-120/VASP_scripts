import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = ['Al', 'Si', 'Ga_d', 'Ge_d', 'In_d', 'Sn']

for i in range(len(z_atoms)):

    if not os.path.isdir('band/regular/step2/{}'.format(z_atoms[i])):
        os.makedirs('band/regular/step2/{}'.format(z_atoms[i]))

    os.link('band/regular/step1/{}/CHGCAR'.format(z_atoms[i]), 'band/regular/step2/{}/CHGCAR'.format(z_atoms[i]))

    cell = Structure.from_file('band/regular/step1/{}/POSCAR'.format(z_atoms[i]))
    poscar = inputs.Poscar(cell)

    potcar = inputs.Potcar(symbols=['Rh_pv', 'Fe', z_potentials[i]], functional="PBE_54")

    kpath = HighSymmKpath(cell)
    kpoints = inputs.Kpoints().automatic_linemode(100, kpath)

    incar = inputs.Incar({'LCHARG' : False,
                          'LWAVE' :  False,
                          'ICHARG' : 11,
                          'IBRION' : 2,
                          'ISIF' : 3,
                          'POTIM' : 0.1,
                          'PREC' : 'Accurate',
                          'ENCUT' : 460,
                          'EDIFF' : 1e-7,
                          'EDIFFG' : 1e-3,
                          'ISMEAR' : 1,
                          'SIGMA' : 0.01,
                          'ISPIN' : 2,
                          'LMAXMIX' : 6,
                          'ADDGRID' : True,
                          'LORBIT' : 11,
                          'LASPH' : True,
                          'NSW ': 0,
                          'NEDOS' : 1000,
                          'NELM' : 101,
                          'ALGO' : 'Fast',
                          'NCORE' : 6,
                          'KPAR' : 2,
                          'MAGMOM': '2*0.5 3.0 -0.1'})
                       
                       
    VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
    VaspInput.write_input('band/regular/step2/{}'.format(z_atoms[i]))
