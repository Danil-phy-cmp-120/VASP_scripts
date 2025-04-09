import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
import os

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = ['Al', 'Si', 'Ga_d', 'Ge_d', 'In_d', 'Sn']
magmoms = {'FIM1': '8*0.5 4*3.0 4*-0.1',
           'AFM_l': '8*0.5 2*3.0 2*-3.0 4*-0.1',
           'AFM_s': '8*0.5 3.0 -3.0 -3.0 3.0 4*-0.1'}


for i in range(len(z_atoms)):
    for key in magmoms.keys():
    
        if not os.path.isdir('ION/isif7/inverse/{}/{}'.format(z_atoms[i], key)):
            os.makedirs('ION/isif7/inverse/{}/{}'.format(z_atoms[i], key))

        cell = Structure(lattice=np.array([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]])*5.5,
                         species = ['Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Fe', 'Fe', 'Fe', 'Fe', z_atoms[i], z_atoms[i], z_atoms[i], z_atoms[i]],
                         
                         coords=[(0.25, 0.25, 0.25),
                                (0.75, 0.25, 0.75),
                                (0.25, 0.75, 0.75),
                                (0.75, 0.75, 0.25),
                                (0.5, 0.5, 0.0 ),
                                (0.0, 0.0, 0.0 ),
                                (0.5, 0.0, 0.5 ),
                                (0.0, 0.5, 0.5 ),
                                (0.5, 0.0, 0.0 ),
                                (0.0, 0.5, 0.0 ),
                                (0.0, 0.0, 0.5 ),
                                (0.5, 0.5, 0.5 ),
                                (0.25, 0.75, 0.25),
                                (0.75, 0.25, 0.25),
                                (0.25, 0.25, 0.75),
                                (0.75, 0.75, 0.75)], )

        poscar = inputs.Poscar(cell)             
                        
        potcar = inputs.Potcar(symbols=['Rh_pv', 'Fe', z_potentials[i]], functional="PBE_54")

        kpoints = inputs.Kpoints(kpts = [[7, 7, 7]])

        incar = inputs.Incar({'LCHARG' : True,
                               'LWAVE' :  True,
                               'ICHARG' : 1,
                               'IBRION' : 2,
                               'ISIF' : 7,
                               'POTIM' : 0.1,
                               'PREC' : 'Accurate',
                               'ENCUT' : 460,
                               'EDIFF' : 1e-7,
                               'EDIFFG' : 1e-3,
                               'ISMEAR' : 1,
                               'SIGMA' : 0.1,
                               'ISPIN' : 2,
                               'LMAXMIX' : 6,
                               'ADDGRID' : True,
                               'LORBIT' : 10,
                               'LASPH' : True,
                               'NSW ': 100,
                               'NEDOS' : 1000,
                               'NELM' : 101,
                               'ALGO' : 'Fast',
                               'NCORE' : 6,
                               'KPAR' : 2,
                               'MAGMOM': magmoms[key]})
                       
                       
        VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
        VaspInput.write_input('ION/isif7/inverse/{}/{}'.format(z_atoms[i], key))
