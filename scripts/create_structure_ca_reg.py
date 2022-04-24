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
        for j in np.arange(0.8, 1.425, 0.025):
    
            if not os.path.isdir('ca/regular/{}/{}/{:.3f}'.format(z_atoms[i], key, j)):
                os.makedirs('ca/regular/{}/{}/{:.3f}'.format(z_atoms[i], key, j))

            structure = Structure.from_file('{}/ION/isif7/regular/{}/{}/CONTCAR'.format(os.getcwd(), z_atoms[i], key))
            
            D = np.array([[(j)**(-1.0/3.0), 0, 0],
                          [0, (j)**(-1.0/3.0), 0],
                          [0, 0, (j)**(2.0/3.0)]])

            structure.lattice = np.dot(structure.lattice.matrix, D)
            poscar = inputs.Poscar(structure)             
                        
            potcar = inputs.Potcar(symbols=['Rh_pv', 'Fe', z_potentials[i]], functional="PBE_54")

            kpoints = inputs.Kpoints(kpts = [[7, 7, 7]])

            incar = inputs.Incar({'LCHARG' : True,
                                   'LWAVE' :  True,
                                   'ICHARG' : 1,
                                   'IBRION' : 2,
                                   'ISIF' : 2,
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
            VaspInput.write_input('ca/regular/{}/{}/{:.3f}'.format(z_atoms[i], key, j))
