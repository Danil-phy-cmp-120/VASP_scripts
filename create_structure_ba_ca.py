import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
import os
import shutil

potentials = {'Fe':'Fe', 'Ge':'Ge_d'}

path_a0 = input('pathion\n')

for i in np.arange(0.85, 1.15, 0.05):
    for j in np.arange(0.85, 1.15, 0.05):
    
        if not os.path.isdir('ba_ca_{}/{:.3f}/{:.3f}'.format(path_a0, i, j)):
            os.makedirs('ba_ca_{}/{:.3f}/{:.3f}'.format(path_a0, i, j))

        structure = Structure.from_file('{}/CONTCAR'.format(path_a0))
           
        D = np.array([[(i*j)**(-1.0/3.0), 0, 0],
                      [0, j**(-1.0/3.0)*i**(2.0/3.0), 0],
                      [0, 0, i**(-1.0/3.0)*j**(2.0/3.0)]])

        structure.lattice = np.dot(structure.lattice.matrix, D)
        structure.sort()
        poscar = inputs.Poscar(structure)
                        
        paw = [potentials[s.symbol] for s in structure.types_of_species]
        potcar = inputs.Potcar(symbols=paw, functional="PBE_54")
                      
        kpoints = inputs.Kpoints.automatic_density(structure, 4000)

        incar = inputs.Incar({'LCHARG' : False,
                                   'LWAVE' :  False,
                                   'ICHARG' : 1,
                                   'IBRION' : 2,
                                   'ISIF' : 2,
                                   'POTIM' : 0.1,
                                   'PREC' : 'High',
                                   'ENCUT' : 410,
                                   'EDIFF' : 1e-8,
                                   'EDIFFG' : -1e-2,
                                   'ISMEAR' : 1,
                                   'SIGMA' : 0.1,
                                   'ISPIN' : 2,
                                   'LMAXMIX' : 6,
                                   'ADDGRID' : False,
                                   'LORBIT' : 10,
                                   'LASPH' : True,
                                   'NSW ': 100,
                                   'NEDOS' : 1000,
                                   'NELM' : 101,
                                   'ALGO' : 'Normal',
                                   'NCORE' : 12,
                                   '#KPAR' : 4,
                                   'MAGMOM': '14*3 2*0'})
                       
                       
        VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
        VaspInput.write_input('ba_ca_{}/{:.3f}/{:.3f}'.format(path_a0, i, j))
