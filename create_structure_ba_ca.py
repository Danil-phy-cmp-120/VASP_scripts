import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
import os
import shutil

potentials = {'In':'In_d', 'Mn':'Mn_pv', 'Ni':'Ni_pv'}

path_a0 = input('pathion\n')

for i in np.arange(0.85, 1.15, 0.05):
    for j in np.arange(0.85, 1.15, 0.05):
    
        if not os.path.isdir('ba_ca/{:.3f}/{:.3f}'.format(i, j)):
            os.makedirs('ba_ca/{:.3f}/{:.3f}'.format(i, j))

        structure = Structure.from_file('{}/CONTCAR'.format(path_a0))
           
        D = np.array([[(i*j)**(-1.0/3.0), 0, 0],
                      [0, j**(-1.0/3.0)*i**(2.0/3.0), 0],
                      [0, 0, i**(-1.0/3.0)*j**(2.0/3.0)]])

        structure.lattice = np.dot(structure.lattice.matrix, D)
        structure.sort()
        poscar = inputs.Poscar(structure)
                        
        paw = [potentials[s.symbol] for s in structure.types_of_species]
        potcar = inputs.Potcar(symbols=paw, functional="PBE_54")
                      
        kpoints = inputs.Kpoints.automatic_density(structure, 5000)

        incar = inputs.Incar({'LCHARG' : False,
                                   'LWAVE' :  False,
                                   'ICHARG' : 1,
                                   'IBRION' : 2,
                                   'ISIF' : 7,
                                   'POTIM' : 0.1,
                                   'PREC' : 'High',
                                   'EDIFF' : 1e-8,
                                   'EDIFFG' : 1e-3,
                                   'ISMEAR' : 1,
                                   'SIGMA' : 0.1,
                                   'ISPIN' : 2,
                                   'LMAXMIX' : 6,
                                   'ADDGRID' : False,
                                   'LORBIT' : 10,
                                   'LASPH' : False,
                                   'NSW ': 100,
                                   'NEDOS' : 1000,
                                   'NELM' : 101,
                                   'ALGO' : 'Normal',
                                   'NCORE' : 12,
                                   '#KPAR' : 4,
                                   'MAGMOM': '16*0.4 8*3.5 4*-3.5 4*0.04'})
                       
                       
        VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
        VaspInput.write_input('ba_ca/{:.3f}/{:.3f}'.format(i, j))
