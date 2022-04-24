import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
import os

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = ['Al', 'Si', 'Ga_d', 'Ge_d', 'In_d', 'Sn']

lattice_parameter = np.loadtxt('ION/lattice_parameter.dat')

for i in range(len(z_atoms)):
    for j in np.arange(2, 12, 2): #давление в GPa
    
        if not os.path.isdir('pressure/regular/{}/{:.0f}'.format(z_atoms[i], j)):
            os.makedirs('pressure/regular/{}/{:.0f}'.format(z_atoms[i], j))

        cell = Structure(lattice=np.array([[lattice_parameter[i,0], 0, 0],
                                           [0, lattice_parameter[i,1], 0],
                                           [0, 0, lattice_parameter[i,2]]]),
                         species = ['Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Rh', 'Fe', 'Fe', 'Fe', 'Fe', z_atoms[i], z_atoms[i], z_atoms[i], z_atoms[i]],
                         coords=[(0.25, 0.25, 0.25),
                                (0.25, 0.75, 0.25),
                                (0.75, 0.25, 0.25),
                                (0.75, 0.75, 0.25),
                                (0.25, 0.25, 0.75),
                                (0.75, 0.25, 0.75),
                                (0.25, 0.75, 0.75),
                                (0.75, 0.75, 0.75),
                                (0.5,  0.0,  0.0),
                                (0.0,  0.5,  0.0),
                                (0.0,  0.0,  0.5),
                                (0.5,  0.5,  0.5),
                                (0.5,  0.5,  0.0),
                                (0.0,  0.0,  0.0),
                                (0.5,  0.0,  0.5),
                                (0.0,  0.5,  0.5)], )
        poscar = inputs.Poscar(cell)             
                        
        potcar = inputs.Potcar(symbols=['Rh_pv', 'Fe', z_potentials[i]], functional="PBE_54")

        kpoints = inputs.Kpoints(kpts = [[7, 7, 7]])

        incar = inputs.Incar({'LCHARG' : True,
                               'LWAVE' :  True,
                               'ICHARG' : 1,
                               'IBRION' : 2,
                               'ISIF' : 3,
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
                               'MAGMOM': '8*0.5 4*3.0 4*-0.1',
                               'PSTRESS': j*10})
                       
                       
        VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
        VaspInput.write_input('pressure/regular/{}/{:.0f}'.format(z_atoms[i], j))
