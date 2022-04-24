import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
import os
from pymatgen.core.surface import SlabGenerator


z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = {'Rh':'Rh_pv', 'Fe':'Fe', 'Al':'Al', 'Si':'Si', 'Ga':'Ga_d', 'Ge':'Ge_d', 'In':'In_d', 'Sn':'Sn'}

structure = Structure.from_file('ION/isif3/regular/Al/FIM1/CONTCAR')

for size in range(6, 7):

    slabgen = SlabGenerator(structure, (0,0,1), size*structure.lattice.c, 10, center_slab=True)
    slabs = slabgen.get_slabs(symmetrize=True)
    print("Number of slabs:", len(slabs))
    
    for i in range(len(slabs)):
    
        if not os.path.isdir('slab/test_distance/regular/Al/{:.0f}/{}'.format(size, i)):
            os.makedirs('slab/test_distance/regular/Al/{:.0f}/{}'.format(size, i))
        
        poscar = inputs.Poscar(slabs[i])
        
        paw = [z_potentials[s.symbol] for s in slabs[i].types_of_species]                
        potcar = inputs.Potcar(symbols=paw, functional="PBE_54")

        kpoints = inputs.Kpoints.automatic_density(slabs[i], 5000)

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
                               'KPAR' : 2})
                       
                       
        VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
        VaspInput.write_input('slab/test_distance/regular/Al/{:.0f}/{}'.format(size, i))
