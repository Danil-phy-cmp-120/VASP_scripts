import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs 
from pymatgen.ext.matproj import MPRester
import os

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = {'Rh':'Rh_pv', 'Fe':'Fe', 'Al':'Al', 'Si':'Si', 'Ga':'Ga_d', 'Ge':'Ge_d', 'In':'In_d', 'Sn':'Sn'}


with MPRester("CLr5q0CgEJ5qvLlWI7") as m:

    for i in range(len(z_atoms)):
    
        for material in m.get_entries_in_chemsys(elements = ['Rh', 'Fe', z_atoms[i]]):
        
            structure = m.get_structure_by_material_id(material.entry_id)
            poscar = inputs.Poscar(structure) 
            
            paw = [z_potentials[s.symbol] for s in structure.types_of_species]                
            potcar = inputs.Potcar(symbols=paw, functional="PBE_54")

            kpoints = inputs.Kpoints.automatic_density(structure, 5000)

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
                               'NCORE' : 12,
                               'KPAR' : 2})

            VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
            if len(structure.species) > 60:
                if not os.path.isdir('hull_energy/{}/big/{}_{}'.format(z_atoms[i], material.composition.reduced_formula, material.entry_id)):
                    os.makedirs('hull_energy/{}/big/{}_{}'.format(z_atoms[i], material.composition.reduced_formula, material.entry_id))
                    VaspInput.write_input(
                        'hull_energy/{}/big/{}_{}'.format(z_atoms[i], material.composition.reduced_formula,
                                                      material.entry_id))
            else:
                if not os.path.isdir('hull_energy/{}/small/{}_{}'.format(z_atoms[i], material.composition.reduced_formula, material.entry_id)):
                    os.makedirs('hull_energy/{}/small/{}_{}'.format(z_atoms[i], material.composition.reduced_formula, material.entry_id))
                    VaspInput.write_input('hull_energy/{}/small/{}_{}'.format(z_atoms[i], material.composition.reduced_formula, material.entry_id))
