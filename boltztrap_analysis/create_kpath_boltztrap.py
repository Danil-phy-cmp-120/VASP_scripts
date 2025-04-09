import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.outputs import Outcar
import os


structure = Structure.from_file('CONTCAR')

#structure = SpacegroupAnalyzer(cell).get_primitive_standard_structure()

kpath = HighSymmKpath(structure)
kpoints = inputs.Kpoints().automatic_linemode(100, kpath)

label = kpoints.labels
kpts = kpoints.kpts

kpts_boltztrap = [kpts[0]]
label_boltztrap = [label[0]]
for k in range(2, len(kpts), 2):
    if label[k] == label[k-1]:
        kpts_boltztrap += [kpts[k]]
        label_boltztrap += [label[k]]
    else:
        kpts_boltztrap += [kpts[k-1]] 
        kpts_boltztrap += [kpts[k]] 
        label_boltztrap += [label[k-1]]
        label_boltztrap += [label[k]] 
kpts_boltztrap += [kpts[-1]]
label_boltztrap += [label[-1]]

kpts_str = ''
for i in kpts_boltztrap:
    kpts_str += '[{}, {}, {}], '.format(i[0], i[1], i[2])
print(kpts_str)
    
