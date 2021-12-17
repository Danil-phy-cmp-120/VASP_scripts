# -*- coding: utf-8 -*-
import numpy as np
import ase.io.vasp
import ase.spacegroup
from ase.visualize.plot import plot_atoms
import shutil
import os

KPPRA = 1000

cell = ase.io.vasp.read_vasp("POSCAR")
KP = KPPRA/len(cell)

reciprocal_cell_lengh = np.array([np.linalg.norm(np.linalg.inv(cell.cell)[0]), np.linalg.norm(np.linalg.inv(cell.cell)[1]), np.linalg.norm(np.linalg.inv(cell.cell)[2])])

reciprocal_cell_lengh_norm = reciprocal_cell_lengh/min(reciprocal_cell_lengh)

k_point = KP/(reciprocal_cell_lengh_norm[0]*reciprocal_cell_lengh_norm[1]*reciprocal_cell_lengh_norm[2])

k_mesh = np.round(k_point**(1/3)*reciprocal_cell_lengh_norm)

f = open('KPOINTS', "w")
f.write('# KPPRA = {:.0f} #\n'.format(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))
f.write('0\n')
f.write('Gamma\n')
f.write('{:.0f} {:.0f} {:.0f}'.format(k_mesh[0], k_mesh[1], k_mesh[2]))
f.close()
