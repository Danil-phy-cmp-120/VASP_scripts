# -*- coding: utf-8 -*-
import numpy as np
import ase.io.vasp
import ase.spacegroup
from ase.visualize.plot import plot_atoms
import shutil
import os

for KPPRA in [100, 200, 300, 500, 1000, 2000, 3000, 4000, 5000]:

    cell = ase.io.vasp.read_vasp("POSCAR")
    KP = KPPRA/len(cell)

    reciprocal_cell_lengh = np.array([np.linalg.norm(np.linalg.inv(cell.cell)[0]), np.linalg.norm(np.linalg.inv(cell.cell)[1]), np.linalg.norm(np.linalg.inv(cell.cell)[2])])

    reciprocal_cell_lengh_norm = reciprocal_cell_lengh/min(reciprocal_cell_lengh)

    k_point = KP/(reciprocal_cell_lengh_norm[0]*reciprocal_cell_lengh_norm[1]*reciprocal_cell_lengh_norm[2])

    k_mesh = np.round(k_point**(1/3)*reciprocal_cell_lengh_norm)
    
    if os.path.exists('{}'.format(int(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))) == False:
        os.mkdir('{}'.format(int(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell))))

    f = open('{}/KPOINTS'.format(int(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell))), "w")
    f.write('# KPPRA = {:.0f} #\n'.format(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))
    f.write('0\n')
    f.write('Gamma\n')
    f.write('{:.0f} {:.0f} {:.0f}'.format(k_mesh[0], k_mesh[1], k_mesh[2]))
    f.close()
    
    shutil.copyfile('POTCAR', '{}/POTCAR'.format(int(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))) # Запись POTCAR
    shutil.copyfile('INCAR', '{}/INCAR'.format(int(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))) # Запись INCAR
    shutil.copyfile('POSCAR', '{}/POSCAR'.format(int(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))) # Запись POSCAR
