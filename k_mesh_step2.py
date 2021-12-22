#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import numpy as np
import ase.io.vasp

def get_params(path, mag_type):
    f = open(path + '/OUTCAR',"r")
    outcar = f.readlines()
    f.close()
    
    for i in np.arange(0, len(outcar)):
        inp = outcar[i].split()
        if len(inp) > 3 and inp[0] == 'number' and inp[1] == 'of' and inp[2] == 'dos' :
            nions = int(inp[11])

    if mag_type == 'std':    
        flag = False
        magmom = []
        for line in outcar:
            inp = line.split()

            if len(inp) > 3 and flag == True and inp[0].isdigit() == True:
                magmom += [float(inp[len(inp)-1])]
            if len(inp) > 1 and inp[0] == 'magnetization':
                flag = True
    else:       
        for i in np.arange(len(outcar)-2000, len(outcar)):
            inp = outcar[i].split()
            if len(inp) > 1 and inp[0] == 'magnetization' and inp[1] == '(x)':
                n = i
            
        magmom_x = []; magmom_y = []; magmom_z = []
        for line in outcar[n+3:n+nions+4]:
            inp = line.split()
            if len(inp) >1:
                magmom_x += [float(inp[4])]   
        for line in outcar[n+nions+13:n+2*nions+13]:
            inp = line.split()
            if len(inp) >1:
                magmom_y += [float(inp[4])]  
        for line in outcar[n+2*nions+22:n+3*nions+22]:
            inp = line.split()
            if len(inp) >1:
                magmom_z += [float(inp[4])]

    if mag_type == 'std': 
        return(magmom)
    else:
        return np.column_stack((magmom_x, magmom_y, magmom_z))



#path0 = input('Enter the path to the step1 calculation:\n')
#directory = os.listdir(path0)

mag_type = input('std or ncl?:\n')
KPPRA = input('Enter the approximate number of k-points per reciprocal atom (default value = 2000):\n')
KPPRA = int(KPPRA)

path0 = 'step1'

if os.path.exists('step2') == False:
    os.mkdir('step2')

shutil.copyfile('{}/POTCAR'.format(path0), 'step2/POTCAR') # Запись POTCAR
shutil.copyfile('{}/POSCAR'.format(path0), 'step2/POSCAR') # Запись POSCAR
os.link('{}/CHGCAR'.format(path0), 'step2/CHGCAR') # Запись CHGCAR
os.link('{}/WAVECAR'.format(path0), 'step2/WAVECAR') # Запись WAWECAR


cell = ase.io.vasp.read_vasp('{}/POSCAR'.format(path0))
KP = KPPRA/len(cell)

reciprocal_cell_lengh = np.array([np.linalg.norm(np.linalg.inv(cell.cell)[0]), np.linalg.norm(np.linalg.inv(cell.cell)[1]), np.linalg.norm(np.linalg.inv(cell.cell)[2])])

reciprocal_cell_lengh_norm = reciprocal_cell_lengh/min(reciprocal_cell_lengh)

k_point = KP/(reciprocal_cell_lengh_norm[0]*reciprocal_cell_lengh_norm[1]*reciprocal_cell_lengh_norm[2])

k_mesh = np.round(k_point**(1/3)*reciprocal_cell_lengh_norm)

f = open('step2/KPOINTS', "w")
f.write('# KPPRA = {:.0f} #\n'.format(k_mesh[0]*k_mesh[1]*k_mesh[2]*len(cell)))
f.write('0\n')
f.write('Gamma\n')
f.write('{:.0f} {:.0f} {:.0f}'.format(k_mesh[0], k_mesh[1], k_mesh[2]))
f.close()


######### Изменение некоторых тегов в INCAR #########

magmom = get_params(path0, mag_type)

f = open('{}/INCAR'.format(path0),"r")
incar = f.readlines()
f.close()

for i in range(len(incar)):
    inp = incar[i].split()
    if len(inp) > 2 and inp[0] == 'LCHARG':
        incar[i] = 'LCHARG = .FALSE.\n'
    if len(inp) > 2 and inp[0] == 'LWAVE':
        incar[i] = 'LWAVE = .FALSE.\n'
    if len(inp) > 2 and inp[0] == 'NSW':
        incar[i] = 'NSW = 0\n'
        
    if len(inp) > 2 and inp[0] == 'MAGMOM':
        if mag_type == 'std':
            incar[i] = 'MAGMOM ='
            for j in range(len(magmom)):
                incar[i] = incar[i] + ' {}'.format(magmom[j])
        else:
            incar[i] = 'MAGMOM ='
            for j in range(len(magmom)):
                incar[i] = incar[i] + ' {} {} {}'.format(magmom[j,0], magmom[j,1], magmom[j,2])
        incar[i] = incar[i] + '\n'

f = open('step2/INCAR', "w") # Запись INCAR
for i in range(len(incar)):
    f.write(incar[i])
f.close

