#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import numpy as np

def get_params(path):
    f = open(path + '/OUTCAR',"r")
    outcar = f.readlines()
    f.close()

    flag = False
    magmom = []
    for line in outcar[len(outcar)-200:len(outcar)]:
        inp = line.split()

        if len(inp) > 3 and flag == True and inp[0].isdigit() == True:
            magmom += [float(inp[len(inp)-1])]
        if len(inp) > 1 and inp[0] == 'magnetization':
            flag = True

    for line in outcar:
        inp = line.split()

        if len(inp) > 13 and inp[13] == 'NBANDS=':
            nbands = int(inp[14])

    return(magmom, nbands)



path0 = input('Enter the path to the SCF calculation:\n')
directory = os.listdir(path0)

quantization_axis = []
print("Enter the spin quantization axis or 'end'")
while True:
    ans = input()
    if ans == 'end':
        break
    else:
        quantization_axis += [ans]
 
for q in quantization_axis:
    if os.path.exists('NSCF_{}_test'.format(q.replace(' ', ''))) == False:
        os.mkdir('NSCF_{}_test'.format(q.replace(' ', '')))

    for d in directory:
        if os.path.isdir('{}/{}'.format(path0, d)) == True:

            if os.path.exists('NSCF_{}_test/{}'.format(q.replace(' ', ''), d)) == False:
                os.mkdir('NSCF_{}_test/{}'.format(q.replace(' ', ''), d))

            shutil.copyfile('{}/{}/POTCAR'.format(path0, d), 'NSCF_{}_test/{}/POTCAR'.format(q.replace(' ', ''), d)) # Запись POTCAR
            shutil.copyfile('{}/{}/POSCAR'.format(path0, d), 'NSCF_{}_test/{}/POSCAR'.format(q.replace(' ', ''), d)) # Запись POSCAR
            shutil.copyfile('{}/{}/KPOINTS'.format(path0, d), 'NSCF_{}_test/{}/KPOINTS'.format(q.replace(' ', ''), d)) # Запись KPOINTS
            os.symlink('{}/{}/CHGCAR'.format(path0, d), 'NSCF_{}_test/{}/CHGCAR'.format(q.replace(' ', ''), d)) # Запись CHGCAR


            ######### Изменение некоторых тегов в INCAR #########

            magmom, nbands = get_params('{}/{}'.format(path0, d))

            f = open('{}/{}/INCAR'.format(path0, d),"r")
            incar = f.readlines()
            f.close()

            for i in range(len(incar)):
                inp = incar[i].split()
                if len(inp) > 2 and inp[0] == 'LCHARG':
                    incar[i] = 'LCHARG = .FALSE.\n'
                if len(inp) > 2 and inp[0] == 'MAGMOM':
                    incar[i] = 'MAGMOM ='
                    for j in range(len(magmom)):
                        incar[i] = incar[i] + ' 0 0 {}'.format(magmom[j])
                    incar[i] = incar[i] + '\n'

            incar = incar + ['LSORBIT = .TRUE.\n']
            incar = incar + ['GGA_COMPAT =.FALSE.\n']
            incar = incar + ['SAXIS = {}\n'.format(q)]
            incar = incar + ['ICHARG = 11\n']
            incar = incar + ['NBANDS = {}\n'.format(int(nbands)*2)]
            incar = incar + ['ISYM = -1']

            print(incar)

            f = open('NSCF_{}_test/{}/INCAR'.format(q.replace(' ', ''), d), "w") # Запись INCAR
            for i in range(len(incar)):
                f.write(incar[i])
            f.close

