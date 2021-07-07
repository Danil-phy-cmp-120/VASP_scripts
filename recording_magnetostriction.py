#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import subprocess
import shutil

def toFixed(f, n=0):
    a, b = str(f).split('.')
    return '{}.{}{}'.format(a, b[:n], '0'*(n-len(b)))

Name = 'FeGa'
EPSILON_Z = np.array([-0.015, -0.010, -0.005, 0.0, 0.005, 0.010, 0.015])

# Считываем вектора трансляций и пааметр решетки #
f = open('POSCAR',"r")
lines = f.readlines()
f.close()

primitive_vectors = np.array([[float(lines[2].split()[0]), float(lines[2].split()[1]), float(lines[2].split()[2])],
                             [float(lines[3].split()[0]), float(lines[3].split()[1]), float(lines[3].split()[2])],
                             [float(lines[4].split()[0]), float(lines[4].split()[1]), float(lines[4].split()[2])]])

a0 = float(lines[1].split()[0])
V0 = primitive_vectors[0,0]*a0 * primitive_vectors[1,1]*a0 * primitive_vectors[2,2]*a0

os.mkdir('Distortion')
for ez in EPSILON_Z:
    os.mkdir('Distortion/' + str(ez))
    T = np.array([[1-0.5*ez, 0, 0],
                 [0,1 -0.5*ez, 0],
                 [0, 0, 1 + ez]])

    newcell = np.dot(primitive_vectors, T)
    V = newcell[0,0]*a0 * newcell[1,1]*a0 * newcell[2,2]*a0
    #print V0, V

    lines[2] = toFixed(newcell[0,0], 5) + '  ' + toFixed(newcell[0,1], 5) + '  ' + toFixed(newcell[0,2], 5) + '\n'
    lines[3] = toFixed(newcell[1,0], 5) + '  ' + toFixed(newcell[1,1], 5) + '  ' + toFixed(newcell[1,2], 5) + '\n'
    lines[4] = toFixed(newcell[2,0], 5) + '  ' + toFixed(newcell[2,1], 5) + '  ' + toFixed(newcell[2,2], 5) + '\n'

    file = open('Distortion/' + str(ez) + '/POSCAR', "wb")
    for line in lines:
        file.write(line)
    file.close()             


    shutil.copyfile(os.getcwd()+'/POTCAR', 'Distortion/' + str(ez) + '/POTCAR')
    shutil.copyfile(os.getcwd()+'/KPOINTS', 'Distortion/' + str(ez) + '/KPOINTS')
    shutil.copyfile(os.getcwd()+'/INCAR', 'Distortion/' + str(ez) + '/INCAR')
    shutil.copyfile(os.getcwd()+'/vasp_qsub_1', 'Distortion/' + str(ez) + '/vasp_qsub_1')

    

       


