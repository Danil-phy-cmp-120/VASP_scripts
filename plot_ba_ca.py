import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.inputs import Poscar
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt 

if os.path.exists('Relax.dat'):
    data = np.loadtxt('Relax.dat')
else:
    data = []

    for i in sorted([s for s in os.listdir('ba_ca_a0') if s.replace(".", "").isdigit()]):
        for j in sorted([p for p in os.listdir('ba_ca_a0/{}'.format(i)) if p.replace(".", "").isdigit()]):
            print(i, j)
            try:
                outcar = Outcar('ba_ca_a0/{}/{}/OUTCAR'.format(i, j))
                data += [[i, j, outcar.final_energy_wo_entrp]]
            except:
                continue

    data = np.array(data, dtype='f')
    np.savetxt('Relax.dat', data, fmt='%.3f %.3f %.6f')

N = 100
xi = np.linspace(data[:,0].min(), data[:,0].max(), N) 
yi = np.linspace(data[:,1].min(), data[:,1].max(), N) 
zi = griddata((data[:,0], data[:,1]), data[:,2], (xi[None,:], yi[:,None]), method='cubic') 


fig = plt.figure() 
ax = fig.add_subplot(111)

ax.set_xlabel(r"b/a", size=16)
ax.set_ylabel(r"c/a", size = 16) 

lev = np.linspace(np.floor(min(data[:,2])), np.around(max(data[:,2])), 30)
cp = ax.contourf(xi, yi, zi, cmap = 'hot', levels = lev)
cbar = plt.colorbar(cp, extend='both', shrink=0.9)#, ticks=ticks)

#  Настраиваем вид вспомогательных тиков:
ax.tick_params(axis = 'both',    #  Применяем параметры к обеим осям
               which = 'major',    #  Применяем параметры к вспомогательным делениям
               direction = 'in',    #  Рисуем деления внутри и снаружи графика
               #length = 10,    #  Длинна делений
               #width = 2,     #  Ширина делений
               #color = 'm',    #  Цвет делений
               #pad = 10,    #  Расстояние между черточкой и ее подписью
               labelsize = 12,    #  Размер подписи
               #labelcolor = 'r',    #  Цвет подписи
               bottom = True,    #  Рисуем метки снизу
               top = True,    #   сверху
               left = True,    #  слева
               right = True)    #  и справа
               
fig.savefig('ba_ca.png', dpi=300)
