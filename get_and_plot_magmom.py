# -*- coding: utf-8 -*-

#Подключение библиотеки для рисования
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import pylab
import os

def get_magmom(p, num):
    f = open(p + '/OUTCAR', "r")
    outcar = f.readlines()
    f.close()

    for i in np.arange(len(outcar)-500, len(outcar)):
        inp = outcar[i].split()
        if len(inp) > 1 and inp[0] == 'magnetization' and inp[1] == '(x)':
            n = i

    for line in outcar[n:len(outcar)]:
        inp = line.split()
        if len(inp) > 1 and inp[0] == str(num):
            magmom = float(inp[4])

    return magmom

def get_atoms(p):
    f = open(p, "r")
    poscar = f.readlines()
    f.close()

    atoms = [poscar[5].split()]
    atoms += [poscar[6].split()]

    return atoms

def test_convergence(p):
    f = open(p + '/OSZICAR', "r")
    oszicar = f.readlines()
    f.close()

    test = oszicar[len(oszicar)-2].split()[1]
    if int(test) >= 101:
        return False
    else:
        return True


fu = 4
paths = sorted(os.listdir(os.getcwd()))
atoms = get_atoms(paths[0] + '/POSCAR')

fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111)
ax.set_xlabel(r"Lattice parameter ($\mathrm{\AA}$)", size=16)
ax.set_ylabel(r"$\mu$ ($\mu_B$/f.u.)", size=16)

colors = ['r', 'g', 'b', 'aqua', 'magenta']

i = 0
j = 0
for n in range(1, sum(np.array(atoms[1], dtype = np.int))+1):
    magmom = []
    for p in paths:
        if os.path.isdir(p):
            if test_convergence(p) == True:
                magmom += [[float(p), get_magmom(p, n)]]
            

    magmom = np.array(magmom)
            
    ax.plot(magmom[:,0], magmom[:,1]/fu, linestyle='-', label = atoms[0][i], color = colors[i])

    j = j + 1
    if j > int(atoms[1][i]):
        i = i + 1
        j = 1 

#ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))


#Добавить легенду
legend = ax.legend(fontsize = 10,
          ncol = 1,    #  количество столбцов
          loc='best',
          #bbox_to_anchor=(0.12, 1.0),
          facecolor = 'white',    #  цвет области
          framealpha = 1,
          #edgecolor = 'None',    #  цвет крайней линии
          #title = 'External pressure:',    #  заголовок
          #title_fontsize = 20   #  размер шрифта заголовка
          )
legend.get_title().set_fontsize('10')

ax.grid(which='major',
        color = 'lightgray',
        linestyle = ':',
        linewidth = 1)

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


#Наиболее подходящеий интервал между графиками
plt.tight_layout()
#Настройка диапазонов по осям
#ax.set_xlim(5.7,6.6)
#ax.set_ylim(-0.002,0.027)
#Отобразить окно
#plt.show()

#Сохранения изображения
fig.savefig('magmom.png', dpi=300)
