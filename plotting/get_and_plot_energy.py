# -*- coding: utf-8 -*-

#Подключение библиотеки для рисования
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import pylab
import os


def get_energy(p):
    f = open(p + '/OUTCAR', "r")
    outcar = f.readlines()
    f.close()

    for line in outcar:
        inp = line.split()
        if len(inp) > 4 and inp[4] == 'energy(sigma->0)':
            energy = float(inp[6])

    return energy


def test_convergence(p):
    f = open(p + '/OSZICAR', "r")
    oszicar = f.readlines()
    f.close()

    test = oszicar[len(oszicar)-2].split()[1]
    if int(test) >= 101:
        return False
    else:
        return True

paths = sorted(os.listdir(os.getcwd()))

fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111)
ax.set_xlabel(r"Lattice parameter ($\mathrm{\AA}$)", size=16)
ax.set_ylabel(r"$\Delta E$ (eV/atom)", size=16)

energy = []
for p in paths:
    if os.path.isdir(p):
        if test_convergence(p) == True:
            energy += [[float(p), get_energy(p)]]
        else:
            print('{} not converged'.format(p))

energy = np.array(energy)

ax.plot(energy[:,0], energy[:,1]/16, linestyle='-', marker = 'o')
np.savetxt('Relax.dat', energy, fmt = '%.3f %.8f')

#ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

'''
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
'''
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
fig.savefig('Relax.png', dpi=300)
