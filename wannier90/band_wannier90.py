# -*- coding: utf-8 -*-

#Подключение библиотеки для рисования
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import pylab

e_min = -2.5
e_max = 2.5

#####  Read OUTCAR  #####
outcar = open('OUTCAR',"r")
lines = outcar.readlines()
outcar.close()

for line in lines:
    inp = line.split()

    if len(inp) > 3 and inp[0] == 'E-fermi' and inp[1] == ':':
        Efermi = float(inp[2])

band_up = np.loadtxt('wannier90.1_band.dat')
band_down = np.loadtxt('wannier90.2_band.dat')
labels = np.loadtxt('wannier90.1_band.labelinfo.dat', usecols=[0], dtype=str)
positions = np.loadtxt('wannier90.1_band.labelinfo.dat', usecols=[2])
band_up[:,0] = band_up[:,0]/max(positions)
band_down[:,0] = band_down[:,0]/max(positions)
positions = positions/max(positions)

labels = list(labels)

for i in range(len(labels)):
    labels[i] = labels[i].replace('GAMMA', 'Gamma')
    labels[i] = labels[i].replace('SIGMA', 'Sigma')

labels_new = []
i = 0
while i < len(labels)-1:

    if positions[i] != positions[i+1]:
        labels_new += [labels[i]]
    else:
        labels_new += [labels[i] +'|' + labels[i+1]] 
        i += 1
    i += 1
labels_new += [labels[len(labels)-1]]

positions = np.unique(positions)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Wave Vector", size=18, labelpad = 0.0)
ax.set_ylabel(r"$E - E_{F}$, eV", size=18, labelpad = -1.0)

ax.plot([positions[0], positions[positions.size-1]], [0.0, 0.0], color ='dimgray', linestyle = '-', linewidth = 1)
for p in positions:
    ax.plot([p, p], [e_min, e_max], color ='dimgray', linestyle = '-', linewidth = 1)

ax.plot(band_up[:,0], band_up[:,1]-Efermi, color ='#fb8500', marker = 'o', markersize = 2.5, linestyle = 'None', alpha = 1.0, markeredgewidth = 0)

ax.set_xticks(positions)
ax.set_xticklabels(labels_new)

#ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

'''
#Добавить легенду
legend = ax.legend(fontsize = 11,
          ncol = 2,    #  количество столбцов
          loc='upper left',
          bbox_to_anchor=(0.10, 0.99),
          facecolor = 'white',    #  цвет области
          framealpha = 1,
          #edgecolor = 'None',    #  цвет крайней линии
          title = ' ',    #  заголовок
          #title_fontsize = 20   #  размер шрифта заголовка
          )
legend.get_title().set_fontsize('14')
'''

#Название графика
#ax.set_title(u'')
#Добавление сетки

ax.grid(which='major',
        color = 'lightgray',
        linestyle = ':',
        linewidth = 1)
'''
ax.minorticks_on()

ax.grid(which='minor',
        color = 'gray',
        linestyle = ':')
'''
ax.tick_params(axis = 'both',    #  Применяем параметры к обеим осям
               which = 'major',    #  Применяем параметры к вспомогательным делениям
               direction = 'in',    #  Рисуем деления внутри и снаружи графика
               #length = 10,    #  Длинна делений
               #width = 2,     #  Ширина делений
               #color = 'm',    #  Цвет делений
               #pad = 10,    #  Расстояние между черточкой и ее подписью
               #labelsize = 15,    #  Размер подписи
               #labelcolor = 'r',    #  Цвет подписи
               bottom = True,    #  Рисуем метки снизу
               top = True,    #   сверху
               left = True,    #  слева
               right = True)    #  и справа

plt.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlim(0.0, 1.0)
ax.set_ylim(e_min, e_max)

plt.tight_layout()

#Сохранения изображения
fig.savefig('band_up.png', dpi=300, transparent=False)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Wave Vector", size=18, labelpad = 0.0)
ax.set_ylabel(r"$E - E_{F}$, eV", size=18, labelpad = -1.0)

ax.plot([positions[0], positions[positions.size-1]], [0.0, 0.0], color ='dimgray', linestyle = '-', linewidth = 1)
for p in positions:
    ax.plot([p, p], [e_min, e_max], color ='dimgray', linestyle = '-', linewidth = 1)

ax.plot(band_down[:,0], band_down[:,1]-Efermi, color ='#219ebc', marker = 'o', markersize = 2.5, linestyle = 'None', alpha = 1.0, markeredgewidth = 0)

ax.set_xticks(positions)
ax.set_xticklabels(labels_new)

#ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

'''
#Добавить легенду
legend = ax.legend(fontsize = 11,
          ncol = 2,    #  количество столбцов
          loc='upper left',
          bbox_to_anchor=(0.10, 0.99),
          facecolor = 'white',    #  цвет области
          framealpha = 1,
          #edgecolor = 'None',    #  цвет крайней линии
          title = ' ',    #  заголовок
          #title_fontsize = 20   #  размер шрифта заголовка
          )
legend.get_title().set_fontsize('14')
'''

#Название графика
#ax.set_title(u'')
#Добавление сетки

ax.grid(which='major',
        color = 'lightgray',
        linestyle = ':',
        linewidth = 1)
'''
ax.minorticks_on()

ax.grid(which='minor',
        color = 'gray',
        linestyle = ':')
'''
ax.tick_params(axis = 'both',    #  Применяем параметры к обеим осям
               which = 'major',    #  Применяем параметры к вспомогательным делениям
               direction = 'in',    #  Рисуем деления внутри и снаружи графика
               #length = 10,    #  Длинна делений
               #width = 2,     #  Ширина делений
               #color = 'm',    #  Цвет делений
               #pad = 10,    #  Расстояние между черточкой и ее подписью
               #labelsize = 15,    #  Размер подписи
               #labelcolor = 'r',    #  Цвет подписи
               bottom = True,    #  Рисуем метки снизу
               top = True,    #   сверху
               left = True,    #  слева
               right = True)    #  и справа

plt.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlim(0.0, 1.0)
ax.set_ylim(e_min, e_max)

plt.tight_layout()

#Сохранения изображения
fig.savefig('band_down.png', dpi=300, transparent=False)
