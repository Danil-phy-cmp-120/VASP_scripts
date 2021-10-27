# -*- coding: utf-8 -*-

#Подключение библиотеки для рисования
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import pylab

fig = plt.figure()#figsize=(15.5/2, 9.5/2), dpi= 300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111)
ax.set_xlabel(r"$c/a$", size=20)
ax.set_ylabel(r"$\Delta E$ (eV/atom)", size=20)

a = np.loadtxt('Relax.dat')
a[:,1] = (a[:,1] + 0.12327923E+03)/16
ax.plot(a[:,0],a[:,1],  color='blue', linestyle='-', marker = '^', markersize = 10)

#ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

'''
#Добавить легенду
legend = ax.legend(fontsize = 10,
          ncol = 1,    #  количество столбцов
          loc='lower left',
          #bbox_to_anchor=(0.8, 0.75),
          facecolor = 'None',    #  цвет области
          #edgecolor = 'None',    #  цвет крайней линии
          #title = 'External pressure:',    #  заголовок
          #title_fontsize = 20   #  размер шрифта заголовка
          )
legend.get_title().set_fontsize('11')
'''
'''
legend1 = ax.legend([p1, p2],['D03', 'L12'],loc='upper left',
                     ncol = 1,    #  количество столбцов
                     #bbox_to_anchor=(0.0, 0.46),
                     edgecolor = 'None',    #  цвет крайней линии
                     fontsize = 11,
                     facecolor = 'None')    #  цвет области
                     #title = 'External pressure:'

legend1.get_title().set_fontsize('11')

legend2 = ax.legend([T1, T2, T3],[r'$E$(Fe$_{3}$Al$_{x}$Cr$_{1-x}$) - 3$E$(Fe)', 'SCAN'], loc='upper left',
                     ncol = 1,    #  количество столбцов
                     bbox_to_anchor=(0.15, 0.75),
                     edgecolor = 'None',    #  цвет крайней линии
                     fontsize = 11,
                     facecolor = 'None')   #  цвет области
                     #title = 'External magnetic:'
ax.add_artist(legend1)
legend2.get_title().set_fontsize('11')
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

#  Настраиваем вид вспомогательных тиков:
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


#Наиболее подходящеий интервал между графиками
plt.tight_layout()
#Настройка диапазонов по осям
#ax.set_xlim(-0.05,1.05)
#ax.set_ylim(-0.002,0.027)
#Отобразить окно
#plt.show()

#Сохранения изображения
fig.savefig('ca.png', dpi=300)


