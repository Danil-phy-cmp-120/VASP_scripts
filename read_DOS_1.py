import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.electronic_structure.core import Spin
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


e_min = -2
e_max = 2
n = 4

plot_total = bool(input('Do plot TDOS (True/False)\n').split())

sites = input('Enter the numbers of atoms for pDOS\n').split()
sites = np.array(sites, dtype = int)

colors = ["#ed0062","#4300ed","#00d37b","#8400e2","#c5f700"]

xml = Vasprun('vasprun.xml')

if xml.converged:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"$E-E_F$", size=18, labelpad = 0.0)
    ax.set_ylabel(r"DOS (st./eV/atom)", size=18, labelpad = -1.0)

    if plot_total == True:
        t_energies = xml.tdos.energies - xml.tdos.efermi
        t_dos_up = xml.tdos.get_densities(Spin(1))/n
        t_dos_down = -xml.tdos.get_densities(Spin(-1))/n
        t_data = np.column_stack((t_energies, t_dos_up, t_dos_down))

        t_condition = ((t_data[:,0]>e_min) & (t_data[:,0]<e_max))
        ax.fill_between(t_data[t_condition,0], t_data[t_condition,1], t_data[t_condition,2], alpha=0.5, label = 'Total', color = colors[0])
        ax.plot(t_data[t_condition,0], t_data[t_condition,1], alpha=1, color = colors[0])
        ax.plot(t_data[t_condition,0], t_data[t_condition,2], alpha=1, color = colors[0])
        
        ax.plot([e_min, e_max], [0.0, 0.0], '--', linewidth=0.8, color='black')#Ноли
        ax.plot([0.0, 0.0], [min(t_data[t_condition,2]), max(t_data[t_condition,1])], '--', linewidth=1, color='black') #Ноли
        
 
    ax.text(0.1, 0.1, '{:.0f} %'.format(xml.complete_dos.spin_polarization*100), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = 18)
        
    n_colors = 1
    for s in sites:
        energies = xml.complete_dos.get_site_dos(xml.structures[-1][s]).energies - xml.complete_dos.get_site_dos(xml.structures[-1][s]).efermi
        dos_up = xml.complete_dos.get_site_dos(xml.structures[-1][s]).get_densities(Spin(1))
        dos_down = -xml.complete_dos.get_site_dos(xml.structures[-1][s]).get_densities(Spin(-1))
            
        data = np.column_stack((energies, dos_up, dos_down))
                        
        condition = ((data[:,0]>e_min) & (data[:,0]<e_max))
                        
        ax.fill_between(data[condition,0], data[condition,1], data[condition,2], alpha=0.5, label = xml.structures[-1][s].specie, color = colors[n_colors])
        ax.plot(data[condition,0], data[condition,1], alpha=1, color = colors[n_colors])
        ax.plot(data[condition,0], data[condition,2], alpha=1, color = colors[n_colors])
        n_colors += 1
        
        if plot_total == False:
            ax.plot([e_min, e_max], [0.0, 0.0], '--', linewidth=0.8, color='black')#Ноли
            ax.plot([0.0, 0.0], [min(data[condition,2]), max(data[condition,1])], '--', linewidth=1, color='black') #Ноли
                        
    #Добавить легенду
    legend = ax.legend(fontsize = 11,
                       ncol = 1,    #  количество столбцов
                       loc='best',
                       #bbox_to_anchor=(0.10, 0.99),
                       facecolor = 'white',    #  цвет области
                       framealpha = 1,
                       #edgecolor = 'None',    #  цвет крайней линии
                       #title = ' ',    #  заголовок
                       #title_fontsize = 20   #  размер шрифта заголовка
                            )
    legend.get_title().set_fontsize('14')
    
    ax.tick_params(axis='both',  # Применяем параметры к обеим осям
                   which='major',  # Применяем параметры к вспомогательным делениям
                   direction='in',  # Рисуем деления внутри и снаружи графика
                   # length = 10,    #  Длинна делений
                   # width = 2,     #  Ширина делений
                   # color = 'm',    #  Цвет делений
                   #pad=10,  # Расстояние между черточкой и ее подписью
                   labelsize=14,  # Размер подписи
                   labelcolor='k',  # Цвет подписи
                   bottom=True,  # Рисуем метки снизу
                   top=True,  # сверху
                   left=True,  # слева
                   right=True,  # и справа                   
                   labelbottom = True,    #  Рисуем подписи снизу
                   labeltop = False,    #  сверху
                   labelleft = True,    #  слева
                   labelright = False,    #  и справа
                   labelrotation = 0)    #  Поворот подписей)  # и справа
                        
    fig.savefig('DOS.png', transparent = False, bbox_inches = 'tight')
                        
else:
        print('not converged')         
