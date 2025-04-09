import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.electronic_structure.core import Spin
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

e_min = -2
e_max = 2

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
mag = ['FIM1', 'AFM_l', 'AFM_s']
lattice = ['regular', 'inverse']

sites = input('Enter the numbers of atoms for pDOS\n').split()
sites = np.array(sites, dtype = int)

#if not os.path.isdir('DOS'):
#    os.makedirs('DOS')

for m in range(len(lattice)):
    for i in range(len(z_atoms)):

        print(lattice[m], z_atoms[i])     
                
        xml = Vasprun('/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Rh2FeZ/PBE/band/{}/step1/{}/vasprun.xml'.format(lattice[m], z_atoms[i]))
            
        if xml.converged:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel(r"Z atom", size=18, labelpad = 0.0)
            ax.set_ylabel(r"$\Delta$ Energy (eV/atom)", size=18, labelpad = -1.0)

            t_energies = xml.tdos.energies - xml.tdos.efermi
            t_dos_up = xml.tdos.get_densities(Spin(1))
            t_dos_down = -xml.tdos.get_densities(Spin(-1))
            t_data = np.column_stack((t_energies, t_dos_up, t_dos_down))

            t_condition = ((t_data[:,0]>e_min) & (t_data[:,0]<e_max))
            ax.fill_between(t_data[t_condition,0], t_data[t_condition,1], t_data[t_condition,2], alpha=0.5, label = 'Total')
                
                    
            ax.plot([e_min, e_max], [0.0, 0.0], '--', linewidth=0.8, color='black')#Ноли
            ax.plot([0.0, 0.0], [min(t_data[t_condition,2]), max(t_data[t_condition,1])], '--', linewidth=1, color='black') #Ноли
                    
            ax.text(0.1, 0.1, '{:.0f} %'.format(xml.complete_dos.spin_polarization*100), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = 18)
                    
            ax.text(0.9, 0.1, z_atoms[i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = 18)


            for s in sites:
                energies = xml.complete_dos.get_site_dos(xml.structures[-1][s]).energies - xml.complete_dos.get_site_dos(xml.structures[-1][s]).efermi
                dos_up = xml.complete_dos.get_site_dos(xml.structures[-1][s]).get_densities(Spin(1))
                dos_down = -xml.complete_dos.get_site_dos(xml.structures[-1][s]).get_densities(Spin(-1))
                        
                data = np.column_stack((energies, dos_up, dos_down))
                        
                condition = ((data[:,0]>e_min) & (data[:,0]<e_max))
                        
                ax.fill_between(data[condition,0], data[condition,1], data[condition,2], alpha=0.5, label = xml.structures[-1][s].species)
                        
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
                        
            fig.savefig('DOS_{}_{}.png'.format(lattice[m], z_atoms[i]), transparent = False, bbox_inches = 'tight')
                        
        else:
            print('{}/{} not converged'.format(lattice[m], z_atoms[i]))        
