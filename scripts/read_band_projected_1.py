import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun, BSVasprun, Procar
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.electronic_structure.core import Spin
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from pymatgen.electronic_structure.plotter import BSPlotter, BSPlotterProjected

plt.rcParams["figure.figsize"] = (8,6)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
    
    
def save_scatter(sc, name, kpoints_ticks, kpoints_label):
    plt.colorbar(sc)
    
    for i in kpoints_ticks:
        plt.plot([i, i], [e_min, e_max], linestyle = '-', linewidth = 1, color = 'black')
    plt.plot([min(contrib[:,0]), max(contrib[:,0])], [0, 0], linestyle = '-', linewidth = 1, color = 'black')
        
    plt.xticks(kpoints_ticks)
    plt.gca().set_xticklabels(kpoints_label)
    
    #plt.xlim(min(contrib[:,0]), max(contrib[:,0]))
    plt.xlim(kpoints_ticks[0], kpoints_ticks[-1])
    plt.ylim(e_min, e_max)
    
    plt.tick_params(axis='both',  # Применяем параметры к обеим осям
                   which='major',  # Применяем параметры к вспомогательным делениям
                   direction='in',  # Рисуем деления внутри и снаружи графика
                   # length = 10,    #  Длинна делений
                   # width = 2,     #  Ширина делений
                   # color = 'm',    #  Цвет делений
                   #pad=10,  # Расстояние между черточкой и ее подписью
                   labelsize=14,  # Размер подписи
                   labelcolor='k',  # Цвет подписи
                   bottom=False,  # Рисуем метки снизу
                   top=True,  # сверху
                   left=True,  # слева
                   right=True)  # и справа
    
    plt.savefig('{}.png'.format(name), transparent = False, bbox_inches = 'tight')

    plt.close()
    

e_min = -3
e_max = 3

sites = input('Enter the numbers of atoms for pDOS\n').split()
sites = np.array(sites, dtype = int)

for site in sites:

    try:
        xml = Vasprun('vasprun.xml', parse_projected_eigen=True)
    except:
        print('vasprun.xml not found')
                    
    if xml.converged:
        bands = xml.get_band_structure(kpoints_filename = 'KPOINTS', line_mode = True)

        kpoints_label = ['${}$'.format(bands.kpoints[0].label)]
        for k in range(1, len(bands.kpoints)):
            if bands.kpoints[k].label is not None:
                kpoints_label += [bands.kpoints[k].label]
    
        discontinuous = np.zeros(len(kpoints_label), dtype=bool)
        for i in range(2, len(kpoints_label), 2):
            if kpoints_label[i] != kpoints_label[i-1]:
                discontinuous[i] = True
            
        for spin in [Spin.up, Spin.down]: 
   
            contrib = []
            k_n = 1
            kpoints_ticks = np.zeros(len(kpoints_label))
            kpoints_coords = 0
        
            for k in range(1, len(bands.kpoints)):
                kpoints_coords += np.linalg.norm( np.dot(bands.lattice_rec.matrix/(2*np.pi), bands.kpoints[k].frac_coords - bands.kpoints[k-1].frac_coords ))

                if bands.kpoints[k].label is not None:
                    if discontinuous[k_n] == False:
                        kpoints_ticks[k_n] = kpoints_coords    
                    else:
                        kpoints_coords = kpoints_ticks[k_n-1]
                        kpoints_ticks[k_n] = kpoints_coords
                    k_n += 1 
                
                for b in range(bands.nb_bands):

                    eg = bands.projections[spin][b, k, 6, site] + bands.projections[spin][b, k, 8, site]
                    t2g = bands.projections[spin][b, k, 4, site] + bands.projections[spin][b, k, 5, site] + bands.projections[spin][b, k, 7, site]
        
                    contrib += [[kpoints_coords, bands.bands[spin][b, k] - bands.efermi, t2g, eg]]   # 0:s 1:py 2:pz 3:px 4:dxy 5:dyz 6:dz2 7:dxz 8:dx2_y2 
            contrib = np.array(contrib)
        
        
            kpoints_label_cut = [kpoints_label[0]]
            kpoints_ticks_cut = [kpoints_ticks[0]]
            for k in range(1, len(kpoints_label)):
                if kpoints_label[k] != kpoints_label[k-1]:
                    if discontinuous[k] == False:
                        kpoints_label_cut += ['${}$'.format(kpoints_label[k])]
                    else:
                        kpoints_label_cut += ['${}/{}$'.format(kpoints_label[k-1], kpoints_label[k])]                    
                    kpoints_ticks_cut += [kpoints_ticks[k]]
                                         
    

            sc = plt.scatter(contrib[:,0], contrib[:,1], c = contrib[:,2] , cmap=truncate_colormap(plt.get_cmap('hot'), maxval=0.8), s = 5)
            save_scatter(sc, 't2g_{}_{}'.format(spin, bands.structure[site].species), kpoints_ticks_cut, kpoints_label_cut)

            sc = plt.scatter(contrib[:,0], contrib[:,1], c = contrib[:,3] , cmap=truncate_colormap(plt.get_cmap('hot'), maxval=0.8), s = 5)
            save_scatter(sc, 'eg_{}_{}'.format(spin, bands.structure[site].species), kpoints_ticks_cut, kpoints_label_cut)

    else:
        print('not converged')
