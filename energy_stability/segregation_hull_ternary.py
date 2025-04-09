# -*- coding: utf-8 -*-
# version: 2.0.0 #

from chempy import balance_stoichiometry
import numpy as np
import itertools
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
import matplotlib.pyplot as plt
import mpltern
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16)/255 for i in range(0, lv, lv // 3))

def generate_cmap(color_finish, color_middle, color_start, points = [0.0, 0.5, 1.0]):
    RGB1 = hex_to_rgb(color_start)
    RGB2 = hex_to_rgb(color_middle)
    RGB3 = hex_to_rgb(color_finish)

    cdict = {'red' :   [( points[0] ,  RGB1[0],  RGB1[0]),
                        ( points[1] ,  RGB2[0],  RGB2[0]),
                        ( points[2],  RGB3[0] , RGB3[0] )],

             'green' : [( points[0] ,  RGB1[1] , RGB1[1]),
                        ( points[1] ,  RGB2[1],  RGB2[1]),
                        ( points[2] ,  RGB3[1] , RGB3[1])],

             'blue' :  [( points[0] ,  RGB1[2] , RGB1[2]),
                        ( points[1],  RGB2[2],  RGB2[2]),
                        ( points[2] ,  RGB3[2] , RGB3[2])],}
    return(colors.LinearSegmentedColormap('my_colormap', cdict))


num_prod = [2, 3]
atoms_order = ['Fe', 'Rh', 'Ga']

path = '/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Rh2FeZ/no_SOC/SCAN/hull_energy/Ga/small/data.dat'
data = np.loadtxt(path, dtype='str')

studied_structures = {'Rh2FeGa': '/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Rh2FeZ/no_SOC/SCAN/ION/isif3/regular/min2/Ga/FIM1', 'Fe2RhGa': '/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Fe2RhZ/no_SOC/SCAN/ION/isif3/XA/FIM1/Ga/min1'}
for k in studied_structures.keys():
    comp = Structure.from_file('{}/CONTCAR'.format(studied_structures[k])).composition.get_el_amt_dict()
    data = np.row_stack((data, [k, comp[atoms_order[0]]/sum(comp.values()), comp[atoms_order[1]]/sum(comp.values()), comp[atoms_order[2]]/sum(comp.values()), Outcar('{}/OUTCAR'.format(studied_structures[k])).final_energy/sum(comp.values())] ))

data[:,0] = [Composition(i.split('_')[0]).formula.replace(' ', '') for i in data[:,0]] # Удалить отметки Material Project
data = np.array(sorted(data, key=lambda x: (x[4])))

min_energy_index = []
for a in np.unique(data[:,0]):
    min_energy_index += [np.where(data[:,0]==a)[0][-1]]

data_min_energy = []
for i in min_energy_index:
    data_min_energy += [data[i]]
data_min_energy = np.array(data_min_energy)

coords = np.array(data_min_energy[:,1:4], dtype = float)
all_prod = dict(zip(data_min_energy[:,0], [float(i) for i in data_min_energy[:,4]]))

for r in all_prod.keys():
    atoms = Composition(r).get_el_amt_dict()
    for a in atoms_order:
        all_prod[r] -= (atoms[a]/sum(atoms.values())) * all_prod[a+'1']


reactions_all = []
num_reactions = np.zeros((len(all_prod.keys()), 2))
for ri, r in enumerate(all_prod.keys()):
    for n in num_prod:
        for p in itertools.combinations(all_prod.keys(), n):
            try:
                reac, prod = balance_stoichiometry([r], p)
                if (sum(1 for number in prod.values() if number < 0)) == 0:
                    e_formation = reac[r]*all_prod[r]
                    reactions = '{} {} ->'.format(reac[r], r)
                    for j in prod.keys():
                        e_formation -= prod[j] * all_prod[j]
                        reactions += ' {} {} '.format(prod[j], j)
                    reactions += '{}'.format(e_formation)
                    
                    if e_formation > 0:
                        reactions += ' {}\n'.format(False)
                        num_reactions[ri, 0] += 1 
                    else:
                        reactions += '{}\n'.format(True)
                        num_reactions[ri, 1] += 1 
                    
                    reactions_all += [reactions]
            except:
                continue

reaction_percent = np.zeros(len(all_prod.keys()))
for i in range(len(num_reactions)):
    if num_reactions[i,0] == 0 and num_reactions[i,1] == 0:
        reaction_percent[i] = 100
    else:
        reaction_percent[i] = 100*num_reactions[i,1]/(num_reactions[i,0] + num_reactions[i,1])
print(reaction_percent)


fig = plt.figure(figsize=(10.8, 4.8))
#fig.subplots_adjust(left=0.075, right=0.85, wspace=0.3)
ax = fig.add_subplot(projection='ternary')

ax.set_tlabel('Fe', fontsize = 16)
ax.set_llabel('Rh', fontsize = 16)
ax.set_rlabel('Ga', fontsize = 16)

ax.grid()

pc = ax.scatter(coords[:,0], coords[:,1], coords[:,2], c=reaction_percent, cmap = generate_cmap("#06d6a0","#ffd166","#ef476f"), s = 100)
#for i in range(len(data_min_energy[:,0])):
#    ax.text(coords[i,0], coords[i,1], coords[i,2], data_min_energy[i,0], ha='center', va='center')

cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9], transform=ax.transAxes)
colorbar = fig.colorbar(pc, cax=cax)
colorbar.set_label('Stability (%)', rotation=270, va='baseline', fontsize = 1)

ax.set_ternary_min(-0.08, -0.08, -0.08)
ax.set_ternary_max(1.08, 1.08, 1.08)


plt.savefig('ternary.png')

