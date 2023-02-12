import numpy as np
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
import matplotlib.pyplot as plt
import mpltern
from matplotlib.patches import ArrowStyle, FancyArrowPatch
import matplotlib.colors as colors


def barycentric_to_cartesian(a,b,c):
    x = 0.5 * (2.0*b + c) / (a + b + c)
    y = 0.5*np.sqrt(3)*c / (a + b + c)
    return x, y

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


atoms_order = ['Fe', 'Rh', 'Al']

path = '/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Rh2FeZ/no_SOC/SCAN/hull_energy/Al/small/data.dat'
data = np.loadtxt(path, dtype='str')

studied_structures = {
    'Rh2FeAl': '/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Rh2FeZ/no_SOC/SCAN/ION/isif3/regular/min2/Al/FIM1',
    'Fe2RhAl': '/home/buche/VaspTesting/Danil/X2YZ_Half_metall/Fe2RhZ/no_SOC/SCAN/ION/isif3/XA/FIM1/Al/min1'}

studied_structures_normalize = []
for k in studied_structures.keys():
    comp = Structure.from_file('{}/CONTCAR'.format(studied_structures[k])).composition.get_el_amt_dict()
    data = np.row_stack((data, [k, comp[atoms_order[0]] / sum(comp.values()), comp[atoms_order[1]] / sum(comp.values()),
                                comp[atoms_order[2]] / sum(comp.values()),
                                Outcar('{}/OUTCAR'.format(studied_structures[k])).final_energy / sum(comp.values())]))
    studied_structures_normalize += [Composition(k).formula.replace(' ', '')]

data[:, 0] = [Composition(i.split('_')[0]).formula.replace(' ', '') for i in
              data[:, 0]]  # Удалить отметки Material Project
data = np.array(sorted(data, key=lambda x: (x[4])))

min_energy_index = []
for a in np.unique(data[:, 0]):
    min_energy_index += [np.where(data[:, 0] == a)[0][-1]]


data_min_energy = []
for i in min_energy_index:
    data_min_energy += [data[i]]
data_min_energy = np.array(data_min_energy)
print(data_min_energy)


coords = np.array(data_min_energy[:, 1:4], dtype=float)
labels = np.array(data_min_energy[:, 0], dtype=str)
energies = []
for i in range(len(data_min_energy[:, 4])): #Энергия на формульную единицу
    comp = Composition(data_min_energy[i, 0])
    energies += [float(data_min_energy[i, 4])*sum(comp.values())]

all_prod = dict(zip(labels, energies))

pdes = []
for k, v in all_prod.items():
    pdes += [PDEntry(k, v)]


energy_above_hull = []
pd = PhaseDiagram(pdes)
for p in pdes:
    energy_above_hull += [pd.get_e_above_hull(p)]



fig = plt.figure(figsize=(10.8, 4.8))
#fig.subplots_adjust(left=0.075, right=0.85, wspace=0.3)
ax = fig.add_subplot(projection='ternary')

arrowstyle = ArrowStyle('simple', head_length=10, head_width=5)
kwargs_arrow = {
    'transform': ax.transAxes,  # Used with ``ax.transAxesProjection``
    'arrowstyle': arrowstyle,
    'linewidth': 1,
    'clip_on': False,  # To plot arrows outside triangle
    'zorder': -10,  # Very low value not to hide e.g. tick labels.
}

# Start of arrows in barycentric coordinates.
ta = np.array([ 0.0, -0.1,  1.1])
la = np.array([ 1.1,  0.0, -0.1])
ra = np.array([-0.1,  1.1,  0.0])

# End of arrows in barycentric coordinates.
tb = np.array([ 1.0, -0.1,  0.1])
lb = np.array([ 0.1,  1.0, -0.1])
rb = np.array([-0.1,  0.1,  1.0])

# This transforms the above barycentric coordinates to the original Axes
# coordinates. In combination with ``ax.transAxes``, we can plot arrows fixed
# to the Axes coordinates.
f = ax.transAxesProjection.transform

tarrow = FancyArrowPatch(f(ta), f(tb), ec='C0', fc='C0', **kwargs_arrow)
larrow = FancyArrowPatch(f(la), f(lb), ec='C1', fc='C1', **kwargs_arrow)
rarrow = FancyArrowPatch(f(ra), f(rb), ec='C2', fc='C2', **kwargs_arrow)
ax.add_patch(tarrow)
ax.add_patch(larrow)
ax.add_patch(rarrow)

# To put the axis-labels at the positions consistent with the arrows above, it
# may be better to put the axis-label-text directly as follows rather than
# using e.g.  ax.set_tlabel.
kwargs_label = {
    'transform': ax.transTernaryAxes,
    'backgroundcolor': 'w',
    'ha': 'center',
    'va': 'center',
    'rotation_mode': 'anchor',
    'zorder': -9,  # A bit higher on arrows, but still lower than others.
}

# Put axis-labels on the midpoints of arrows.
tpos = (ta + tb) * 0.5
lpos = (la + lb) * 0.5
rpos = (ra + rb) * 0.5

ax.text(*tpos, 'Fe'  , color='C0', rotation=-60, **kwargs_label, size = 16)
ax.text(*lpos, 'Rh' , color='C1', rotation= 60, **kwargs_label, size = 16)
ax.text(*rpos, 'Al', color='C2', rotation=  0, **kwargs_label, size = 16)

'''
ax.set_tlabel('Fe', fontsize = 16)
ax.set_llabel('Rh', fontsize = 16)
ax.set_rlabel('Al', fontsize = 16)
'''
ax.grid(True, lw = 1)
ax.set_axisbelow(True)

pc = ax.scatter(coords[:,0], coords[:,1], coords[:,2], c=energy_above_hull, cmap = generate_cmap("#ef476f","#ffd166","#06d6a0", points = [0.0, 0.2, 1.0]), s = 130, alpha = 1)

#for i in range(len(labels)):
#    ax.text(coords[i,0], coords[i,1], coords[i,2], labels[i], ha='center', va='top', fontsize = 7)

cax = ax.inset_axes([1.1, 0.1, 0.05, 0.8], transform=ax.transAxes)
colorbar = fig.colorbar(pc, cax=cax)
colorbar.set_label('Energy above hull (meV)', rotation=270, va='baseline', fontsize = 16)

ax.set_ternary_min(-0.08, -0.08, -0.08)
ax.set_ternary_max(1.08, 1.08, 1.08)

plt.savefig('ternary.png', transparent = False, bbox_inches = 'tight')