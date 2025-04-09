# -*- coding: utf-8 -*-
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen.io.vasp.inputs as inputs
from pymatgen.core.structure import Structure
import os
from scipy.constants import pi

def naruto_start_picture():
        print('⠀⠀⠀⠀⠀⠀⠀   ⠀⠀     ⠀⠀⠘⢦⠐⠒⠒⠒⠀⢠⣸⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠉⠒⠒⠒⠢⢤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣘⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡤⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⢿⡉⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡔⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠲⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠳⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⠀⠈⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠖⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⣇⠒⠒⠒⠶⠾⠶⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⢁⣀⡠⡴⠀⡾⠀⠀⠀⣠⠀⠀⠀⠀⣴⣶⣶⣶⣶⣶⣿⣿⣿⣿⣿⣿⣷⣶⣶⣶⡆⠀⠀⠀⢠⠀⠀⠀⣄⠀⠘⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠞⠛⠉⠁⠀⣸⢁⣴⠇⠀⠀⣼⡏⠀⠀⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠁⠀⠀⢀⣼⠀⠀⢠⣿⠛⠢⢬⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣰⡷⢋⡜⠀⢠⣾⣿⡇⠀⠀⡼⠯⠿⠿⠿⠿⢿⠿⠿⠿⠿⠭⠍⠙⠙⠋⠙⢻⠀⠀⢀⣾⣿⠀⠀⣼⣿⡄⠀⠀⠉⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⣤⣤⣤⣤⣤⣶⣤⣶⣶⣶⣶⣾⣷⣶⣶⠇⣠⣿⣿⣿⡇⠀⣸⠁⠀⠀⠀⠀⠀⠀⠀⠀⣠⣴⣦⣴⡖⠀⢀⠀⡇⠀⡴⠛⣿⣿⠀⢠⣿⣿⣁⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⣿⣿⣿⣿⡿⠟⠋⠁⠀⠀⠉⠉⠛⠿⣿⠟⣿⣿⣿⣿⡇⢠⡟⠃⠀⠀⢆⠀⠀⠀⣠⡾⢡⡶⢶⣌⠀⠀⠤⢸⢁⠞⠁⢛⣿⣿⢀⣾⣿⣿⢉⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⡿⠟⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⢿⠷⣿⠙⢿⣿⣷⣿⡇⣄⠀⠀⠀⠀⠀⣴⣏⠻⣜⣛⣻⡟⠀⠀⠀⡿⠋⠀⢀⣸⣿⣿⣿⢿⠋⢸⠟⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣾⣿⣿⠀⣹⡀⢸⠈⠛⠿⣿⣶⣤⣄⣀⣀⠀⠀⠈⠉⠉⠉⠁⠉⠀⣀⣀⣀⣤⣤⣾⣿⠿⠛⠁⢸⠀⣾⡀⣼⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣾⣿⠟⠁⢿⣀⡏⡇⡇⠀⠀⠀⢀⣩⣿⣟⠿⣿⣿⣿⣷⣶⣶⣶⣾⣿⣿⣿⣿⣿⣟⣿⣿⣅⡀⠀⡀⢸⣄⣷⣷⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⢀⣴⡿⠛⠁⠀⠀⠈⣯⢳⣿⡇⠀⠀⠘⡏⢽⣿⣿⢹⣮⡷⣄⠈⣹⠉⠉⠹⣍⣀⣴⣾⠿⣯⣿⢹⡍⢻⡇⠁⠀⣿⡝⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⢀⣠⣶⠟⠉⠀⠀⠀⠀⠀⠀⠹⡄⠛⡇⠀⠀⠀⠳⣄⠻⠿⠞⠁⢳⠈⠉⠁⠀⠀⠀⠨⠟⠁⣞⠀⠳⠿⣞⣡⡞⠀⠀⠀⡿⢀⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⣤⣶⣶⡿⠟⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡿⣄⣽⠀⠀⠀⠀⠈⠛⠛⠛⠛⠙⠠⠴⠆⠀⠀⠀⠀⠀⠀⠈⠛⠛⠛⠋⠁⠀⠀⠀⢠⣷⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠿⠛⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣼⣶⠃⠘⡄⠀⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⡜⠁⠀⠀⡼⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢉⣟⣰⠄⢳⡤⠤⠚⠓⠊⠙⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠑⢲⣇⣿⠀⠀⠀⠀⠀⠀⠀⣀⣀⣤⣄⡀⠀⠀⠀⠀\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠼⠋⢸⠀⣨⣧⠀⠀⣀⣬⠶⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡘⠒⠦⢤⣀⠀⢠⣿⣯⣛⠀⠀⠀⣶⣶⡗⣉⣁⣸⣉⡿⠮⣍⡒⠦⠄\n        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣾⣾⣾⣿⣷⠉⠁⠀⠀⡀⠀⠀⠀⠀⠀⠉⠀⠈⠀⠀⠀⠀⠱⣄⠀⠀⠀⢩⣿⣿⣿⣿⣷⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣌⣙⠲\n        ⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⣀⣀⣀⣤⣶⣿⣿⣿⣿⣿⣿⣿⣧⣀⠤⠚⠁⠀⠀⠀⠀⠀⠀⣀⣀⣀⣀⠀⠀⠀⠀⠙⠢⣠⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿\n        ⠀⣀⣠⣤⣶⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⡀⠀⠀⠀⠐⠒⠉⢉⣭⣉⠉⠉⠉⠁⠀⠀⢀⡴⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⡀⠀⠀⠀⣤⣼⡟⠈⡇⠀⠀⠀⢀⣴⣯⣶⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣾⣷⣄⣠⡇⡾⠁⠀⡀⠀⣀⡴⣿⣵⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⣿⠟⣸⠀⠀⠀⣿⣿⣿⣷⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠿⠿\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⠧⠭⡽⠀⡇⠀⠀⢠⠇⠐⠒⠒⠒⠒⠚⠛⠉⠀⢀⡤⠄⠀⠀⠠⣄⠀⠀⠀⠈⠉⠉⠉⠁⠀⠀⠀⠀⠀\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⡀⠀⢰⠇⣸⠃⠠⠀⣼⠀⠀⠀⠀⠀⠐⢤⣠⠖⠒⠋⠀⠀⠀⠀⠀⠈⢀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡏⢠⠇⠀⠀⢀⣿⣶⣦⣤⣤⣤⣤⠜⠁⠀⠀⠀⡆⢠⠞⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⢀⡎⠀⠀⠀⣸⣷⣿⣿⣿⡿⣟⠁⠀⠀⢀⣴⡾⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⡞⠀⠀⠀⢀⡟⠉⠉⠀⠀⠀⢈⠆⠀⢴⣿⠏⠀⠀⠀⣰⡶⠖⠚⠛⣓⡦⠤⠤⠤⠤⢀⣠⢤⣂⣀⣀\n        ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣏⡜⠀⠠⡀⢀⣾⡛⠋⠙⣆⣀⣠⠊⠀⠀⣘⡇⠀⠀⣤⠞⣧⣴⡂⣀⢀⡀⢀⣀⡀⣀⣨⣵⣾⣿⣿⣿⣿\n')
        
        
def print_paralel_data(node_name, nnodes, ncore, kpar, memory):
   print(node_name)
   print('------------')
   print('nodes: {}'.format(nnodes))
   print('NCORE: {}'.format(ncore))
   print('KPAR: {}'.format(kpar))
   print('Required Memory = {:.3f} Gb'.format(memory))
   print()


#naruto_start_picture()

nodes = {'regular12':[12, 12, 11.5],'commom20':[20, 4, 63.0],'regular20':[20, 8, 63.0], 'low24':[24, 4, 62.9],
         'big24':[24, 12, 126.0], 'new24':[24, 8, 188.0], 'common32':[32, 8, 188.0], 'regular32':[32, 4, 188.0]}

paws = {'H':'H', 'He':'He', 'Li':'Li_sv', 'Be':'Be', 'B':'B', 'C':'C', 'N':'N', 'O':'O', 'F':'F', 'Ne':'Ne', 'Na':'Na_pv', 'Mg':'Mg',
        'Al':'Al', 'Si':'Si', 'P':'P','S':'S', 'Cl':'Cl', 'Ar':'Ar', 'K':'K_sv', 'Ca':'Ca_sv', 'Sc':'Sc_sv', 'Ti':'Ti_sv', 'V':'V_sv',
        'Cr':'Cr_pv', 'Mn':'Mn_pv', 'Fe':'Fe', 'Co':'Co', 'Ni':'Ni', 'Cu':'Cu', 'Zn':'Zn', 'Ga':'Ga_d', 'Ge':'Ge_d', 'As':'As', 'Se':'Se',
        'Br':'Br', 'Kr':'Kr', 'Rb':'Rb_sv', 'Sr':'Sr_sv', 'Y':'Y_sv', 'Zr':'Zr_sv', 'Nb':'Nb_sv', 'Mo':'Mo_sv', 'Tc':'Tc_pv', 'Ru':'Ru_pv',
        'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'In':'In_d', 'Sn':'Sn_d', 'Sb':'Sb', 'Te':'Te', 'I':'I', 'Xe':'Xe', 'Cs':'Cs_sv',
        'Ba':'Ba_sv', 'La':'La', 'Ce':'Ce', 'Pr':'Pr_3', 'Nd':'Nd_3', 'Pm':'Pm_3', 'Sm':'Sm_3', 'Eu':'Eu_2', 'Gd':'Gd_3', 'Tb':'Tb_3',
        'Dy':'Dy_3', 'Ho':'Ho_3', 'Er':'Er_3', 'Tm':'Tm_3', 'Yb':'Yb_2', 'Lu':'Lu_3', 'Hf':'Hf_pv', 'Ta':'Ta_pv', 'W':'W_sv', 'Re':'Re',
        'Os':'Os', 'Ir':'Ir', 'Pt':'Pt', 'Au':'Au', 'Hg':'Hg', 'Tl':'Tl_d', 'Pb':'Pb_d', 'Bi':'Bi_d', 'Po':'Po_d', 'At':'At', 'Rn':'Rn',
        'Fr':'Fr_sv', 'Ra':'Ra_sv', 'Ac':'Ac', 'Th':'Th', 'Pa':'Pa', 'U':'U', 'Np':'Np', 'Pu':'Pu', 'Am':'Am', 'Cm':'Cm'}



structure = Structure.from_file('POSCAR')
poscar = inputs.Poscar(structure)
incar = inputs.Incar.from_file('INCAR')

if os.path.exists('KPOINTS'):
   kpoints = inputs.Kpoints.from_file('KPOINTS')
else:
   kppra = float(input('Input required KPPRA\n'))
   kpoints = inputs.Kpoints.automatic_density(structure, kppra)
   kpoints.write_file('KPOINTS')
    
if os.path.exists('POTCAR'):
   potcar = inputs.Potcar.from_file('POTCAR') 
   potentials = {p.split('_')[0]:p for p in potcar.symbols}
else:
   potentials = [paws[s.symbol] for s in structure.types_of_species]
   potcar = inputs.Potcar(symbols=potentials, functional="PBE_54")

nelectrons = np.array([int(potcar[i].nelectrons) for i in range(len(potentials))])
encut = max([potcar[i].keywords['ENMAX'] for i in range(len(potentials))])


if incar['PREC'][0].lower() == 'a' or incar['PREC'][0].lower() == 's':
        prec_coefficient = 2
else:
        prec_coefficient = 3/2

rytoev = 13.605826
aB = 0.529177249
cutof = np.array([np.sqrt(encut/rytoev) / (2*pi/(anorm/aB)) for anorm in poscar.structure.lattice.abc])
ng = prec_coefficient * cutof + 0.5

nbands = int(sum(nelectrons*poscar.natoms)*0.6 + sum(poscar.natoms))
nkdim = len(SpacegroupAnalyzer(structure).get_ir_reciprocal_mesh(kpoints.kpts))

memory_per_core = np.prod(ng) * nbands * nkdim * 16 + 4*(ng[0]*prec_coefficient/2 + 1) * ng[1]*prec_coefficient * ng[2]*prec_coefficient * 16 
memory_per_core /= 1024**3

max_nodes = int(input('Input max nodes for task (0 - no limit)\n')) 
for k in nodes.keys():
        nbands_round =  nodes[k][0] * round(nbands/nodes[k][0]) # Округлить до кратного числу cores
        npar_min = nodes[k][0]
        for npar in range(nodes[k][0], 0, -2): # Цикл по npar
           if max_nodes == 0:
              if int(np.ceil(nbands_round/npar/nodes[k][0])) <= nodes[k][1]:
                 npar_min = npar
           else:
              if int(np.ceil(nbands_round/npar/nodes[k][0])) <= max_nodes: 
                 npar_min = npar           

        nnodes = int(np.ceil(nbands_round/npar_min/nodes[k][0]))
        
        memory = memory_per_core * nnodes * nodes[k][0]
        if memory < nnodes * nodes[k][2]:
            print_paralel_data(k, nnodes, int(nodes[k][0]/npar_min), np.gcd(nkdim, nnodes*nodes[k][0]), memory)
        else:
            print('Not enough memory on these nodes')



'''
if nnodes < node_core[node][1]:
    print("NODES : {}\nNCORE = {}\nKPAR = {}\ndon't use KPAR if required memory bigger than memory per one node".format(nnodes, int(node_core[node][0]/npar), nnodes))
else:
    print('Not enough nodes to efficiently parallelize the task\n\t* select a larger NPAR\n\t* use other nodes\n\t* use non-paralleling settings (remove NPAR, NCORE, KPAR tags from INCAR)')
if memory < nnodes * node_core[node][2]:
    print('Required Memory = {:.3f} Gb'.format(memory))
else:
    print('Required Memory = {:.3f} Gb. Not enough memory on the nodes,\n increase the number of nodes to {} and use non-paralleling settings (remove NPAR, NCORE, KPAR tags from INCAR)'.format(memory, int(memory/node_core[node][2])))
'''
