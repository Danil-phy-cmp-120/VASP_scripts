import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import  Element
import pymatgen.io.vasp.inputs as inputs 
import os

paws = {'H':'H', 'He':'He', 'Li':'Li_sv', 'Be':'Be', 'B':'B', 'C':'C', 'N':'N', 'O':'O', 'F':'F', 'Ne':'Ne', 'Na':'Na_pv', 'Mg':'Mg',
        'Al':'Al', 'Si':'Si', 'P':'P','S':'S', 'Cl':'Cl', 'Ar':'Ar', 'K':'K_sv', 'Ca':'Ca_sv', 'Sc':'Sc_sv', 'Ti':'Ti_sv', 'V':'V_sv',
        'Cr':'Cr_pv', 'Mn':'Mn_sv', 'Fe':'Fe', 'Co':'Co', 'Ni':'Ni', 'Cu':'Cu', 'Zn':'Zn', 'Ga':'Ga_d', 'Ge':'Ge_d', 'As':'As_d', 'Se':'Se',
        'Br':'Br', 'Kr':'Kr', 'Rb':'Rb_sv', 'Sr':'Sr_sv', 'Y':'Y_sv', 'Zr':'Zr_sv', 'Nb':'Nb_sv', 'Mo':'Mo_sv', 'Tc':'Tc_pv', 'Ru':'Ru_pv',
        'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'In':'In_d', 'Sn':'Sn_d', 'Sb':'Sb', 'Te':'Te', 'I':'I', 'Xe':'Xe', 'Cs':'Cs_sv',
        'Ba':'Ba_sv', 'La':'La', 'Ce':'Ce', 'Pr':'Pr_3', 'Nd':'Nd_3', 'Pm':'Pm_3', 'Sm':'Sm_3', 'Eu':'Eu_2', 'Gd':'Gd_3', 'Tb':'Tb_3',
        'Dy':'Dy_3', 'Ho':'Ho_3', 'Er':'Er_3', 'Tm':'Tm_3', 'Yb':'Yb_2', 'Lu':'Lu_3', 'Hf':'Hf_pv', 'Ta':'Ta_pv', 'W':'W_sv', 'Re':'Re',
        'Os':'Os', 'Ir':'Ir', 'Pt':'Pt', 'Au':'Au', 'Hg':'Hg', 'Tl':'Tl_d', 'Pb':'Pb_d', 'Bi':'Bi_d', 'Po':'Po_d', 'At':'At', 'Rn':'Rn',
        'Fr':'Fr_sv', 'Ra':'Ra_sv', 'Ac':'Ac', 'Th':'Th', 'Pa':'Pa', 'U':'U', 'Np':'Np', 'Pu':'Pu', 'Am':'Am', 'Cm':'Cm'}


path = input('Enter the path to CONTCAR\n')

for j in np.arange(0.8, 1.425, 0.025):       
    if not os.path.isdir('ca/{:.3f}'.format(j)):
        os.makedirs('ca/{:.3f}'.format(j))
            
    structure = Structure.from_file('{}/CONTCAR'.format(path))
            
    D = np.array([[(j)**(-1.0/3.0), 0, 0],
                  [0, (j)**(-1.0/3.0), 0],
                  [0, 0, (j)**(2.0/3.0)]])

    structure.lattice = np.dot(structure.lattice.matrix, D)
    poscar = inputs.Poscar(structure)                   
                        
    potentials = [paws[s.symbol] for s in structure.types_of_species]    
    potcar = inputs.Potcar(symbols=paw, functional="PBE_54")

    kpoints = inputs.Kpoints.automatic_density(structure, 4000)
    
    incar = inputs.Incar.from_file('{}/INCAR'.format(path))
    incar['ISIF'] = 7
    incar['NSW'] = 50

    '''incar = inputs.Incar({ 'METAGGA': 'SCAN',
                           'LCHARG' : False,
                           'LWAVE' :  False,
                           'ICHARG' : 1,
                           'IBRION' : 2,
                           'ISIF' : 2,
                           'POTIM' : 0.1,
                           'PREC' : 'Accurate',
                           'ENCUT' : 500,
                           'EDIFF' : 1e-7,
                           'EDIFFG' : 1e-3,
                           'ISMEAR' : 1,
                           'SIGMA' : 0.1,
                           'ISPIN' : 2,
                           'LMAXMIX' : 6,
                           'LMIXTAU' : True,
                           'ADDGRID' : True,
                           'LORBIT' : 11,
                           'LASPH' : True,
                           'NSW ': 100,
                           'NEDOS' : 1000,
                           'NELM' : 101,
                           'ALGO' : 'Normal',
                           'NCORE' : 6,    
                           'KPAR' : 2,
                           'MAGMOM': '8*4 4*0.5 4*-0.1'})'''
                       
                       
    VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
    VaspInput.write_input('ca/{:.3f}'.format(j))
