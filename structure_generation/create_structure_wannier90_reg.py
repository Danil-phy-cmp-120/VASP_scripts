import numpy as np
from pymatgen.core.structure import Structure
import pymatgen.io.vasp.inputs as inputs
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os

z_atoms = ['Al', 'Si', 'Ga', 'Ge', 'In', 'Sn']
z_potentials = ['Al', 'Si', 'Ga_d', 'Ge_d', 'In_d', 'Sn']

for i in range(len(z_atoms)):
    
    if not os.path.isdir('wannier90/regular/{}'.format(z_atoms[i])):
        os.makedirs('wannier90/regular/{}'.format(z_atoms[i]))

    structure = Structure.from_file('ION/isif3/regular/{}/FIM1/CONTCAR'.format(z_atoms[i]))
    structure = SpacegroupAnalyzer(structure).find_primitive()

    poscar = inputs.Poscar(structure)
                        
    potcar = inputs.Potcar(symbols=['Rh_pv', 'Fe', z_potentials[i]], functional="PBE_54")

    kpoints = inputs.Kpoints(kpts = [[7, 7, 7]])

    incar = inputs.Incar({'LCHARG': True,
                          'LWAVE': True,
                          'ICHARG': 1,
                          'IBRION': 2,
                          'ISIF': 3,
                          'POTIM': 0.1,
                          'PREC': 'Accurate',
                          'ENCUT': 460,
                          'EDIFF': 1e-7,
                          'EDIFFG': 1e-3,
                          'ISMEAR': 1,
                          'SIGMA': 0.1,
                          'ISPIN': 2,
                          'LMAXMIX': 6,
                          'ADDGRID': True,
                          'LORBIT': 10,
                          'LASPH': True,
                          'NSW ': 0,
                          'NEDOS': 1000,
                          'NELM': 101,
                          'ALGO': 'Fast',
                          'MAGMOM': '2*0.5 3.0 -0.1',
                          
                          'LWANNIER90': True,
                          'LWRITE_UNK': False,
                          'LWRITE_MMN_AMN': True,
                          
                          'NUM_WANN': 27,
                          
                          'WANNIER90_WIN': '"#dis_win_min       = -10.0d0\n#dis_win_max       = 11.0d0\ndis_froz_min      = -5.0d0\ndis_froz_max      =  2.0d0\n\ndis_num_iter      =  500\nnum_iter          = 200\nwrite_xyz = .true.\nwrite_hr  = .true.\ndis_mix_ratio = 0.5\nbegin projections\nRh:s,p,d\nFe:s,p,d\nend projections\n\nbands_plot        = true\nbegin kpoint_path\nG 0.000 0.000 0.000 X 0.000 0.500 0.500\nX 0.000 0.500 0.500 W 0.250 0.750 0.500\nW 0.250 0.750 0.500 K 0.375 0.750 0.375\nK 0.375 0.750 0.375 G 0.000 0.000 0.000\nG 0.000 0.000 0.000 L 0.500 0.500 0.500\nL 0.500 0.500 0.500 U 0.250 0.625 0.625\nU 0.250 0.625 0.625 W 0.250 0.750 0.500\nW 0.250 0.750 0.500 L 0.500 0.500 0.500\nL 0.500 0.500 0.500 K 0.375 0.750 0.375\nU 0.250 0.625 0.625 X 0.000 0.500 0.500\nend kpoint_path"'
                          })

                       
                       
    VaspInput = inputs.VaspInput(incar, kpoints, poscar, potcar)
    VaspInput.write_input('wannier90/regular/{}'.format(z_atoms[i]))
