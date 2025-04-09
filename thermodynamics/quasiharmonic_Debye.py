import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.quasiharmonic import QuasiharmonicDebyeApprox
import os
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
from scipy.constants import physical_constants


t_min = 100
t_max = 1000
t_step = 10

mu = float(input('Enter Poisson ratio\n'))
path = input('Enter path to E(V) calculation\n')
    
dirs = sorted([i for i in os.listdir(path)])

data = []
for d in dirs:
    outcar = Outcar('{}/{}/OUTCAR'.format(path, d))
    structure = Structure.from_file('{}/{}/CONTCAR'.format(path, d))
    data += [[structure.lattice.volume, outcar.final_energy]]
data = np.array(data)
    
spl = splrep(data[:, 0], data[:, 1])
data_interpolate = np.column_stack((np.linspace(min(data[:, 0]), max(data[:, 0]), 100), splev(np.linspace(min(data[:, 0]), max(data[:, 0]), 100), spl)))

structure_0 = Structure.from_file('{}/1.00/CONTCAR'.format(path))
mass = sum(e.atomic_mass for e in structure_0.species)
natoms = structure_0.composition.num_atoms
avg_mass = physical_constants["atomic mass constant"][0] * mass / natoms
    
qda = QuasiharmonicDebyeApprox(data_interpolate[:, 1], data_interpolate[:, 0], structure_0, poisson=mu, t_min=t_min, t_step=t_step,
                                   t_max=t_max, eos = 'murnaghan')
        
gamma = np.mean(qda.gruneisen_parameter(0, data_interpolate[:, 0]))
print(gamma)
        
theta_d = qda.debye_temperature(structure_0.lattice.volume)
theta_a = theta_d * natoms ** (-1.0 / 3.0)   
        
prefactor = (0.849 * 3 * 4 ** (1 / 3)) / (20.0 * np.pi**3)
prefactor = prefactor * (physical_constants["Boltzmann constant in eV/K"][0] / physical_constants["Planck constant over 2 pi in eV s"][0]) ** 3  * avg_mass
kappa = prefactor / (gamma**2 - 0.514 * gamma + 0.228) 
kappa = kappa * theta_a**2 * structure_0.lattice.volume ** (1 / 3) * 1e-10 

thermal_conductivity = []
for t in np.arange(t_min, t_max + t_step, t_step):       
   thermal_conductivity += [[t, kappa * (theta_a / t)]]
thermal_conductivity = np.array(thermal_conductivity)

np.savetxt('thermal_conductivity.dat'), thermal_conductivity)
