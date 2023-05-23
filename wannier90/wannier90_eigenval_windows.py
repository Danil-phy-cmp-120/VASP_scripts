import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun, BSVasprun, Procar
from pymatgen.electronic_structure.core import Spin

def get_eigenvalues_extremum(path_eigenval='vasprun.xml', path_kpoint='KPOINTS'):
    xml = Vasprun('vasprun.xml')
    eigenvalues = xml.get_band_structure(kpoints_filename='KPOINTS').bands
    eigenvalues_extremum_up = np.zeros((eigenvalues[Spin.up].shape[0], 3))
    eigenvalues_extremum_down = np.zeros((eigenvalues[Spin.down].shape[0], 3))
    for i in range(eigenvalues[Spin.up].shape[0]):
        eigenvalues_extremum_up[i, 0] = i
        eigenvalues_extremum_up[i, 1] = min(eigenvalues[Spin.up][i, :])  # -xml.efermi
        eigenvalues_extremum_up[i, 2] = max(eigenvalues[Spin.up][i, :])  # -xml.efermi

        eigenvalues_extremum_down[i, 0] = i
        eigenvalues_extremum_down[i, 1] = min(eigenvalues[Spin.down][i, :])  # -xml.efermi
        eigenvalues_extremum_down[i, 2] = max(eigenvalues[Spin.down][i, :])  # -xml.efermi
    return eigenvalues_extremum_up, eigenvalues_extremum_down

def eigenval(procar,bands_number,ions_number):
    eigenval = []
    m = 0
    n = 0
    i = 8
    while i < len(procar)-1:
        if i == len(procar)//2-1:
            i += 9
            m = 0
            n = 0
        if m == bands_number*ions_number:
            i += 8
            m = 0
            n = 0
        if n == ions_number:
            i += 5
            n = 0
        eigenval += [[float(procar[i].split()[1]), float(procar[i].split()[2]), float(procar[i].split()[3])]]
        m += 1
        i += 1
        n += 1
    eigenval = np.array(eigenval)
    return eigenval[0:len(eigenval)//2, :], eigenval[len(eigenval)//2:len(eigenval), :]

def eigenval_n_type(eigenval,n_type,num_atoms):
    eigenval_n_type = []
    for n in range(n_type):
        i = 0
        while i*sum(num_atoms)+sum(num_atoms[0:n]) <= len(eigenval)-num_atoms[-1]:
            eigenval_n_type += [np.sum(np.array(eigenval[i*sum(num_atoms)+sum(num_atoms[0:n]):i*sum(num_atoms)+sum(num_atoms[0:n+1])]),axis=0)]
            i += 1
    eigenval_n_type = np.array(eigenval_n_type)
    return eigenval_n_type

def get_eigenval_band(eigenval_n_type,bands_number,k_points_number,n_type):
    eigenval_band = np.zeros((n_type*bands_number,3))
    for n in range(n_type):
        for i in range(bands_number):
            for k in range(k_points_number):
                eigenval_band[i+n*bands_number] += eigenval_n_type[int((1/n_type)*n*len(eigenval_n_type)+i+k*bands_number)]
    eigenval_band = np.array(eigenval_band)
    return eigenval_band

def total_band(eigenval_band):
    total_band = np.sum(eigenval_band, axis=1)
    eigenval_band[:, 0] /= total_band
    eigenval_band[:, 1] /= total_band
    eigenval_band[:, 2] /= total_band
    eigenval_band *= 100
    return eigenval_band

#из файла procar выбираем вклады s p d электронов для каждого иона и каждой зоны и записываем новый массив
procar = open('PROCAR').readlines()
k_points_number = int(procar[1].split()[3])
bands_number = int(procar[1].split()[7])
ions_number = int(procar[1].split()[11])

poscar = open('POSCAR').readlines()
n_type = int(len(poscar[5].split()))
num_atoms = np.array(poscar[6].split(),dtype = int)
type_atoms = poscar[5].split()

eigenval_spin_up, eigenval_spin_down = eigenval(procar,bands_number,ions_number)

spin_up = eigenval_n_type(eigenval_spin_up,n_type,num_atoms)
spin_down = eigenval_n_type(eigenval_spin_down,n_type,num_atoms)

band_up, band_down = get_eigenvalues_extremum()

eigenval_band_up = get_eigenval_band(spin_up,bands_number,k_points_number,n_type)
eigenval_band_down = get_eigenval_band(spin_down,bands_number,k_points_number,n_type)
for i in range(n_type):
    eigenval_band_up_atom = eigenval_band_up[i*bands_number:(i+1)*bands_number, :]
    eigenval_band_down_atom = eigenval_band_down[i*bands_number:(i+1)*bands_number, :]
    np.savetxt('wannier90_eigenval_windows_{}_up.dat'.format(type_atoms[i]),
                np.column_stack((band_up,eigenval_band_up_atom)), header='N\tmin\tmax\t\ts\tp\td', fmt='%.0f\t%.3f\t%.3f\t\t%.1f\t%.1f\t%.1f')
    np.savetxt('wannier90_eigenval_windows_{}_down.dat'.format(type_atoms[i]),
                np.column_stack((band_down,eigenval_band_down_atom)), header='N\tmin\tmax\t\ts\tp\td', fmt='%.0f\t%.3f\t%.3f\t\t%.1f\t%.1f\t%.1f')


