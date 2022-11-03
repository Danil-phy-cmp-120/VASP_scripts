# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate, interpolate
import matplotlib.pyplot as plt
import os

kB = 1.380658E-23
ec = 1.60217733E-19

def f_Fermi_Dirac(E, mu, T):
    return 1/(np.exp((E-mu)/(T*kB/ec)) + 1)


dos_up = np.loadtxt('wannier90.1_boltzdos.dat')
dos_down = np.loadtxt('wannier90.2_boltzdos.dat')

dos_up_interp = interpolate.interp1d(dos_up[:,0], dos_up[:,1])
dos_down_interp = interpolate.interp1d(dos_down[:,0], dos_down[:,1])
energy_interp = np.linspace(-15, max(dos_up[:,0]), 1000)


dos = np.column_stack(( energy_interp, (dos_up_interp(energy_interp) + dos_down_interp(energy_interp))/2 ))


vol = float(open('wannier90.1_boltzdos.dat').readlines()[3].split(':')[1])
N_e = 29

vol *= 10e-30
dos[:,1] /= vol
np.savetxt('DOS_boltz_all.dat', dos)


T_min = 0
T_max = 600
T_step = 2

mu_T = np.zeros( (int((T_max-T_min)/T_step)+1, 2) )

int_dos = np.column_stack(  (dos[:,0], integrate.cumtrapz(dos[:,1], dos[:,0], initial=0)) ) 
int_dos_interpolate = interpolate.interp1d(int_dos[:,1], int_dos[:,0])
mu_0 = int_dos_interpolate(N_e/vol) # Ef (T = 0K)
#mu_T[0,0] = 0
#mu_T[0,1] = mu_0

if not os.path.isdir('DOS_T'):
    os.makedirs('DOS_T')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dos[:,0], dos[:,1], marker = 'None', linestyle = '-', linewidth = 2) #*f_Fermi_Dirac(dos[:,0], mu_0, 0)
fig.savefig('DOS_T/0.png', dpi=100, transparent=False)


delta = 1e-12
n = 0
for t in np.arange(T_min, T_max+T_step, T_step):

    mu_T[n, 0] = t
    mu_1 =  mu_0 * 0.9 # left edge 
    mu_2 =  mu_0 * 1.15 # right edge
    while abs(mu_2 - mu_1) > delta:

        #int_dos_1 = np.column_stack(  (dos[:,0], integrate.cumtrapz(dos[:,1]*f_Fermi_Dirac(dos[:,0], mu_1, t), dos[:,0], initial=0)) )
        #int_dos_interpolate_1 = interpolate.interp1d(int_dos_1[:,0], int_dos_1[:,1])    
        int_dos_1 = integrate.simpson(dos[:,1]*f_Fermi_Dirac(dos[:,0], mu_1, t), dos[:,0])
        
        x = (mu_1 + mu_2)/2
        #int_dos_x = np.column_stack(  (dos[:,0], integrate.cumtrapz(dos[:,1]*f_Fermi_Dirac(dos[:,0], x, t), dos[:,0], initial=0)) )
        #int_dos_interpolate_x = interpolate.interp1d(int_dos_x[:,0], int_dos_x[:,1]) 
        int_dos_x = integrate.simpson(dos[:,1]*f_Fermi_Dirac(dos[:,0], x, t), dos[:,0])

        if (int_dos_x - N_e/vol < 0 and int_dos_1 - N_e/vol < 0) or (int_dos_x - N_e/vol > 0 and int_dos_1 - N_e/vol > 0):
            mu_1 = x
        else:
            mu_2 = x
            
    mu_0 = mu_2    
    mu_T[n, 1] = mu_2
    n += 1
'''    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(dos[:,0], dos[:,1]*f_Fermi_Dirac(dos[:,0], mu_1, t), marker = 'None', linestyle = '-', linewidth = 2)
    fig.savefig('DOS_T/{}.png'.format(t), dpi=100, transparent=False)
'''
np.savetxt('mu_T', mu_T, fmt='%.1f\t%.18e')

