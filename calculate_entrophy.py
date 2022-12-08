import numpy as np
import os
from scipy.interpolate import splrep, splev, splint, CubicSpline, BSpline
from scipy.integrate import cumtrapz
from scipy.misc import derivative
import matplotlib.pyplot as plt

energy = np.loadtxt('E_T.dat')
energy[:,1] /= 1000


t, c, k = splrep(energy[:,0], energy[:,1], s=0.0035, k=3)

#interp_function = CubicSpline(energy[:,0], energy[:,1], bc_type='not-a-knot')
interp_function = BSpline(t, c, k, extrapolate=False)
energy_new = np.column_stack( (energy[:,0], interp_function(energy[:,0])) )

energy_derivative = interp_function.derivative(1)


heat_capacity = np.column_stack( (energy_new[:,0], energy_derivative(energy_new[:,0])) )
c_t = np.column_stack( (heat_capacity[:,0], heat_capacity[:,1]/energy[:,0]) )

entrophy = np.column_stack( (c_t[:,0], cumtrapz(c_t[:,1], c_t[:,0], initial = 0)) )

np.savetxt('out.dat', np.column_stack( (energy_new, heat_capacity[:,1], entrophy[:,1], entrophy[:,1]*energy_new[:,0]) ), fmt = '%.0f\t%.18e\t%.18e\t%.18e\t%.18e', header='T (K)\tE (kJ/mol)\tC (kJ/(K*mol))\tS (kJ/(K*mol))\tS*T (kJ/mol)')

plt.plot(energy[:,0], energy[:,1], energy_new[:,0], energy_new[:,1] )
plt.savefig('test_interp.png')
plt.close()

plt.plot(heat_capacity[:,0], heat_capacity[:,1])
plt.savefig('heat_capacity.png')
plt.close()
