# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants 
from scipy import integrate
from scipy.misc import derivative

######### Функция Ферми #########
def f_Fermi(e, mu, T):
    return 1.0/(np.exp(((e-mu)*constants.e) / (constants.k*T)) + 1.0)

######### Считывание энергии Ферми из OUTCAR #########
f = open('OUTCAR',"r")
outcar = f.readlines()
f.close()

for line in outcar:
    inp = line.split()
    if len(inp) > 3 and inp[0] == 'E-fermi' and inp[1] == ':':
        E_Fermi = float(inp[2])
    if len(inp) > 5 and inp[0] == 'number' and inp[1] == 'of' and inp[2] == 'dos' and inp[3] == 'NEDOS':
        NEDOS = int('{:.0f}'.format(float(inp[5])))
        num = int('{:.0f}'.format(float(inp[11])))

#####  Считывание DOSCAR #####
#doscar = open('DOS_cut',"r")
doscar = open('DOSCAR',"r")
lines = doscar.readlines()
doscar.close()

totaldos = []
for i in range(len(lines)):
    if i > 5 and i < (NEDOS + 6):
        inp = lines[i].split()
        totaldos += [float(inp[0]), float(inp[1]), float(inp[2])]

totaldos = np.array(totaldos)
totaldos.shape = (int(totaldos.size/3), 3)

num_prim = float(input('Enter the number of atoms in primitive cell:\n'))
T_min = float(input('Enter the minimum temperature:\n'))
T_max = float(input('Enter the maximum temperature:\n'))
T_step = float(input('Enter the temperature step:\n'))
T = np.arange(T_min, T_max + T_step, T_step)

S_el = np.zeros((T.size, 4))
#### ПЕРВЫЙ СПОСОБ ВЫЧИСЛЕНИЯ ####
for i in range(T.size):
    #####  Вычисление подинтегрального выражения #####
    func = f_Fermi(totaldos[:,0], E_Fermi, T[i])
    Interant = np.zeros(totaldos.shape[0])
    for j in range(func.size):
        if func[j] == 1.0 or func[j] == 0.0:
            Interant[j] =  0.0
        else:
            Interant[j] =  - ((totaldos[j,1]+totaldos[j,2])/(num/num_prim)) * ( (1-func[j])*np.log(1 - func[j]) + func[j]*np.log(func[j]) )

    #####  Вычисление интегралла #####
    S_el[i,2] = integrate.simps(Interant, totaldos[:,0]) # В kB простое интегирование
    S_el[i,0] = S_el[i,2] * (constants.k/constants.e) * 1000 # В meV простое интегирование

#### ВТОРОЙ СПОСОБ ВЫЧИСЛЕНИЯ (Зоммерфельд) ####
val, idx = min((val, idx) for (idx, val) in enumerate(abs(totaldos[:,0]-E_Fermi))) # Поиск DOS ближе всего к уровню Ферми
for i in range(T.size):
    S_el[i,3] = ((constants.pi**2.0)/3.0) * T[i] * ((totaldos[idx,1] + totaldos[idx,2]) / ((constants.e/constants.k) * (num/num_prim))) # В kB Зоммерфельд
    S_el[i,1] = S_el[i,3] * (constants.k/constants.e) * 1000 # В meV Зоммерфельд

np.savetxt('S_el.dat', np.column_stack((T, S_el))) # Колонки: температура, энтропия в meV (интегрирование), энтропия в meV (приближенный способ), энтропия в kB (интегрирование), энтропия в kB (приближенный способ)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r"$T$ (K)", size=20)
ax.set_ylabel(r"$S_{el}$ (meV)", size=20)
ax.plot(T, S_el[:,0], label = 'integrall')
ax.plot(T, S_el[:,1], label = 'approx')
fig.savefig('S_el.png', dpi=300)
