# -*- coding: utf-8 -*-

#################### V3 ####################

import numpy as np
import os
from pymatgen.io.vasp.outputs import Outcar, Vasprun, BSVasprun, Procar
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from scipy.interpolate import splrep, BSpline, interp2d
from scipy.integrate import cumtrapz
from scipy.misc import derivative
import matplotlib.pyplot as plt

min_mu = -0.05 #eV
max_mu = 0.15 #eV

flag_plot = bool(input('Do you want to make graphs? (True/False)\n'))

def format_ax(ax):
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
    #ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    #ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))


    ax.grid(which='major',
            color = 'lightgray',
            linestyle = ':',
            linewidth = 1)


    ax.tick_params(axis = 'both',    #  Применяем параметры к обеим осям
                   which = 'major',    #  Применяем параметры к вспомогательным делениям
                   direction = 'in',    #  Рисуем деления внутри и снаружи графика
                   #length = 10,    #  Длинна делений
                   #width = 2,     #  Ширина делений
                   #color = 'm',    #  Цвет делений
                   #pad = 10,    #  Расстояние между черточкой и ее подписью
                   labelsize = 14,    #  Размер подписи
                   #labelcolor = 'r',    #  Цвет подписи
                   bottom = True,    #  Рисуем метки снизу
                   top = True,    #   сверху
                   left = True,    #  слева
                   right = True)    #  и справа


    #ax.set_ylim(-120, 20)
    #ax.set_xlim(85, 615)

    legend= ax.legend(fontsize = 11,
                ncol = 1,    #  количество столбцов
                loc='best',
                #bbox_to_anchor=(0.10, 0.99),
                facecolor = 'white',    #  цвет области
                framealpha = 1,
                #edgecolor = 'None',    #  цвет крайней линии
                #title = ' ',    #  заголовок
                #title_fontsize = 20   #  размер шрифта заголовка
                )
    legend.get_title().set_fontsize('14')
    
    return(ax)


log = open('log_tot.out').read()
log = log.replace('INFO    │ ', '').replace('\n', ' ').replace('[', '').replace(']', '')
log_split = log.split()

mu_T = []
flag_T = False
flag_mu = False
for l in log_split:
    if l == 'Fermi':
        flag_T = False 
    if l == 'Chemical':
        flag_mu = False 

    if flag_T == True:
        mu_T += [float(l)]
    if flag_mu == True:
        mu_T += [float(l)*2]

    if l == 'Temperatures:':
        flag_T = True 
    if l == 'T:':
        flag_mu = True 

        
mu_T = np.array(mu_T)
mu_T = np.column_stack( (mu_T[0:int(mu_T.size/2)], mu_T[int(mu_T.size/2)::]) )
mu_T[:,1] *= 13.605693122994
print(mu_T)
np.savetxt('mu_T.dat', mu_T)

t, c, k = splrep(mu_T[:,0], mu_T[:,1], s=0.0, k=3)
chempot_function = BSpline(t, c, k, extrapolate=True)


for p in os.listdir(os.getcwd()):
    if p.find('_up.trace') != -1:
        trace_up = np.loadtxt(p)
        trace_up[:,0] *= 13.605693122994
        trace_up[:,3] /= 2*13.605693122994
        trace_up = np.column_stack( (trace_up, trace_up[:,4]**2 * trace_up[:,5], np.nan_to_num((trace_up[:,4]**2 * trace_up[:,5] * trace_up[:,1])/trace_up[:,7]), trace_up[:,5]/trace_up[:,7]) )
        trace_up = np.array(sorted(trace_up, key=lambda x: (x[1])))
    elif p.find('_dw.trace') != -1:
        trace_dn = np.loadtxt(p)
        trace_dn[:,0] *= 13.605693122994
        trace_dn[:,3] /= 2*13.605693122994
        trace_dn = np.column_stack( (trace_dn, trace_dn[:,4]**2 * trace_dn[:,5], np.nan_to_num((trace_dn[:,4]**2 * trace_dn[:,5] * trace_dn[:,1])/trace_dn[:,7]), trace_dn[:,5]/trace_dn[:,7]) )
        trace_dn = np.array(sorted(trace_dn, key=lambda x: (x[1])))
    elif p.find('_tot.trace') != -1:
        trace_tot = np.loadtxt(p)
        trace_tot[:,0] *= 13.605693122994
        trace_tot[:,3] /= 2*13.605693122994
        trace_tot = np.column_stack( (trace_tot, trace_tot[:,4]**2 * trace_tot[:,5], np.nan_to_num((trace_tot[:,4]**2 * trace_tot[:,5] * trace_tot[:,1])/trace_tot[:,7]), trace_tot[:,5]/trace_tot[:,7]) )
        trace_tot = np.array(sorted(trace_tot, key=lambda x: (x[1])))

T_for_fig = np.array(input('Enter temperature for X(mu) graphs (min value: {} ; max value: {} ; step: {})\n'. format(min(trace_up[:,1]), max(trace_up[:,1]), trace_up[1,1] - trace_up[0,1])).split(), dtype = float)

if not os.path.isdir('data/from_T'):
    os.makedirs('data/from_T')
if not os.path.isdir('data/from_mu'):
    os.makedirs('data/from_mu')

comments = ['N [e/uc]', 'DOS($E_F$) [1/(eV*uc)]', 'S [V/K]', '$\sigma$/$\\tau_0$ [1/(ohm*m*s)]', 'RH [m$^3$/C]', '$\kappa_e$/$\\tau_0$ [W/(m*K*s)]', 'cv [J/(mol*K)]', '$\chi$ [m$^3$/mol]', 'PF (W/[K m$^2$])', 'ZT', '$\sigma$/$\kappa_e$ [K/(Om*W)]']
names = ['N', 'DOS', 'S', 'sigma', 'RH', 'kappae', 'cv', 'chi', 'PF', 'ZT', 'sigma_per_kappae']

for i in range(len(names)):
    print(names[i])
    trace_t = []
    for t in np.unique(trace_up[:,1]):
        
        a1, b1, c1 = splrep(trace_up[trace_up[:,1]==t, 0], trace_up[trace_up[:,1]==t, i+2], s=0.0, k=3)
        trace_function_1 = BSpline(a1, b1, c1, extrapolate=True)
        a2, b2, c2 = splrep(trace_dn[trace_dn[:,1]==t, 0], trace_dn[trace_dn[:,1]==t, i+2], s=0.0, k=3)
        trace_function_2 = BSpline(a2, b2, c2, extrapolate=True)
        a3, b3, c3 = splrep(trace_tot[trace_tot[:,1]==t, 0], trace_tot[trace_tot[:,1]==t, i+2], s=0.0, k=3)
        trace_function_3 = BSpline(a3, b3, c3, extrapolate=True)
    
        trace_t += [[t, trace_function_1(chempot_function(t)), trace_function_2(chempot_function(t)), trace_function_3(chempot_function(t))]]
        mu_interp = np.linspace(min_mu+chempot_function(t), max_mu+chempot_function(t), 5000)
        
        np.savetxt('data/from_mu/{}_mu_{}.dat'.format(names[i], t), np.column_stack((mu_interp- chempot_function(t), trace_function_1(mu_interp), trace_function_2(mu_interp), trace_function_3(mu_interp))))

        if (flag_plot == True) and (t in T_for_fig):      
            fig, axs = plt.subplots(1, 3, figsize = tuple(np.array([15, 5])))
            plt.subplots_adjust(hspace=0.3,wspace=0.3)
            axs[0].set_ylabel("{}".format(comments[i]), size=18, labelpad = 0.0)
        
            for j in range(3):
                axs[j].set_xlabel("Chemical potential (eV)", size=18, labelpad = -1.0) 
        
            axs[0].plot(mu_interp- chempot_function(t), trace_function_1(mu_interp), color = '#f26430', marker = 'None', markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#f26430', markeredgecolor = '#f26430', label = 'spin-up')
            axs[1].plot(mu_interp- chempot_function(t), trace_function_2(mu_interp), color = '#009ddc', marker = 'None', markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#009ddc', markeredgecolor = '#009ddc', label = 'spin-down')
            axs[2].plot(mu_interp- chempot_function(t), trace_function_3(mu_interp), color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#2a2d34', markeredgecolor = '#2a2d34', label = 'spin-tot')
  
            axs[0].plot([0, 0], [min(trace_function_1(mu_interp)), max(trace_function_1(mu_interp))], color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '--', linewidth=1, alpha = 1)  
            axs[1].plot([0, 0], [min(trace_function_2(mu_interp)), max(trace_function_2(mu_interp))], color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '--', linewidth=1, alpha = 1)  
            axs[2].plot([0, 0], [min(trace_function_3(mu_interp)), max(trace_function_3(mu_interp))], color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '--', linewidth=1, alpha = 1)  
          
            for j in range(3):
                format_ax(axs[j])
        
            fig.savefig('data/from_mu/{}_mu_{}.png'.format(names[i], t), dpi=300, transparent=False, bbox_inches = 'tight')
            plt.close()
        
    trace_t = np.array(trace_t)
    np.savetxt('data/from_T/{}_T.dat'.format(names[i]), trace_t)   


    if flag_plot == True: 
        fig, axs = plt.subplots(1, 3, figsize = tuple(np.array([15, 5])))
        plt.subplots_adjust(hspace=0.3,wspace=0.3)
        axs[0].set_ylabel("{}".format(comments[i]), size=18, labelpad = 0.0)
        for j in range(3):
            axs[j].set_xlabel("Temperature (K)", size=18, labelpad = -1.0)    


        axs[0].plot(trace_t[:,0], trace_t[:,1], color = '#f26430', marker = 'o', markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#f26430', markeredgecolor = '#f26430', label = 'spin-up')
        axs[1].plot(trace_t[:,0], trace_t[:,2], color = '#009ddc', marker = 'o', markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#009ddc', markeredgecolor = '#009ddc', label = 'spin-down')
        axs[2].plot(trace_t[:,0], trace_t[:,3], color = '#2a2d34', marker = 'o', markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#2a2d34', markeredgecolor = '#2a2d34', label = 'spin-tot')

        for j in range(3):
            format_ax(axs[j])

        fig.savefig('data/from_T/{}_T.png'.format(names[i]), dpi=300, transparent=False, bbox_inches = 'tight')
