# -*- coding: utf-8 -*-

import numpy as np
import os
from pymatgen.io.vasp.outputs import Vasprun
from scipy.interpolate import splrep, BSpline
import matplotlib.pyplot as plt

min_mu = -0.1 #eV
max_mu = 0.1 #eV

def fermi_dirac(t, e, ef):
    k = 8.617333262e-5 #eV/K
    f = 1/(np.exp((e - ef)/(k*t)) + 1)
    return(f)

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



chempot = np.loadtxt('mu_T.dat')
#fermi_level_dft = Vasprun('vasprun.xml').efermi
fermi_level_dft = 6.95235578
t, c, k = splrep(chempot[:,0], fermi_level_dft+chempot[:,1], s=0.0, k=3)
chempot_function = BSpline(t, c, k, extrapolate=False)


mu_T = np.loadtxt('wannier90.1_seebeck.dat')[:,[0,1]]

kappal_path = '/home/buche/DHH/Korneva/Ti2MnNiSi2/wannie90/Kappal.dat'
#kappal_path = '/run/user/1000/gvfs/sftp:host=195.54.14.162,user=buche/home/buche/DHH/Korneva/Ti2MnNiSi2/wannie90/Kappal.dat'

kappal = np.loadtxt(kappal_path)
i1, i2, i3 = splrep(kappal[:,0], kappal[:,1], s=0.0, k=3)
kappal_interp = BSpline(i1, i2, i3, extrapolate=True)
kappal = np.transpose([kappal_interp(mu_T[:,1])]*3)


T_for_fig = np.array(input('Enter temperature for X(mu) graphs (min value: {} ; max value: {} ; step: {})\n'. format(min(mu_T[:,1]), max(mu_T[:,1]), mu_T[1,1] - mu_T[0,1])).split(), dtype = float)

seebeck_up = np.loadtxt('wannier90.1_seebeck.dat')[:,[2,6,10]]
elcond_up = np.loadtxt('wannier90.1_elcond.dat')[:,[2,4,7]]
kappae_up = np.loadtxt('wannier90.1_kappa.dat')[:,[2,4,7]]

data_up = np.column_stack((mu_T,
                           seebeck_up,
                           elcond_up, 
                           kappae_up,
                           seebeck_up**2 * elcond_up,
                           (seebeck_up**2 * elcond_up * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / kappae_up,
                           (seebeck_up**2 * elcond_up * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / kappal,
                           (seebeck_up**2 * elcond_up * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / (kappae_up + kappal)))
                           
seebeck_dn = np.loadtxt('wannier90.2_seebeck.dat')[:,[2,6,10]]   
elcond_dn = np.loadtxt('wannier90.2_elcond.dat')[:,[2,4,7]]
kappae_dn = np.loadtxt('wannier90.2_kappa.dat')[:,[2,4,7]]
                 
data_dn = np.column_stack((mu_T,
                           seebeck_dn,
                           elcond_dn, 
                           kappae_dn,
                           seebeck_dn**2 * elcond_dn,
                           (seebeck_dn**2 * elcond_dn * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / kappae_dn,
                           (seebeck_dn**2 * elcond_dn * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / kappal,
                           (seebeck_dn**2 * elcond_dn * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / (kappae_dn + kappal)))

seebeck_tot = (seebeck_up*elcond_up + seebeck_dn*elcond_dn) / (elcond_up + elcond_dn)
elcond_tot = elcond_up + elcond_dn
kappae_tot = kappae_up + kappae_dn

data_tot = np.column_stack((mu_T,
                           seebeck_tot,
                           elcond_tot, 
                           kappae_tot,
                           seebeck_tot**2 * elcond_tot,
                            (seebeck_tot**2 * elcond_tot * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / kappae_tot,
                            (seebeck_tot**2 * elcond_tot * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / kappal,
                           (seebeck_tot**2 * elcond_tot * np.column_stack((mu_T[:,1], mu_T[:,1], mu_T[:,1]))) / (kappae_tot + kappal)))

if not os.path.isdir('boltzwann_data/from_T'):
    os.makedirs('boltzwann_data/from_T')
if not os.path.isdir('boltzwann_data/from_mu'):
    os.makedirs('boltzwann_data/from_mu')

#comments = ['S [V/K]', '$\sigma$ [1/(ohm*m)]', '$\kappa_e$ [W/(m*K)]', 'PF (W/[K m$^2$])', 'ZT', '$\sigma$/$\kappa_e$ [K/(Om*W)]']
comments = ['S [V/K]', '$\sigma$ [1/(ohm*m)]', '$\kappa_e$ [W/(m*K)]', 'PF (W/[K m$^2$])', 'ZT$_e$', 'ZT$_l$', 'ZT']
names = ['S', 'sigma', 'kappae', 'PF', 'ZT_e', 'ZT_l', 'ZT']




num = np.arange(2, data_tot.shape[1]+3, 3)
markers = ['o', '^', 's']
linestyles = ['-', '--', ':']
label_proj = ['xx', 'yy', 'zz']
for i in range(len(names)):

    print(names[i])
    
    fig, axs = plt.subplots(1, 3, figsize = tuple(np.array([15, 5])))
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    axs[0].set_ylabel("{}".format(comments[i]), size=18, labelpad = 0.0)
    for j in range(3):
        axs[j].set_xlabel("Temperature (K)", size=18, labelpad = -1.0)  
    axs[0].set_title('spin-up')
    axs[1].set_title('spin-down')
    axs[2].set_title('spin-tot')        
    
    fig_t, axs_t = plt.subplots(1, 3, figsize = tuple(np.array([15, 5])))
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    axs_t[0].set_ylabel("{}".format(comments[i]), size=18, labelpad = 0.0)
    axs_t[0].set_title('spin-up')
    axs_t[1].set_title('spin-down')
    axs_t[2].set_title('spin-tot')  
    
    for num_proj in range(3):
        data_t = []
        for t in np.unique(data_up[:,1]):
            a1, b1, c1 = splrep(data_up[data_up[:,1]==t, 0], data_up[data_up[:,1]==t, num[i]+num_proj], s=0.0, k=3)
            data_function_1 = BSpline(a1, b1, c1, extrapolate=False)
            a2, b2, c2 = splrep(data_dn[data_dn[:,1]==t, 0], data_dn[data_dn[:,1]==t, num[i]+num_proj], s=0.0, k=3)
            data_function_2 = BSpline(a2, b2, c2, extrapolate=False)
            a3, b3, c3 = splrep(data_tot[data_tot[:,1]==t, 0], data_tot[data_tot[:,1]==t, num[i]+num_proj], s=0.0, k=3)
            data_function_3 = BSpline(a3, b3, c3, extrapolate=False)
        
            mu_interp = np.linspace(min_mu+chempot_function(t), max_mu+chempot_function(t), 1000)
 
            if t in T_for_fig:          
                curent_t = t
                for j in range(3):
                    axs_t[j].set_xlabel("Chemical potential (eV)", size=18, labelpad = -1.0) 

                axs_t[0].plot(mu_interp- chempot_function(t), data_function_1(mu_interp), color = '#f26430', marker = 'None', markersize = 9,  linestyle = linestyles[num_proj], linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#f26430', markeredgecolor = '#f26430', label = label_proj[num_proj])
                axs_t[1].plot(mu_interp- chempot_function(t), data_function_2(mu_interp), color = '#009ddc', marker = 'None', markersize = 9, linestyle = linestyles[num_proj], linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#009ddc', markeredgecolor = '#009ddc', label = label_proj[num_proj])
                axs_t[2].plot(mu_interp- chempot_function(t), data_function_3(mu_interp), color = '#2a2d34', marker = 'None', markersize = 9, linestyle = linestyles[num_proj], linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#2a2d34', markeredgecolor = '#2a2d34', label = label_proj[num_proj])
  
                axs_t[0].plot([0, 0], [min(data_function_1(mu_interp)), max(data_function_1(mu_interp))], color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '--', linewidth=1, alpha = 1)  
                axs_t[1].plot([0, 0], [min(data_function_2(mu_interp)), max(data_function_2(mu_interp))], color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '--', linewidth=1, alpha = 1)  
                axs_t[2].plot([0, 0], [min(data_function_3(mu_interp)), max(data_function_3(mu_interp))], color = '#2a2d34', marker = 'None', markersize = 9, linestyle = '--', linewidth=1, alpha = 1)  
          
                for j in range(3):
                    format_ax(axs_t[j])
                
            data_t += [[t, data_function_1(chempot_function(t)), data_function_2(chempot_function(t)), data_function_3(chempot_function(t))]]

        data_t = np.array(data_t)
        np.savetxt('boltzwann_data/from_T/{}_{}_T.dat'.format(names[i], num_proj), data_t)  


        axs[0].plot(data_t[:,0], data_t[:,1], color = '#f26430', marker = markers[num_proj], markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#f26430', markeredgecolor = '#f26430', label = label_proj[num_proj])
        axs[1].plot(data_t[:,0], data_t[:,2], color = '#009ddc', marker = markers[num_proj], markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#009ddc', markeredgecolor = '#009ddc', label = label_proj[num_proj])
        axs[2].plot(data_t[:,0], data_t[:,3], color = '#2a2d34', marker = markers[num_proj], markersize = 9, linestyle = '-', linewidth=3, alpha = 0.5, markeredgewidth=1, markerfacecolor = '#2a2d34', markeredgecolor = '#2a2d34', label = label_proj[num_proj])

        for j in range(3):
            format_ax(axs[j])

    fig_t.savefig('boltzwann_data/from_mu/{}_mu_{}.png'.format(names[i], curent_t), dpi=300, transparent=False, bbox_inches = 'tight')
    fig.savefig('boltzwann_data/from_T/{}_T.png'.format(names[i]), dpi=300, transparent=False, bbox_inches = 'tight')


#Отдельно для DOS
print("DOS")

dos_up = np.loadtxt('wannier90.1_boltzdos.dat')
min_index = np.where(dos_up[:,1]>0)[0][0]
max_index = np.where(dos_up[min_index:,1]==0)[0][0] + min_index
dos_up = dos_up[min_index-5:max_index+5, :]

dos_dn = np.loadtxt('wannier90.2_boltzdos.dat')
min_index = np.where(dos_dn[:,1]>0)[0][0]
max_index = np.where(dos_dn[min_index:,1]==0)[0][0] + min_index
dos_dn = dos_dn[min_index-5:max_index+5, :]


for t in T_for_fig:

    np.savetxt('boltzwann_data/from_mu/DOS_mu_{}.dat'.format(t), np.column_stack((dos_up[:,0] - chempot_function(t), dos_up[:,1]*fermi_dirac(t,dos_up[:,0],chempot_function(t)))))

    fig, axs = plt.subplots(1, 2, figsize=tuple(np.array([15, 5])))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    axs[0].set_ylabel("DOS", size=18, labelpad=0.0)

    for j in range(2):
        axs[j].set_xlabel("$E-E_F$ (eV)", size=18, labelpad=-1.0)

    axs[0].plot(dos_up[:,0] - chempot_function(t), dos_up[:,1]*fermi_dirac(t,dos_up[:,0],chempot_function(t)), color='#f26430', marker='None',
                markersize=9, linestyle='-', linewidth=3, alpha=0.5, markeredgewidth=1, markerfacecolor='#f26430',
                markeredgecolor='#f26430', label='spin-up')
    axs[1].plot(dos_dn[:,0] - chempot_function(t), dos_dn[:,1]*fermi_dirac(t,dos_dn[:,0],chempot_function(t)), color='#009ddc', marker='None',
                markersize=9, linestyle='-', linewidth=3, alpha=0.5, markeredgewidth=1, markerfacecolor='#009ddc',
                markeredgecolor='#009ddc', label='spin-down')

    for j in range(2):
        format_ax(axs[j])

    fig.savefig('boltzwann_data/from_mu/DOS_mu_{}.png'.format(t), dpi=300, transparent=False,
                bbox_inches='tight')
    plt.close()

#Отдельно для kappal
print("kappal")

np.savetxt('boltzwann_data/from_T/kappal_T.dat', np.column_stack((np.unique(mu_T[:,1]), kappal[0:len(np.unique(mu_T[:,1])),0])))

fig, ax = plt.subplots(figsize=tuple(np.array([15, 5])))
plt.subplots_adjust(hspace=0.3, wspace=0.3)
ax.set_ylabel("$\kappa_e$ [W/(m*K)]", size=18, labelpad=0.0)
ax.set_xlabel("Temperature (K)", size=18, labelpad=-1.0)

ax.plot(np.unique(mu_T[:,1]), kappal[0:len(np.unique(mu_T[:,1])),0], color='#f26430', marker='None',
                markersize=9, linestyle='-', linewidth=3, alpha=0.5, markeredgewidth=1, markerfacecolor='#f26430',
                markeredgecolor='#f26430', label='spin-up')

format_ax(ax)

fig.savefig('boltzwann_data/from_T/kappal_T.png', dpi=300, transparent=False,
                bbox_inches='tight')
plt.close()
