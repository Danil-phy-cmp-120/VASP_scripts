# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def Fermi_Dirac(T):
    E = np.linspace(0, 3, 100)
    Ef = 1
    k = 8.617333262e-5	
    f = np.exp((E - Ef)/(k*T)) + 1
    return(np.column_stack((E, 1/f)))
        
    
T = float(input('Enter the temperature (K)\n'))

fig, axs = plt.subplots(1, 2, figsize = tuple(np.array([1.5, 0.5])*9))
plt.subplots_adjust(hspace=0.05,wspace=0.15)
axs[0].set_ylabel("n", size=18, labelpad = 0.0)
axs[0].set_xlabel("$\epsilon$ / $\mu$", size=18, labelpad = -1.0)
axs[1].set_xlabel("$\epsilon$ / $\mu$", size=18, labelpad = -1.0)

data = Fermi_Dirac(T)
axs[0].plot(data[:,0], data[:,1], marker = 'None', linewidth=3, alpha = 1)
axs[1].plot(data[1:,0], np.diff(data[:,1])/np.diff(data[:,0]), marker = 'None', linewidth=3, alpha = 1)

for ax in axs:
    ax.tick_params(axis = 'both',    #  РџСЂРёРјРµРЅСЏРµРј РїР°СЂР°РјРµС‚СЂС‹ Рє РѕР±РµРёРј РѕСЃСЏРј
                   which = 'major',    #  РџСЂРёРјРµРЅСЏРµРј РїР°СЂР°РјРµС‚СЂС‹ Рє РІСЃРїРѕРјРѕРіР°С‚РµР»СЊРЅС‹Рј РґРµР»РµРЅРёСЏРј
                   direction = 'in',    #  Р РёСЃСѓРµРј РґРµР»РµРЅРёСЏ РІРЅСѓС‚СЂРё Рё СЃРЅР°СЂСѓР¶Рё РіСЂР°С„РёРєР°
                   #length = 10,    #  Р”Р»РёРЅРЅР° РґРµР»РµРЅРёР№
                   #width = 2,     #  РЁРёСЂРёРЅР° РґРµР»РµРЅРёР№
                   #color = 'm',    #  Р¦РІРµС‚ РґРµР»РµРЅРёР№
                   #pad = 10,    #  Р Р°СЃСЃС‚РѕСЏРЅРёРµ РјРµР¶РґСѓ С‡РµСЂС‚РѕС‡РєРѕР№ Рё РµРµ РїРѕРґРїРёСЃСЊСЋ
                   labelsize = 14,    #  Р Р°Р·РјРµСЂ РїРѕРґРїРёСЃРё
                   #labelcolor = 'r',    #  Р¦РІРµС‚ РїРѕРґРїРёСЃРё
                   bottom = True,    #  Р РёСЃСѓРµРј РјРµС‚РєРё СЃРЅРёР·Сѓ
                   top = True,    #   СЃРІРµСЂС…Сѓ
                   left = True,    #  СЃР»РµРІР°
                   right = True)    #  Рё СЃРїСЂР°РІР°

fig.savefig('Fermi_Dirac.png', dpi=300, transparent=False, bbox_inches = 'tight')
