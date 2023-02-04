# -*- coding: utf-8 -*-

from chempy import balance_stoichiometry
import numpy as np
import itertools 

#################################################################################################################################################################################################

num = 3 # Число соединений на которые распадается исходный сплав

all_reac = {'Ni7Co1Mn7In1': -113.98884} # Все энергии должны быть приведены на число атомов указанном ключах словарей
all_prod = {'Co1': -7.036345475, 'Co3Ni1': -26.69953137, 'Ga1': -2.90738149125, 'Ga3Co1': -16.7982716125, 'Ga3Ni2': -21.72045611, 'Ga3Ni5': -39.108281075, 'Ga7Ni3': -40.050000875, 'Ga9Ni13': -105.990912015, 'Ga1Co1': -10.47789089, 'Ga1Ni3': -20.52036508, 'Mn1': -8.9993800177586206896551724137931, 'Mn3Co1': -34.10424484, 'Mn1Ga1': -12.16317968, 'Mn1Ga2Ni9': -66.95526976, 'Mn1Ga4': -21.47390636, 'Mn1Ga1Co2': -26.9066561825, 'Mn1Ga1Ni2': -24.04144694, 'Mn1Ni3': -25.91876489, 'Ni1': -5.4921941975, 'Mn1Co1': -16.07801210}
 # Даже если атом один после него необходимо ставить "1"

#################################################################################################################################################################################################


def find_sum_str(S):
    ans = 0
    last_number = 0
    for char in S:
        if '0' <= char <= '9':
            last_number = 10 * last_number + int(char)
        else:
            ans += last_number
            last_number = 0
    return(ans + last_number)
    
print('{} atoms in reac'.format(find_sum_str(list(all_reac.keys())[0])))

data = []
for i in itertools.combinations(all_prod.keys(), num):
    try:
        reac, prod = balance_stoichiometry(all_reac.keys(), i)
        if (sum(1 for number in  prod.values() if number < 0)) == 0:
        
            E_formation = all_reac[list(all_reac.keys())[0]]*reac[(list(reac.keys())[0])]
            reactions = '{} {} ->'.format(reac[list(all_reac.keys())[0]], list(all_reac.keys())[0])

            for j in prod.keys():
                E_formation -= prod[j] * all_prod[j] 
                reactions += ' {} {} +'.format(prod[j], j)    
            
            E_formation = E_formation/(reac[list(reac.keys())[0]] * find_sum_str(list(all_reac.keys())[0]))
            data += [[reactions[:len(reactions)-2], E_formation*1000]] # в мэВ
                 
    except:
        continue
 
 
with open('reactions.dat', 'w') as filehandle:  
    for d in data:
        filehandle.write('{}\t{}\n'.format(d[0], d[1]))
