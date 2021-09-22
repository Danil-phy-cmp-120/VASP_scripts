# -*- coding: utf-8 -*-

from chempy import balance_stoichiometry
import numpy as np
import itertools 

#################################################################################################################################################################################################

num = 3 # Число соединений на которые распадается исходный сплав

all_reac = {'Ni7Co1Mn7In1': -113.98884} # Все энергии должны быть приведены на число атомов указанном ключах словарей
all_prod = {'Co1': -7.1083, 'Co3Ni1': -27.2084, 'In1': -2.7517, 'In3Co1': -15.416, 'In7Ni3': -38.02, 'In1Ni1': -8.8814, 'Mn1': -9.1617, 'Mn3Co1': -34.6588, 'Mn1Co1': -16.322, 'Mn1Ni3': -26.9488, 'Ni1': -5.7798, 'Mn2In1Co1': -27.3588, 'Mn1In1Ni2': -23.8236, 'Mn8Ni8': -122.427}

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

data = []
for i in itertools.combinations(all_prod.keys(), num):
    try:
        reac, prod = balance_stoichiometry(all_reac.keys(), i)
        if (sum(1 for number in  prod.values() if number < 0)) == 0:
        
            E_formation = all_reac[list(all_reac.keys())[0]]
            reactions = '{} ->'.format(list(all_reac.keys())[0])

            for j in prod.keys():
                E_formation -= prod[j] * all_prod[j] 
                reactions += ' {} +'.format(j)
                
            reactions[:len(reactions)-2]    
            E_formation = E_formation/find_sum_str(list(all_reac.keys())[0]) 
            
            data += [[reactions[:len(reactions)-2], E_formation]]
                 
    except:
        continue
 
 
with open('reactions.dat', 'w') as filehandle:  
    for d in data:
        filehandle.write('{}\t{}\n'.format(d[0], d[1]))
