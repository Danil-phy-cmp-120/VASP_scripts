# -*- coding: utf-8 -*-

from chempy import balance_stoichiometry
import numpy as np
import itertools 

#################################################################################################################################################################################################

num = 3 # Число соединений на которые распадается исходный сплав

all_reac = {'Mn2Sc1Si1': -61.9275297725} # Все энергии должны быть приведены на число атомов указанном ключах словарей
all_prod = {'Mn2Sc1Si2': -74.20357883, 'Mn3Si1': -66.33381478, 'Mn4Si4': -116.95593121, 'Mn8Sc4': -203.82809295, 'Mn29': -528.15173, 'Sc2': -28.4978581, 'Sc2Si2': -51.90588478, 'Sc10Si6': -215.2864648, 'Si2': -20.0152677, 'Mn16Si28': -589.01413457, 'Mn6Sc4Si2': -190.6439151, 'Mn12Sc12Si24': -661.56152365} # Даже если атом один после него необходимо ставить "1"

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
            data += [[reactions[:len(reactions)-2], E_formation]]
                 
    except:
        continue
 
 
with open('reactions.dat', 'w') as filehandle:  
    for d in data:
        filehandle.write('{}\t{}\n'.format(d[0], d[1]))
