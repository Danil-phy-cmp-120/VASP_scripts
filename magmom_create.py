#!/usr/bin/python
# -*- coding: utf-8 -*-

i = 1
a = ''
magmom = ''
while a != 'end':
    a = input('Enter magmom for {} atom (3 coordinates or end):\n'.format(i))

    if a != 'end':
        n = int(input('Enter the number of {} atoms:\n'.format(i)))
        magmom = magmom + (a + ' ')*n

    i = i + 1

print('MAGMOM = {}'.format(magmom))



