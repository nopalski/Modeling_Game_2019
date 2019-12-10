# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 20:05:33 2019

@author: Nic
"""
import random
import sys
import pandas as pd
from modeling_game import set_mRNA_data, set_prot_data, make_model_combos, make_model, min_full

set_mRNA_data('mRNA.csv')
set_prot_data('prot.csv')

#guesses = {
#    'mRNA1': ['iA4A'],
#    'mRNA2': [],
#    'mRNA3': [],
#    'mRNA4': [],
#    'mRNA5': [],
#    'mRNA6': ['1R7A'],
#    'mRNA7': [],
#    'mRNA8': []
#}

guesses = {
    'mRNA1': ['iA4A'],
    'mRNA2': ['2R4A', '4A5R', '4A', '3R4A', '4A6R'],
    'mRNA3': ['6A', '2R6A'],
    'mRNA4': ['2R5R', '2R3R', '2R6R'],
    'mRNA5': ['2R6A', '6A', '1R6A', '4R6A'],
    'mRNA6': ['1R7A'],
    'mRNA7': ['8R', '2R8R', '4R8R'],
    'mRNA8': ['1R7R']
}

#guesses = {
#    'mRNA1': ['iA4A'],
#    'mRNA2': ['2R4A', '4A5R', '4A', '3R4A'],
#    'mRNA3': ['6A', '2R6A', '1R6A'],
#    'mRNA4': ['2R5R', '2R6R'],
#    'mRNA5': ['2R6A', '6A', '1R6A', '4R6A', '5R6A'],
#    'mRNA6': ['1R7A'],
#    'mRNA7': ['8R', '2R8R', '4R8R'],
#    'mRNA8': ['1R7R']
#}

splits = [[0, 60], [60, 120], [120, 180], [180, 240], [240, 300], [300, 360]]
split = splits[int(sys.argv[1])]
#split = [0, 3]
results = []
all_combos = make_model_combos(guesses)
all_combos = all_combos[split[0]:split[1]]
#all_combos = random.sample(all_combos, 5)
count = split[0]
print('started fitting ' + str(len(all_combos)) + ' models...')
for combo in all_combos:
    model = make_model(combo)
    params, chisqr = min_full(model)
    results.append([str(combo), chisqr, str(dict(params))])
    print('fitted ' + str(count) + ' of ' + str(len(all_combos)))
    count = count + 1
df = pd.DataFrame(results, columns=['regulators', 'chisqr', 'params'])
df.to_csv('modeling_' + sys.argv[1] + '.csv', index=False)