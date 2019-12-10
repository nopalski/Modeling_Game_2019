# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:45:01 2019

@author: Nic
"""
import random
import pandas as pd
from modeling_game import set_mRNA_data, set_prot_data, generate_regs, generate_mRNA_regs, generate_sub_model, min_sub

set_mRNA_data('mRNA.csv')
set_prot_data('prot.csv')

results = []

print('starting...')
mRNA_nums = [2, 3, 4, 5, 7, 8]
#mRNA_nums = [2]
for num in mRNA_nums:
    all_regs = generate_mRNA_regs(num)
    #all_regs = generate_regs(None)
    #all_regs = random.sample(all_regs, 10)
    #all_regs = [all_regs[0]]
    for reg in all_regs:
        sub_model = generate_sub_model(num, reg)
        #print(sub_model)
        params, chisqr = min_sub(sub_model, num)
        results.append([num, reg, chisqr, str(dict(params))])
        print('fitted mRNA' + str(num) + ' for ' + reg)

df = pd.DataFrame(results, columns=['mRNA_num', 'regulators', 'chisqr', 'params'])
df.to_csv('sub_modeling.csv', index=False)
