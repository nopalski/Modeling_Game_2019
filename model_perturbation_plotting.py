# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 13:59:04 2019

@author: Nic
"""

import ast
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from modeling_game import get_dims, make_model

data = pd.read_csv('csvs/modeling_full_1.csv')
data = data.sort_values('chisqr')

results = pd.read_csv('csvs/perturbations_1.csv')
results = results.sort_values('sum_sqr')



model_nums = [0, 1, 2, 3]

perturb_num = 5

#for i in range(50):
#    print(results.iloc[i][0])
#    print(results.iloc[i][1])

rnd = False
places = 3
plus = 0

plt.figure()
sim_time = 1200
for i, num in enumerate(model_nums):
    regs = ast.literal_eval(data.iloc[num][0])
    chisqr = data.iloc[num][1]
    params = ast.literal_eval(data.iloc[num][2])
    model = make_model(regs)
    rr = te.loada(model)
    for p in params:
        rr[p] = params[p]
    sim_data = rr.simulate(0, sim_time, 5000)
    rr.reset()
    
    per = ast.literal_eval(results.iloc[perturb_num][0])
    sum_sqr = results.iloc[perturb_num][1]
    dict_str = results.iloc[perturb_num][2]
    if type(dict_str) == str:
        per_dict = ast.literal_eval(dict_str)
    else:
        per_dict = None
        
    if per_dict is not None:
        for p in per_dict:
            if rnd:
                per_dict[p] = round(per_dict[p], places)
            rr[p] = rr[p]*(per_dict[p] + plus)
            if i == 0:
                print('upregulated mRNA' + p[-1] + ' by ' + str(per_dict[p]) + 'x')
    for p in per:
        if p[0] == 'K':
            L = 'L' + p[-1]
            Vm = 'Vm' + p[-1]
            rr[L] = 0
            rr[Vm] = 0
            if i == 0:
                print('knocked out mRNA' + p[-1])
    
    mut_data = rr.simulate(0, sim_time, 5000)
    
    rows, cols = get_dims(len(model_nums))
    plt.subplot(rows, cols, i + 1)
    for col in sim_data.colnames:
        if col.count('mRNA') > 0:
            # plt.plot(sim_data['time'], sim_data[col], color='C' + col[-2], label=col[1:-1], linestyle='--')
            plt.plot(mut_data['time'], mut_data[col], color='C' + col[-2])
    plt.title('model ' + str(num))

# results.to_csv('derp.csv', index=False)