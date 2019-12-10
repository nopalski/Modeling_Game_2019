# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:22:48 2019

@author: Nic
"""

import ast
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from modeling_game import get_dims, make_model

results = pd.read_csv('csvs/modeling_full_1.csv')
results = results.sort_values('chisqr')

mRNA_data = pd.read_csv('mRNA.csv')
KO = pd.read_csv('KO_1_5_P3.csv')
DN1 = pd.read_csv('DOWN_1_6_8_P3.csv')
UP1 = pd.read_csv('UP_4_5_6_P7.csv')
DN2 = pd.read_csv('DOWN_1_5_8_P6.csv')
UP2 = pd.read_csv('UP_7_P7.csv')

more_data = UP2
#model_nums = range(4)
model_nums = [1, 6]
species = ['P7']
#species = []
sim_time = 1200
plt.figure()
for i, num in enumerate(model_nums):
    regs = ast.literal_eval(results.iloc[num][0])
    chisqr = results.iloc[num][1]
    params = ast.literal_eval(results.iloc[num][2])
    model = make_model(regs)
    rr = te.loada(model)
    if len(species) > 0:
        rr.timeCourseSelections = ['time'] + species
    else:
        rr.timeCourseSelections = ['time', 'mRNA1', 'mRNA2', 'mRNA3', 'mRNA4', 'mRNA5', 'mRNA6', 'mRNA7', 'mRNA8']
        # rr.timeCourseSelections = ['time', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8']
    for p in params:
        rr[p] = params[p]
    sim_data = rr.simulate(0, sim_time, 5000)
    rr.reset()

#    rr['Vm1'] = 0
#    rr['L1'] = 0
#    rr['Vm5'] = 0
#    rr['L5'] = 0
#    rr['Vm8'] = 0
#    rr['L8'] = 0

    v = 1.25
    
    
#    rr['Vm1'] = rr['Vm1']*v
    
#    rr['Vm2'] = rr['Vm2']*v
#    
#    rr['Vm4'] = rr['Vm4']*v
#    
#    rr['Vm5'] = rr['Vm5']*v
##    
#    rr['Vm6'] = rr['Vm6']*v
    
    rr['Vm7'] = rr['Vm7']*v
    
#    rr['Vm8'] = rr['Vm8']*v
    
    mut_data = rr.simulate(0, sim_time, 5000)
    
    rows, cols = get_dims(len(model_nums))
    #plt.subplot(rows, cols, i + 1)
    for col in sim_data.colnames[1:]:
        if i == 0:
            plt.plot(mut_data['time'], mut_data[col], color='C' + col[-1], label=str(num), linewidth=1)
        else:
            plt.plot(mut_data['time'], mut_data[col], color='C' + col[-1], label=str(num), linewidth=1, linestyle='--')
        
        # plt.plot(sim_data['time'], sim_data[col], color='C' + col[-1], linewidth=0.5, linestyle='--')
        #plt.plot(mRNA_data['time'], mRNA_data[col], color='C' + col[-1], linewidth=0.75, linestyle='--')
        #plt.xlim([0, 200])
    # plt.legend()
#    if len(species) == 0:
#        for col in sim_data.colnames:
#            if col.count('mRNA') > 0:
#                # plt.plot(sim_data['time'], sim_data[col], color='C' + col[-2], label=col[1:-1], linestyle='--')
#                plt.plot(mut_data['time'], mut_data[col], color='C' + col[-2])
#    else:
#        for s in species:
#            name = '[' + s + ']'
#            # plt.plot(sim_data['time'], sim_data[name], color='C' + name[-2], label=name[1:-1], linestyle='--')
#            plt.plot(mut_data['time'], mut_data[name], color='C' + name[-2])
if more_data is not None:
    for col in more_data.columns[1:]:
        plt.plot(more_data['time'], more_data[col], 'k.', linewidth=1)
plt.title('model ' + str(model_nums[0]) + ' vs model ' + str(model_nums[1]))
plt.show()
plt.legend()
    
    
