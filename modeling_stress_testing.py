# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 16:00:22 2019

@author: Nic
"""

import ast
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from modeling_game import get_dims, make_model

results = pd.read_csv('csvs/modeling_full_1.csv')
prot = pd.read_csv('prot.csv')
KO1 = pd.read_csv('KO_1_5_P3.csv')
DN1 = pd.read_csv('DOWN_1_6_8_P3.csv')
UP1 = pd.read_csv('UP_4_5_6_P7.csv')
DN2 = pd.read_csv('DOWN_1_5_8_P6.csv')
UP2 = pd.read_csv('UP_7_P7.csv')

# model_num = 15
#plt.ioff()
stress_tests = [
#        {'species': 'P1', 'params':{'L1': 1.5}, 'ylim': (0, 6)}, 
        {'species': 'P1', 'params':{'K1_1': 1.5}, 'ylim': (0, 6)},
        {'species': 'P1', 'params':{'K1_1': 1.48}, 'ylim': (0, 6)},
#        {'species': 'P7', 'params':{'Vm2': 1.5}, 'ylim': (0, 9)},
#        {'species': 'P3', 'params':{'Vm1': 1.5}, 'ylim': (0, 9)},
#        {'species': 'P8', 'params':{'Vm5': 0, 'L5': 0}, 'ylim': (0, 6)},
#        {'species': 'P4', 'params':{'Vm2': 1.5, 'Vm8': 1.5}, 'ylim': (0, 3)}
        ]
#stress_tests = [
#                {'species': 'P3', 'params':{'L1': 0, 'Vm1': 0, 'L5': 0, 'Vm5': 0}, 'ylim': (0, 9)}, 
#                {'species': 'P3', 'params':{'Vm1': 0.39, 'Vm6': 0.39, 'Vm8': 0.39}, 'ylim': (0, 5)},
#                {'species': 'P7', 'params':{'Vm4': 1.7, 'Vm5': 1.7, 'Vm6': 1.7}, 'ylim': (0, 10)},
#                {'species': 'P6', 'params':{'Vm1': 0.55, 'Vm5': 0.55, 'Vm8': 0.55}, 'ylim': (0, 8)},
#                {'species': 'P7', 'params':{'Vm7': 1.25}, 'ylim': (0, 10)}
#                ]

#stress_tests = [
#                {'species': 'P1', 'params':{'Vm8': 0.5, 'Vm3': 1.5}, 'ylim': (0, 6)}, 
#                {'species': 'P3', 'params':{'Vm3': 0.25, 'Vm5': 1.25}, 'ylim': (0, 8)},
#                {'species': 'P4', 'params':{'Vm1': 0.34, 'Vm5': 1.25}, 'ylim': (0, 3)},
#                {'species': 'P4', 'params':{'L5': 0, 'Vm5': 0}, 'ylim': (0, 3)},
#                {'species': 'P5', 'params':{'K2_6': 1.33}, 'ylim': (0, 6)},
#                {'species': 'P8', 'params':{'Vm4': 2}, 'ylim': (0, 6)}
#                ]

results = results.sort_values('chisqr')
#plt.figure()
for idx, model_num in enumerate([0, 1, 2, 3]):
    sim_time = 1200
    regs = ast.literal_eval(results.iloc[model_num][0])
    params = ast.literal_eval(results.iloc[model_num][2])
    model = make_model(regs)
    rr = te.loada(model)
    for p in params:
        rr[p] = params[p]
    sim_data_wild = rr.simulate(0, sim_time, 5000)
    plt.figure('model ' + str(model_num), figsize=(10, 5), dpi=100)
    rows, cols = get_dims(len(stress_tests))
#    row = 1
#    cols = 3
    #rows, cols = get_dims(4)
    #plt.subplot(rows, cols, idx + 1)
    for i, test in enumerate(stress_tests):
        rr.resetToOrigin()
        for p in params:
            rr[p] = params[p]
        species = '[' + test['species'] + ']'
        title = species + ', '
        for p in test['params']:
            rr[p] = rr[p]*test['params'][p]
            title = title + p + '=' + str(test['params'][p]) + 'x, '
        sim_data_mut = rr.simulate(0, sim_time, 5000)
        plt.subplot(rows, cols, i + 1)
        #plt.ylim(test['ylim'])
        #plt.plot(sim_data_wild['time'], sim_data_wild[species], color='C' + species[-2], label='wt')
        plt.plot(sim_data_wild['time'], sim_data_wild[species], linewidth=2)
        # plt.plot(sim_data_mut['time'], sim_data_mut[species], color='black', label='mut')
        # plt.plot(sim_data_wild['time'], sim_data_wild[species], color='C0')
        #plt.plot(perf['time'], perf[species], color='C0', ls='--')
        # plt.plot(prot['time'], prot[test['species']], color='C0')
#        plt.plot(sim_data_mut['time'], sim_data_mut[species], color='C' + species[-2], linewidth=2, label=test['species'])
        plt.plot(sim_data_mut['time'], sim_data_mut[species], linewidth=2)
#        if i == 0:
#            True
#            plt.plot(KO1['time'], KO1['P3'], 'k.', label='')
#        elif i == 1:
#            True
#            plt.plot(DN1['time'], DN1['P3'], 'k.', label='')
#        elif i == 2:
#            True
#            plt.plot(UP1['time'], UP1['P7'], 'k.', label='')
#        elif i == 3:
#            True
#            plt.plot(DN2['time'], DN2['P6'], 'k.',label='')
#        elif i == 4:
#            True
#            plt.plot(UP2['time'], UP2['P7'], 'k.',label='')
        #plt.title(title[0:-2].replace('[', '').replace(']', ''))
        
        plt.legend()
        plt.title('K1_1 = ' + str(test['params']['K1_1']) + 'x')
    prefix = 'figs/stress_2/'
    #plt.suptitle('model ' + str(model_num))
    if model_num < 10:
        prefix = prefix + '00'
    elif model_num < 100:
        prefix = prefix + '0'
    plt.savefig(prefix + str(model_num) + '.png')
    plt.close()
    #plt.show()

    
    
#rr.L1 = rr.L1 * 1.5

#rr.K1_1 = rr.K1_1 * 1.5

#rr.Vm2 = rr.Vm2*1.5

#rr.Vm1 = rr.Vm1 * 1.5

#rr.Vm5 = 0;
#rr.L5 = 0;

#rr.Vm8 = rr.Vm8*1.5
#rr.Vm2 = rr.Vm2*1.5

#sim_data_mut = rr.simulate(0, sim_time, 5000)
#plt.figure()
#
#
#plt.title(species)
#plt.legend()

