# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 22:26:52 2019

@author: Nic
"""

import ast
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from modeling_game import get_dims, make_model, get_sim_time, get_num_points

#results = pd.read_csv('csvs/modeling_full_1.csv')
results = pd.read_csv('csvs/modeling_full_1.csv')

redchi = []

for idxx in range(len(results)):
    redchi.append(results.iloc[idxx]['chisqr'] / len(ast.literal_eval(results.iloc[idxx]['params'])))

results['r_chisqr'] = redchi



#mRNA_data = pd.read_csv('mRNA.csv')
mRNA_data = pd.read_csv('mRNA.csv')
prot_data = pd.read_csv('prot.csv')

top_plots = 45

results = results.sort_values('r_chisqr')
time = mRNA_data['time']
t = get_sim_time()
tp = get_num_points()

#tp = 1000
sum_sqrs = []
print('top models:')
for i in range(top_plots):
    sum_sqr = 0
    regs = ast.literal_eval(results.iloc[i][0])
    chisqr = results.iloc[i][1]
    params = ast.literal_eval(results.iloc[i][2])
    model = make_model(regs)
    rr = te.loada(model)
    for p in params:
        rr[p] = params[p]
    sim_data = rr.simulate(0, t, tp)
    rows, cols = get_dims(top_plots)
    for name_type in ['mRNA', 'P']:
        plt.figure(name_type)
        plt.subplot(rows, cols, i + 1)
        for j in range(1, 9):
            name = name_type + str(j)
            name_te = '[' + name + ']'
            if name_type == 'P':
                1
                plt.plot(time, prot_data[name], color='C' + str(j), linewidth=1, linestyle='', marker='.', markersize=2)
                # sum_sqr = sum_sqr + ((prot_data[name] - sim_data[name_te]) ** 2).sum()
            else:
                plt.plot(time, mRNA_data[name], color='C' + str(j), linewidth=1, linestyle='', marker='.', markersize=2)
                #sum_sqr = sum_sqr + ((mRNA_data[name] - sim_data[name_te]) ** 2).sum()
            plt.plot(time, sim_data[name_te], color='C' + str(j))
            #plt.legend()
        #plt.title(str(regs) + ', chisqr: ' + str(round(chisqr, 2)), fontdict = {'fontsize' : 7})
        #plt.title('chisqr: ' + str(round(chisqr, 2)), fontdict = {'fontsize' : 12})
        #plt.xticks([])
        #plt.yticks([])
        #plt.title('model ' + str(i))
    if i < 10:
        s = str(i) + ' : '
    else:
        s = str(i) + ': '
    print(s + str(regs) + ', chisqr: ' + str(chisqr) + ', r_chisqr: ' + str(results.iloc[i][3]))
    #plt.legend()
    #sum_sqrs.append(sum_sqr)
plt.show()

    
#results['sumsqr'] = sum_sqrs
#results = results.sort_values('sumsqr')
#for i in range(9):
#    print(results.iloc[i][0], results.iloc[i][1], results.iloc[i][3])
    
#plt.figure()
#plt.plot(time, mRNA_data['mRNA1'], color='black', marker='.', linestyle=None)
#plt.plot(sim_data['time'], sim_data['[mRNA1]'], color='C4')
