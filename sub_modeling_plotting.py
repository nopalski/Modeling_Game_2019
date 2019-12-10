# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 17:39:26 2019

@author: Nic
"""

import ast
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from modeling_game import get_dims, generate_sub_model, get_sim_time, get_num_points, set_mRNA_data, set_prot_data

set_mRNA_data('mRNA.csv')
set_prot_data('prot.csv')

#results = pd.read_csv('sub_modeling_1.csv')
#mRNA_data = pd.read_csv('mRNA.csv')

results = pd.read_csv('csvs/sub_modeling_1.csv')
mRNA_data = pd.read_csv('mRNA.csv')

mRNA_num = 8
top_plots = 9

mRNA_name = 'mRNA' + str(mRNA_num)
mRNA_result = results[results['mRNA_num'] == mRNA_num]
mRNA_result = mRNA_result.sort_values('chisqr')
time = mRNA_data['time']
t = get_sim_time()
tp = get_num_points()
plt.figure()
print('top regs for ' + mRNA_name + ':')
for i in range(top_plots):
    reg = mRNA_result.iloc[i][1]
    chisqr = mRNA_result.iloc[i][2] 
    params = ast.literal_eval(mRNA_result.iloc[i][3])
    sub_model = generate_sub_model(mRNA_num, reg)
    rr = te.loada(sub_model)
    for p in params:
        rr[p] = params[p]
    sim_data = rr.simulate(0, t, tp)
    rows, cols = get_dims(top_plots)
    plt.subplot(rows, cols, i + 1)
    plt.plot(time, mRNA_data[mRNA_name], color='C' + str(mRNA_num), label='exp')
    plt.plot(time, sim_data['[mRNA]'], color='black', label='sim')
    plt.title('reg: ' + reg + ', chisqr: ' + str(round(chisqr, 2)))
    plt.legend()
    print(reg + ', chisqr: ' + str(chisqr))
plt.suptitle('mRNA' + str(mRNA_num))
plt.show()













def p(df):
    print(df[['mRNA_num', 'regulators', 'chisqr']])
