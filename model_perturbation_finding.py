# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 23:32:05 2019

@author: Nic
"""

import ast
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
import lmfit
from modeling_game import get_dims, make_model
from itertools import combinations, chain
import numpy as np

data = pd.read_csv('csvs/modeling_full_1.csv')
top_models = [1, 6]
tps = 4800








def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def make_params(typ):
    parameters = lmfit.Parameters()
    if typ == 'up_reg':
        parameters.add('times', min=1, max=1.75)
    elif typ == 'down_reg':
        parameters.add('times', min=0.25, max=1)
    else:
        raise ValueError('invalid typ: ', typ)
    return parameters

def make_resid_fnc(poss):
    rr_models, orig_params = make_all_rrs(poss)
    last = [None]
    def residuals(pp):
        f_params = pp.valuesdict()
        all_sim_data = []
        for i in top_models:
            all_sim_data.append(sim_model(rr_models[i], orig_params[i], f_params))
        sum_sqr = calc_sum_sqr(all_sim_data)
        last[0] = sum_sqr
        return sum_sqr
    return last, residuals

def reset_rr(rr, orig_params):
    rr.resetToOrigin()
    for p in orig_params:
        rr[p] = orig_params[p]

def make_all_rrs(models):
    all_rrs = []
    all_orig_params = []
    for m in models:
        regs = ast.literal_eval(data.iloc[m][0])
        model_params = ast.literal_eval(data.iloc[m][2])
        model = make_model(regs)
        rr = te.loada(model)
        rr.timeCourseSelections = ['mRNA1', 'mRNA2', 'mRNA3', 'mRNA4', 'mRNA5', 'mRNA6', 'mRNA7', 'mRNA8']
        for p in model_params:
            rr[p] = model_params[p]
        all_rrs.append(rr)
        all_orig_params.append(model_params)
    return (all_rrs, all_orig_params)

def sim_model(rr, o_params, f_params):
    rr.resetToOrigin()
    for o in o_params:
        rr[o] = o_params[o]
    for f in f_params:
        rr[f] = rr[f]*f_params[f]
    return rr.simulate(0, 1200, tps)

def calc_sum_sqr(all_sim_data):
    sum_sqr = ((all_sim_data[0] - all_sim_data[1]) ** 2).sum()
    return 0 - sum_sqr





data = data.sort_values('chisqr')
gene_list = ['1', '2', '3', '4', '5', '6', '7', '8']
vals = [round(n, 2) for n in np.arange(0.25, 1.750001, 0.05) if round(n, 2) != 1]
vals.append(0)

gene_combos = list(powerset(gene_list))[1:-1]
model_combos = list(combinations(top_models, 2))
gene_combos = [gc for gc in gene_combos if len(gc) < 4]

for models in model_combos:
    all_rrs, all_orig_params = make_all_rrs(models)
    results = []
    for genes in gene_combos:
        for v in vals:
            all_sim_data = []
            for i, rr in enumerate(all_rrs):
                reset_rr(rr, all_orig_params[i])
                for g in genes:
                    if v == 0:
                        rr['L' + g] =  0
                    rr['Vm' + g] = rr['Vm' + g]*v
                all_sim_data.append(rr.simulate(0, 1200, tps))
            for col in all_sim_data[0].colnames:
                sum_sqr = ((all_sim_data[0][col] - all_sim_data[1][col]) ** 2).sum()
                results.append([str(list(genes)), sum_sqr, v, col])
            #sum_sqr = calc_sum_sqr(all_sim_data)
                
    df = pd.DataFrame(results, columns=['genes', 'sumsqr', 'value', 'species'])
    df = df.sort_values('sumsqr', ascending=False)
    df.to_csv('perturbations_' + str(models[0]) + '_' + str(models[1]) + '.csv', index=False)
    print('finished models:', models)

    
    

