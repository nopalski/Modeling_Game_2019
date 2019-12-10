# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 18:04:58 2019

@author: Nic
"""

import matplotlib.pyplot as plt
import ast
import numpy as np
from itertools import combinations


Vm5_list = []
Vm7_list = []
K2_2_list = []
K2_4_list = []
K3_4_list = []
H5_list = []
K1_5_list = []
H7_list = []
K1_7_list = []
K2_8_list = []
K3_8_list = []


for i in range(1,7):
    f = open('plist_' + str(i) + '.txt', 'r')
    lines = f.read().split('\n')
    Vm5_list = Vm5_list + ast.literal_eval(lines[0])
    Vm7_list = Vm7_list + ast.literal_eval(lines[1])
    K2_2_list = K2_2_list + ast.literal_eval(lines[2])
    K2_4_list = K2_4_list + ast.literal_eval(lines[3])
    K3_4_list = K3_4_list + ast.literal_eval(lines[4])
    H5_list = H5_list + ast.literal_eval(lines[5])
    K1_5_list = K1_5_list + ast.literal_eval(lines[6])
    H7_list = H7_list + ast.literal_eval(lines[7])
    K1_7_list = K1_7_list + ast.literal_eval(lines[8])
    K2_8_list = K2_8_list + ast.literal_eval(lines[9])
    K3_8_list = K3_8_list + ast.literal_eval(lines[10])
    
print('Vm5', np.mean(Vm5_list), np.percentile(Vm5_list, 2.5), np.percentile(Vm5_list, 97.5))
print()
print('Vm7', np.mean(Vm7_list), np.percentile(Vm7_list, 2.5), np.percentile(Vm7_list, 97.5))
print()
print('K2_2', np.mean(K2_2_list), np.percentile(K2_2_list, 2.5), np.percentile(K2_2_list, 97.5))
print()
print('K2_4', np.mean(K2_4_list), np.percentile(K2_4_list, 2.5), np.percentile(K2_4_list, 97.5))
print()
print('K3_4', np.mean(K3_4_list), np.percentile(K3_4_list, 2.5), np.percentile(K3_4_list, 97.5))
print()
print('H5', np.mean(H5_list), np.percentile(H5_list, 2.5), np.percentile(H5_list, 97.5))
print()
print('K1_5', np.mean(K1_5_list), np.percentile(K1_5_list, 2.5), np.percentile(K1_5_list, 97.5))
print()
print('H7', np.mean(H7_list), np.percentile(H7_list, 2.5), np.percentile(H7_list, 97.5))
print()
print('K1_7', np.mean(K1_7_list), np.percentile(K1_7_list, 2.5), np.percentile(K1_7_list, 97.5))
print()
print('K2_8', np.mean(K2_8_list), np.percentile(K2_8_list, 2.5), np.percentile(K2_8_list, 97.5))
print()
print('K3_8', np.mean(K3_8_list), np.percentile(K3_8_list, 2.5), np.percentile(K3_8_list, 97.5))
print()

alll = [{'name': 'Vm5', 'list': Vm5_list}, 
        {'name': 'Vm7', 'list': Vm7_list}, 
        {'name': 'K2_2', 'list': K2_2_list}, 
        {'name': 'K2_4', 'list': K2_4_list}, 
        {'name': 'K3_4', 'list': K3_4_list}, 
        {'name': 'H5', 'list': H5_list}, 
        {'name': 'K1_5', 'list': K1_5_list}, 
        {'name': 'H7', 'list': H7_list}, 
        {'name': 'K1_7', 'list': K1_7_list},
        {'name': 'K2_8', 'list': K2_8_list}, 
        {'name': 'K3_8', 'list': K3_8_list}]

combos = list(combinations(alll, 2))

combos.reverse()

plt.figure()
cur = combos[0][0]['name']
row = 1
col = 0
for i, c in enumerate(combos):
    if c[0]['name'] != cur:
        row = row + 1
        cur = c[0]['name']
        col = 1
    else:
        col = col + 1
    # print(row, col, c[0]['name'], c[1]['name'])
    plt.subplot(10, 10, (row - 1)*10 + col)
    plt.plot(c[1]['list'], c[0]['list'], 'm.', markersize=1)
    if col == 1:
        plt.ylabel(c[0]['name'])
    else:
        plt.yticks([])
    if row == 10:
        plt.xlabel(c[1]['name'])
    else:
        plt.xticks([])
plt.show()
    

    
plt.figure()
plt.subplot(4, 3, 1)
plt.hist(Vm5_list, color = "g", bins = 80)
plt.title("Parameter Vm5 Distribution Histogram")
plt.xlabel("Vm5 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 2)
plt.hist(Vm7_list, color = "g", bins = 80)
plt.title("Parameter Vm7 Distribution Histogram")
plt.xlabel("Vm7 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 3)
plt.hist(K2_2_list, color = "g", bins = 80)
plt.title("Parameter K2_2 Distribution Histogram")
plt.xlabel("K2_2 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 4)
plt.hist(K2_4_list, color = "g", bins = 80)
plt.title("Parameter K2_4 Distribution Histogram")
plt.xlabel("K2_4 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 5)
plt.hist(K3_4_list, color = "g", bins = 80)
plt.title("Parameter K3_4 Distribution Histogram")
plt.xlabel("K3_4 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 6)
plt.hist(H5_list, color = "g", bins = 80)
plt.title("Parameter H5 Distribution Histogram")
plt.xlabel("H5 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 7)
plt.hist(K1_5_list, color = "g", bins = 80)
plt.title("Parameter K1_5 Distribution Histogram")
plt.xlabel("K1_5 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 8)
plt.hist(H7_list, color = "g", bins = 80)
plt.title("Parameter H7 Distribution Histogram")
plt.xlabel("H7 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 9)
plt.hist(K1_7_list, color = "g", bins = 80)
plt.title("Parameter K1_7 Distribution Histogram")
plt.xlabel("K1_7 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 10)
plt.hist(K2_8_list, color = "g", bins = 80)
plt.title("Parameter K2_8 Distribution Histogram")
plt.xlabel("K2_8 values")
plt.ylabel("Frequency")

plt.subplot(4, 3, 11)
plt.hist(K3_8_list, color = "g", bins = 80)
plt.title("Parameter K3_8 Distribution Histogram")
plt.xlabel("K3_8 values")
plt.ylabel("Frequency")
plt.show()

        

