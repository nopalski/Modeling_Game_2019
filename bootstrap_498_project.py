import tellurium as te
import numpy as np
import pandas as pd
import lmfit as lmfit
import matplotlib.pyplot as plt
from random import sample
import sys



exp = pd.read_csv('all.csv')
alll = ['mRNA1', 'mRNA2', 'mRNA3', 'mRNA4', 'mRNA5', 'mRNA6', 'mRNA7', 'mRNA8', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8']
f = open('mdl.txt', 'r')
model = f.read()
rr = te.loada(model)
rr.timeCourseSelections = alll
master = rr.simulate(0, 1200, 121)
residuals = (master - exp).to_numpy()
time = np.linspace(0, 1200, 121)


def make_resid_fnc(synth):
    def residuals_fnc(p):
        rr.reset()
        params = p.valuesdict()
        for param in params:
            rr[param] = params[param]
        data = rr.simulate(0, 1200, 121)
        res = (data - synth).ravel()
        return res
    return residuals_fnc

parameters = lmfit.Parameters()
parameters.add('Vm5', value=0.9167186340378028, min=0.5, max=2.0)
parameters.add('Vm7', value=1.4775023746542497, min=0.5, max=2.0)
parameters.add('K2_2', value=0.01, min=0.01, max=0.03)
parameters.add('K2_4', value=0.016845150801881946, min=0.01, max=0.03)
parameters.add('K3_4', value=0.01991907317992471, min=0.01, max=0.03)
parameters.add('H5', value=3.703392608035244, min=2.0, max=8.0)
parameters.add('K1_5', value=0.023025286249435378, min=0.01, max=0.03)
parameters.add('H7', value=3.5230017323544534, min=2.0, max=8.0)
parameters.add('K1_7', value=0.011344375364800022, min=0.01, max=0.03)
parameters.add('K2_8', value=0.012470037814671847, min=0.01, max=0.03)
parameters.add('K3_8', value=0.020519681925885457, min=0.01, max=0.03)


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


for i in range(500):
    synth_data = np.copy(master)
    for j in range(synth_data.shape[0]):
        for k in range(synth_data.shape[1]):
            synth_data[j,k] = max(0, synth_data[j,k] + np.random.choice(residuals[:,k]))
    fitter = lmfit.Minimizer(make_resid_fnc(synth_data), parameters)
    fit = fitter.minimize(method='differential_evolution')
    fit = fitter.minimize(method='leastsq', params=fit.params)
    Vm5_list.append(fit.params.valuesdict()['Vm5'])
    Vm7_list.append(fit.params.valuesdict()['Vm7'])
    K2_2_list.append(fit.params.valuesdict()['K2_2'])
    K2_4_list.append(fit.params.valuesdict()['K2_4'])
    K3_4_list.append(fit.params.valuesdict()['K3_4'])
    H5_list.append(fit.params.valuesdict()['H5'])
    K1_5_list.append(fit.params.valuesdict()['K1_5'])
    H7_list.append(fit.params.valuesdict()['H7'])
    K1_7_list.append(fit.params.valuesdict()['K1_7'])
    K2_8_list.append(fit.params.valuesdict()['K2_8'])
    K3_8_list.append(fit.params.valuesdict()['K3_8'])
    print('completed', i, 'iterations of bootstrapping')


# Compute 5th and 95th percentiles of fitted parameter values for Vm5

percent_Vm5_high = np.percentile(Vm5_list, 97.5)
percent_Vm5_low = np.percentile(Vm5_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for Vm7

percent_Vm7_high = np.percentile(Vm7_list, 97.5)
percent_Vm7_low = np.percentile(Vm7_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K2_2

percent_K2_2_high = np.percentile(K2_2_list, 97.5)
percent_K2_2_low = np.percentile(K2_2_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K2_4

percent_K2_4_high = np.percentile(K2_4_list, 97.5)
percent_K2_4_low = np.percentile(K2_4_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K3_4

percent_K3_4_high = np.percentile(K3_4_list, 97.5)
percent_K3_4_low = np.percentile(K3_4_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for H5

percent_H5_high = np.percentile(H5_list, 97.5)
percent_H5_low = np.percentile(H5_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K1_5

percent_K1_5_high = np.percentile(K1_5_list, 97.5)
percent_K1_5_low = np.percentile(K1_5_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for H7

percent_H7_high = np.percentile(H7_list, 97.5)
percent_H7_low = np.percentile(H7_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K1_7

percent_K1_7_high = np.percentile(K1_7_list, 97.5)
percent_K1_7_low = np.percentile(K1_7_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K2_8

percent_K2_8_high = np.percentile(K2_8_list, 97.5)
percent_K2_8_low = np.percentile(K2_8_list, 2.5)

# Compute 5th and 95th percentiles of fitted parameter values for K3_8

percent_K3_8_high = np.percentile(K3_8_list, 97.5)
percent_K3_8_low = np.percentile(K3_8_list, 2.5)


 Print percentile for K1 and K2

print("For parameter Vm5, 95% of fitted estimates fall between " + str(percent_Vm5_low) + " and " + str(percent_Vm5_high) + ".")
print("For parameter Vm7, 95% of fitted estimates fall between " + str(percent_Vm7_low) + " and " + str(percent_Vm7_high) + ".")
print("For parameter K2_2, 95% of fitted estimates fall between " + str(percent_K2_2_low) + " and " + str(percent_K2_2_high) + ".")
print("For parameter K2_4, 95% of fitted estimates fall between " + str(percent_K2_4_low) + " and " + str(percent_K2_4_high) + ".")
print("For parameter K3_4, 95% of fitted estimates fall between " + str(percent_K3_4_low) + " and " + str(percent_K3_4_high) + ".")
print("For parameter H5, 95% of fitted estimates fall between " + str(percent_H5_low) + " and " + str(percent_H5_high) + ".")
print("For parameter K1_5, 95% of fitted estimates fall between " + str(percent_K1_5_low) + " and " + str(percent_K1_5_high) + ".")
print("For parameter H7, 95% of fitted estimates fall between " + str(percent_H7_low) + " and " + str(percent_H7_high) + ".")
print("For parameter K1_7, 95% of fitted estimates fall between " + str(percent_K1_7_low) + " and " + str(percent_K1_7_high) + ".")
print("For parameter K2_8, 95% of fitted estimates fall between " + str(percent_K2_8_low) + " and " + str(percent_K2_8_high) + ".")
print("For parameter K3_8, 95% of fitted estimates fall between " + str(percent_K3_8_low) + " and " + str(percent_K3_8_high) + ".")

f2 = open("plist_" + sys.argv[1] + ".txt","w+")
f2.write(str(Vm5_list))
f2.write('\n')
f2.write(str(Vm7_list))
f2.write('\n')
f2.write(str(K2_2_list))
f2.write('\n')
f2.write(str(K2_4_list))
f2.write('\n')
f2.write(str(K3_4_list))
f2.write('\n')
f2.write(str(H5_list))
f2.write('\n')
f2.write(str(K1_5_list))
f2.write('\n')
f2.write(str(H7_list))
f2.write('\n')
f2.write(str(K1_7_list))
f2.write('\n')
f2.write(str(K2_8_list))
f2.write('\n')
f2.write(str(K3_8_list))
f2.close()










































