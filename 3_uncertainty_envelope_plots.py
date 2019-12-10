import ast
import numpy as np
import tellurium as te
import matplotlib.pyplot as plt

rr = te.loada(open('mdl.txt', 'r').read())

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


time = np.linspace(0, 1200, 121)


Vm5_hist = np.histogram(Vm5_list, bins = 80)
Vm7_hist = np.histogram(Vm7_list, bins = 80)
K2_2_hist = np.histogram(K2_2_list, bins = 80)
K2_4_hist = np.histogram(K2_4_list, bins = 80)
K3_4_hist = np.histogram(K3_4_list, bins = 80)
H5_hist = np.histogram(H5_list, bins = 80)
K1_5_hist = np.histogram(K1_5_list, bins = 80)
H7_hist = np.histogram(H7_list, bins = 80)
K1_7_hist = np.histogram(K1_7_list, bins = 80)
K2_8_hist = np.histogram(K2_8_list, bins = 80)
K3_8_hist = np.histogram(K3_8_list, bins = 80)

import scipy.stats
Vm5_dist = scipy.stats.rv_histogram (Vm5_hist) 
Vm7_dist = scipy.stats.rv_histogram (Vm7_hist) 
K2_2_dist = scipy.stats.rv_histogram (K2_2_hist) 
K2_4_dist = scipy.stats.rv_histogram (K2_4_hist) 
K3_4_dist = scipy.stats.rv_histogram (K3_4_hist) 
H5_dist = scipy.stats.rv_histogram (H5_hist) 
K1_5_dist = scipy.stats.rv_histogram (K1_5_hist) 
H7_dist = scipy.stats.rv_histogram (H7_hist) 
K1_7_dist = scipy.stats.rv_histogram (K1_7_hist) 
K2_8_dist = scipy.stats.rv_histogram (K2_8_hist) 
K3_8_dist = scipy.stats.rv_histogram (K3_8_hist) 


mRNA1_data = [[] for i in range(121)]
P4_data = [[] for i in range(121)]
mRNA2_data = [[] for i in range(121)]
P2_data = [[] for i in range(121)]
mRNA3_data = [[] for i in range(121)]
P6_data = [[] for i in range(121)]
mRNA4_data = [[] for i in range(121)]
P5_data = [[] for i in range(121)]
mRNA5_data = [[] for i in range(121)]
mRNA6_data = [[] for i in range(121)]
P7_data = [[] for i in range(121)]
P1_data = [[] for i in range(121)]
mRNA7_data = [[] for i in range(121)]
P8_data = [[] for i in range(121)]
mRNA8_data = [[] for i in range(121)]
P3_data = [[] for i in range(121)]


for i in range(200):
    rr.reset()
    rr.Vm5 = Vm5_dist.rvs() 
    rr.Vm7 = Vm7_dist.rvs()
    rr.K2_2 = K2_2_dist.rvs()
    rr.K2_4 = K2_4_dist.rvs()
    rr.K3_4 = K3_4_dist.rvs()
    rr.H5 = H5_dist.rvs()
    rr.K1_5 = K1_5_dist.rvs()
    rr.H7 = H7_dist.rvs()
    rr.K1_7 = K1_7_dist.rvs()
    rr.K2_8 = K2_8_dist.rvs()
    rr.K3_8 = K3_8_dist.rvs()
    data = rr.simulate(0, 1200, 121)
    for j in range(121):
        mRNA1_data[j].append(data[j,1])
        P4_data[j].append(data[j,2])
        mRNA2_data[j].append(data[j,3])
        P2_data[j].append(data[j,4])
        mRNA3_data[j].append(data[j,5])
        P6_data[j].append(data[j,6])
        mRNA4_data[j].append(data[j,7])
        P5_data[j].append(data[j,8])
        mRNA5_data[j].append(data[j,9])
        mRNA6_data[j].append(data[j,10])
        P7_data[j].append(data[j,11])
        P1_data[j].append(data[j,12])
        mRNA7_data[j].append(data[j,13])
        P8_data[j].append(data[j,14])
        mRNA8_data[j].append(data[j,15])
        P3_data[j].append(data[j,16])
        

mRNA1_percentiles = []
mRNA2_percentiles = []
mRNA3_percentiles = []
mRNA4_percentiles = []
mRNA5_percentiles = []
mRNA6_percentiles = []
mRNA7_percentiles = []
mRNA8_percentiles = []

P1_percentiles = []
P2_percentiles = []
P3_percentiles = []
P4_percentiles = []
P5_percentiles = []
P6_percentiles = []
P7_percentiles = []
P8_percentiles = []

mRNA1_mean = []
mRNA2_mean = []
mRNA3_mean = []
mRNA4_mean = []
mRNA5_mean = []
mRNA6_mean = []
mRNA7_mean = []
mRNA8_mean = []

P1_mean = []
P2_mean = []
P3_mean = []
P4_mean = []
P5_mean = []
P6_mean = []
P7_mean = []
P8_mean = []


# Generate means and percentiles of time course data
for i in range(121):
    
    mRNA1_percentiles.append(np.percentile(mRNA1_data[i], [2.5, 97.5]))
    mRNA2_percentiles.append(np.percentile(mRNA2_data[i], [2.5, 97.5]))
    mRNA3_percentiles.append(np.percentile(mRNA3_data[i], [2.5, 97.5]))
    mRNA4_percentiles.append(np.percentile(mRNA4_data[i], [2.5, 97.5]))
    mRNA5_percentiles.append(np.percentile(mRNA5_data[i], [2.5, 97.5]))
    mRNA6_percentiles.append(np.percentile(mRNA6_data[i], [2.5, 97.5]))
    mRNA7_percentiles.append(np.percentile(mRNA7_data[i], [2.5, 97.5]))
    mRNA8_percentiles.append(np.percentile(mRNA8_data[i], [2.5, 97.5]))
    
    P1_percentiles.append(np.percentile(P1_data[i], [2.5, 97.5]))
    P2_percentiles.append(np.percentile(P2_data[i], [2.5, 97.5]))
    P3_percentiles.append(np.percentile(P3_data[i], [2.5, 97.5]))
    P4_percentiles.append(np.percentile(P4_data[i], [2.5, 97.5]))
    P5_percentiles.append(np.percentile(P5_data[i], [2.5, 97.5]))
    P6_percentiles.append(np.percentile(P6_data[i], [2.5, 97.5]))
    P7_percentiles.append(np.percentile(P7_data[i], [2.5, 97.5]))
    P8_percentiles.append(np.percentile(P8_data[i], [2.5, 97.5]))
    
    mRNA1_mean.append(np.mean(mRNA1_data[i]))
    mRNA2_mean.append(np.mean(mRNA2_data[i]))
    mRNA3_mean.append(np.mean(mRNA3_data[i]))
    mRNA4_mean.append(np.mean(mRNA4_data[i]))
    mRNA5_mean.append(np.mean(mRNA5_data[i]))
    mRNA6_mean.append(np.mean(mRNA6_data[i]))
    mRNA7_mean.append(np.mean(mRNA7_data[i]))
    mRNA8_mean.append(np.mean(mRNA8_data[i]))
    
    P1_mean.append(np.mean(P1_data[i]))
    P2_mean.append(np.mean(P2_data[i]))
    P3_mean.append(np.mean(P3_data[i]))
    P4_mean.append(np.mean(P4_data[i]))
    P5_mean.append(np.mean(P5_data[i]))
    P6_mean.append(np.mean(P6_data[i]))
    P7_mean.append(np.mean(P7_data[i]))
    P8_mean.append(np.mean(P8_data[i]))


mRNA1_percentiles = np.array(mRNA1_percentiles)
mRNA2_percentiles = np.array(mRNA2_percentiles)
mRNA3_percentiles = np.array(mRNA3_percentiles)
mRNA4_percentiles = np.array(mRNA4_percentiles)
mRNA5_percentiles = np.array(mRNA5_percentiles)
mRNA6_percentiles = np.array(mRNA6_percentiles)
mRNA7_percentiles = np.array(mRNA7_percentiles)
mRNA8_percentiles = np.array(mRNA8_percentiles)

P1_percentiles = np.array(P1_percentiles)
P2_percentiles = np.array(P2_percentiles)
P3_percentiles = np.array(P3_percentiles)
P4_percentiles = np.array(P4_percentiles)
P5_percentiles = np.array(P5_percentiles)
P6_percentiles = np.array(P6_percentiles)
P7_percentiles = np.array(P7_percentiles)
P8_percentiles = np.array(P8_percentiles)


# # Plotting mRNA data uncertainty envelope

plt.figure()



plt.subplot(4, 4, 1)
plt.plot(time, mRNA1_mean, label='mRNA1', color='C1')
plt.fill_between(time, mRNA1_percentiles[:,0], mRNA1_percentiles[:,1], color='C1', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA1 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 2)
plt.plot(time, mRNA2_mean, label='mRNA2', color='C2')
plt.fill_between(time, mRNA2_percentiles[:,0], mRNA2_percentiles[:,1], color='C2', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA2 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 3)
plt.plot(time, mRNA3_mean, label='mRNA3', color='C3')
plt.fill_between(time, mRNA3_percentiles[:,0], mRNA3_percentiles[:,1], color='C3', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA3 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 4)
plt.plot(time, mRNA4_mean, label='mRNA4', color='C4')
plt.fill_between(time, mRNA4_percentiles[:,0], mRNA4_percentiles[:,1], color='C4', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA4 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 5)
plt.plot(time, mRNA5_mean, label='mRNA5', color='C5')
plt.fill_between(time, mRNA5_percentiles[:,0], mRNA5_percentiles[:,1], color='C5', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA5 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 6)
plt.plot(time, mRNA6_mean, label='mRNA6', color='C6')
plt.fill_between(time, mRNA6_percentiles[:,0], mRNA6_percentiles[:,1], color='C6', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA6 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 7)
plt.plot(time, mRNA7_mean, label='mRNA7', color='C7')
plt.fill_between(time, mRNA7_percentiles[:,0], mRNA7_percentiles[:,1], color='C7', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA7 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 8)
plt.plot(time, mRNA8_mean, label='mRNA8', color='C8')
plt.fill_between(time, mRNA8_percentiles[:,0], mRNA8_percentiles[:,1], color='C8', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of mRNA8 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

# Plotting protein data uncertainty envelope

plt.subplot(4, 4, 9)
plt.plot(time, P1_mean, label='P1', color='C1')
plt.fill_between(time, P1_percentiles[:,0], P1_percentiles[:,1], color='C1', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P1 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 10)
plt.plot(time, P2_mean, label='P2', color='C2')
plt.fill_between(time, P2_percentiles[:,0], P2_percentiles[:,1], color='C2', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P2 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 11)
plt.plot(time, P3_mean, label='P3', color='C3')
plt.fill_between(time, P3_percentiles[:,0], P3_percentiles[:,1], color='C3', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P3 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 12)
plt.plot(time, P4_mean, label='P4', color='C4')
plt.fill_between(time, P4_percentiles[:,0], P4_percentiles[:,1], color='C4', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P4 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 13)
plt.plot(time, P5_mean, label='P5', color='C5')
plt.fill_between(time, P5_percentiles[:,0], P5_percentiles[:,1], color='C5', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P5 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 14)
plt.plot(time, P6_mean, label='P6', color='C6')
plt.fill_between(time, P6_percentiles[:,0], P6_percentiles[:,1], color='C6', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P6 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 15)
plt.plot(time, P7_mean, label='P7', color='C7')
plt.fill_between(time, P7_percentiles[:,0], P7_percentiles[:,1], color='C7', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P7 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(4, 4, 16)
plt.plot(time, P8_mean, label='P8', color='C8')
plt.fill_between(time, P8_percentiles[:,0], P8_percentiles[:,1], color='C8', alpha=0.2)
#plt.title('Envelope graph showing mean, upper and lower bounds of P8 values')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.show()








































