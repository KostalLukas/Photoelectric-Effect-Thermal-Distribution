# -*- coding: utf-8 -*-
"""
Plotting IV Curves Analysis v1.1

Lukas Kostal, 25.10.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt


# specify dataset number to plot
ds_num = 1

# set true to also plot blackground curve
PBG = True

# specify available filters
filt = ['UV', 'V', 'B', 'G', 'Y', 'R']

# constant error in power supply output in V
psu_err_const = 0.005

# linear error in power supply as fraction of output voltage
psu_err_var = 5e-3

# append background to list of curves to plot
if PBG == True:
    filt.append('BG')
 
# arrays to hold arrays of averaged voltage for each dataset
V_arr = np.zeros(len(filt), dtype=object)
V_err = np.zeros(len(filt), dtype=object)
I_arr = np.zeros(len(filt), dtype=object)
I_err = np.zeros(len(filt), dtype=object)   
 
# average over repeated measurements for each filter
for i in range(0, len(filt)):
     
    # load the data voltage V in V photocurrent I in nA and convert to A
    V, I = np.loadtxt(f'Data/SerMes_{ds_num}_{filt[i]}.csv', delimiter=',', unpack=True, skiprows=1)
 
    # get indices at which the voltage changes
    idx = np.where(V[:-1] != V[1:])[0] + 1
     
    # insert 0 at the beggining and len at the of array of indices
    idx = np.insert(idx, 0, 0)
    idx = np.append(idx, len(V))
     
    # get number of different voltages
    n_idx = len(idx) -1
     
    # array to hold applied voltages and error in V
    V_avg = np.zeros(n_idx)
    V_avg_err = np.zeros(n_idx)
     
    # array to hold average photocurrent at each voltage and error in A
    I_avg = np.zeros(n_idx)
    I_avg_err = np.zeros(n_idx)
     
    # loop over all indices at which voltage changes
    for j in range(0, n_idx):
         
        # array of repeated voltages and current measurements
        V_rpt = V[idx[j] : idx[j+1]]
        I_rpt = I[idx[j] : idx[j+1]]
 
        # get average voltage note all voltages should be the same
        V_avg[j] = np.mean(V_rpt)
        V_avg_err[j] = np.abs(V_avg[j] * psu_err_var) + psu_err_const
         
        # get average and sem of repeated photocurrent measurements
        I_avg[j] = np.mean(I_rpt)
        I_avg_err[j] = np.std(I_rpt) / np.sqrt(len(I_rpt))
         
    # write arrays with averages into arrays of objects
    V_arr[i] = V_avg
    V_err[i] = V_avg_err
    I_arr[i] = I_avg
    I_err[i] = I_avg_err

# plot IV curve for each filter
for i in range(0, len(filt)):   
    
    # set grid parameters
    plt.rc('grid', linestyle=':', color='black', alpha=0.8)

    # start subplot with shared x axis
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax3 = ax2.twinx()
    
    # use different title for line curve and background curve
    if filt[i] != 'BG':
        fig.suptitle(f'{filt[i]} Line IV Curve \n dataset: SerMes_{ds_num}_{filt[i]}')
    
    if filt[i] == 'BG':
        fig.suptitle(f'Background IV Curve \n dataset: SerMes_{ds_num}_{filt[i]}')
    
    # plot IV curve
    ax1.set_ylabel(r'photocurrent $I$ (A)')
    ax1.grid()
    ax1.plot(V_arr[i], I_arr[i], c='blue')
    
    # plot error in photocurrent
    ax2.set_xlabel(r'voltage $V$ (V)')
    ax2.set_ylabel(r'photocurrent SEM $\delta \: I$ (A)')
    ax3.set_ylabel(r'voltage uncertainty $\delta V$ (mV)')
    ax2.grid()
    ax2.plot(V_arr[i], I_err[i], c='blue')
    
    # plot error in applied voltage
    ax3.spines['left'].set_color('blue')
    ax3.spines['right'].set_color('darkorange')
    ax2.tick_params(axis='y', colors='blue')
    ax3.tick_params(axis='y', colors='darkorange')
    
    ax3.set_ylim(0, np.amax(V_err[i]* 1e3))
    ax3.plot(V_arr[i], V_err[i] * 1e3, c='darkorange')

    # save the plot
    plt.savefig(f'Output/IV_{ds_num}_{filt[i]}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    