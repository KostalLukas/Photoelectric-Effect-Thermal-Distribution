# -*- coding: utf-8 -*-
"""
DAQ Measuring Sequence v2.0

Lukas Kostal, 15.10.2023, ICL
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
import sys
import glob as gb
from importlib import reload
import alive_progress as ap

import SerLib as sl
sl = reload(sl)


# function to pass arguments from terminal
def parg(*arg_var):
    arg_sys = sys.argv[1:]

    arg_name = []
    arg_type = []
    for i in range(0, len(arg_var)):
        arg_id = id(arg_var[i])

        for key in globals().keys():
            if not(key in arg_name or key[0] == '_'):
                val = globals()[key]
                if id(val) == arg_id:
                    arg_name.append(key)
                    arg_type.append(type(val))

    for i in range(0, len(arg_sys)):
        for j in range(0, len(arg_name)):
            if arg_sys[i].split('=')[0] == arg_name[j]:

                arg_val = arg_sys[i].split('=')[1]

                if arg_type[j] == bool:
                    arg_val = arg_val == 'True'
                 
                globals()[arg_name[j]] = arg_type[j](arg_val)
    return None


# function to get mean of every n elements in array
def get_avg(arr, n):
    n = int(n)
    arr = arr[(len(arr) % n):]
    arr = np.mean(arr.reshape(-1, n), axis=1)
    
    return arr


# function to get sem of every n elements in array
def get_sem(arr, n):
    n = int(n)
    arr = arr[(len(arr) % n):]
    std = np.std(arr.reshape(-1, n), axis=1)
    err = std / np.sqrt(n)
    
    return err


# function to print to console and file simultaneously
def tprint(text=''):
    print(text)
    global file
    with open(file, 'a') as output:
        print(text, file=output)
        
    return None


# specify serial ports for power supply and picoammeter
psu_port = '/dev/tty.usbserial-14230'
pam_port = '/dev/tty.usbserial-14220'

# minimum and maximum voltage and increment in V and no of repeat measurements
Vmin = -10
Vmax = 10
Vinc = 0.1
Nrpt = 10

# get dataset number
ds_num = input('Eneter dataset number: ')

# get =filter type
ds_filt = input('Enter filter type: ')

# generate dataset name
ds_name = f'SerMes_{ds_num}_{ds_filt}'

# get list of all dataset filepaths in the Data directory
ds_arr = gb.glob('Data/SerMes_**.csv')

# check if dataset already exists if yes ask to add _rt extension
if (f'Data/{ds_name}.csv' in ds_arr) == True:
    print('Dataset already exists')
    input('Press enter to continue with _rt extension')
    ds_name += '_rt'
    
# caclulate number of voltages to be set
Nvlt = int((Vmax - Vmin) / Vinc)

# no of total measurements
Ntot = int(Nrpt * Nvlt)

# get timestamp for the dataset
ts = datetime.datetime.now()
ts.replace(microsecond=0)

# write metadata file
file = f'Data/{ds_name}_md.txt'
open(file, 'w')

print()
tprint(f'Dataset: {ds_name}')
tprint(f'Timestamp: {ts}')
tprint()
tprint(f'Vmin = {Vmin:.2f} V')
tprint(f'Vmax = {Vmax:.2f} V')
tprint(f'Vinc = {Vinc:.2f} V')
tprint()
tprint(f'Nrpt = {Nrpt:.0f}')
tprint(f'Nvlt = {Nvlt:.0f}')
tprint(f'Ntot = {Ntot:.0f}')
print()

# initialise measurement instruments
psu = sl.psu(psu_port, 9600, bar=True)
pam = sl.pam(pam_port, 9600)

# ensure serial ports are open
psu.opn()
pam.opn()

# set timeout to 400ms
psu.timeout(200)
pam.timeout(200)

# check psu and pam connected with IDN? query
IDN_psu = psu.IDN()
IDN_pam = pam.IDN().split(',')[0]

# print connected instrument status
print()
print(f'PSU connected: {IDN_psu}')
print(f'PAM connected: {IDN_pam}')

# set psu output to 0V and enable
psu.VSET(0)
psu.OUT(True)

# reset the picoammeter, perform zero calibration, turn filter off
pam.RST()
pam.ZCAL()
pam.FOFF()

# set starting range in nA
rang = 2e-9
pam.RANG(rang)

# array of voltages to be measured
V = np.arange(Vmin, Vmax, Vinc)

# array to store measured photocurrents
I = np.zeros(len(V)*Nrpt)

# create dataset file
with open(f'Data/{ds_name}.csv', 'w') as ds:

    # create column names for the dataset
    print('voltage (V), photocurrent (A)', file=ds)

    # initialise a progress bar
    with ap.alive_bar(len(V) * Nrpt) as bar:
        
        # loop over all voltages
        for i in range(0, len(V)):
        
            # set psu voltage
            psu.VSET(V[i])
        
            # loop to take repeated measurements
            for j in range(0, Nrpt):    
            
                # calculate index of measurement
                idx = i * Nrpt + j    
            
                # take measuremenet and write to dataset file
                I[idx] = pam.READ()
                
                # check for overflow error
                while I[idx] >= 1e+37:
                    # increase range to next highest value and set it
                    rang *= 10
                    pam.RANG(rang)
                    
                    # take another reading
                    I[idx] = pam.READ()
                    
                # write the reading to the dataset
                print(f'{V[i]:.3f}, {I[idx]}', file=ds)
            
                # show current value in terminal and update progress bar
                bar.title(f'I = {I[idx]} A \t')
                bar()

# print status
print() 
print('Data acquisition finished')
print()

# set output voltage to 0V and disable output
psu.VSET(0)
psu.OUT(False)

# close the serial ports
psu.close()
pam.close()

# get mean and sem of measured current for preliminary IV curve
I_avg = get_avg(I, Nrpt)
I_sem = get_sem(I, Nrpt)

#print status
print('Plotting measured IV curve')
print()

# parameters for plotting preliminary IV curve
plt.figure(1)
plt.title(f'IV Curve \n dataset: {ds_name}')
plt.xlabel(r'voltage (V)')
plt.ylabel(r'photocurrent $I$ (A)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

# plot IV curve
plt.errorbar(V, I_avg, yerr=I_sem, fmt='x', c='blue', capsize=5)

plt.show()

