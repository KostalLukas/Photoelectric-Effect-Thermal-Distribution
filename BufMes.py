# -*- coding: utf-8 -*-
"""
DAQ Buffer Measuring Sequence v1.0

Lukas Kostal, 16.10.2023, ICL
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
import sys
import subprocess
import glob as gb
import time as tm
from importlib import reload

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
    std = np.std(arr.reshape(-1, n), axis=11)
    err = std / np.sqrt(n)
    
    return err


# function to print to console and file simultaneously
def tprint(text=''):
    print(text)
    global file
    with open(file, 'a') as output:
        print(text, file=output)
        
    return None


# function to convert time to s or min with units for printing
def smtim(time):
    if time > 60:
        return f'{(time/60):.0f} min'
    else:
        return f'{time:.0f} s'


# specify serial ports for power supply and picoammeter
psu_port = '/dev/tty.usbserial-14230'
pam_port = '/dev/tty.usbserial-14220'

# minimum and maximum voltage and increment in V and no of repeat measurements
Vmin = -10
Vmax = 10
Vinc = 0.01
Nrpt = 10

# delay time between repeated measurements in s
Tdly = 0

# integration time in number of line phase cycles
Nplc = 1

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
ts = ts.replace(microsecond=0)

# calculate integration time in ms
Tint = Nplc / 50 * 1000

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
tprint(f'Tdly = {Tdly:.0f} ms')
tprint(f'Nplc = {Nplc:.2f}')
tprint(f'Tint = {Tint:.0f} ms')
tprint()
tprint(f'Nrpt = {Nrpt:.0f}')
tprint(f'Nvlt = {Nvlt:.0f}')
tprint(f'Ntot = {Ntot:.0f}')
print()

# boolean to sepcify if using macOS
macOS = False

# call caffinate to prevent osx from entering sleep mode
if 'darwin' in sys.platform:
    print('macOS detected running caffinate subprocess')
    print()
    macOS = True
    subprocess.Popen('caffeinate')

# initialise measurement instruments
psu = sl.psu(psu_port, 9600)
pam = sl.pam(pam_port, 57600)

# ensure serial ports are open
psu.opn()
pam.opn()

# set timeout to 300 ms
psu.timeout(300)
pam.timeout(300)

# check psu and pam connected with IDN? query
IDN_psu = psu.IDN()
IDN_pam = pam.IDN().split(',')[0]

# print connected instrument status
print(f'PSU connected: {IDN_psu}')
print(f'PAM connected: {IDN_pam}')
print()

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

# prepare the picoammeter for reading into buffer
pam.BPREP(Nrpt, Tdly, Nplc)

# array of voltages to be measured
V = np.arange(Vmin, Vmax, Vinc)

# array to store measured photocurrents
I = np.zeros((len(V), Nrpt))

# create dataset file
with open(f'Data/{ds_name}.csv', 'w') as ds:

    # create column names for the dataset
    print('voltage (V), photocurrent (A)', file=ds)

    # get time at the start of the data acquisition
    Ti = tm.time()
    
    # loop over all voltages
    for i in range(0, len(V)):
        # set psu voltage
        psu.VSET(V[i])
    
        # prepare the picoammeter for next measurement into buffer
        pam.BARM()
        
        # take set of measurements into buffer and read it
        I_buff = pam.BREAD()
        
        # check for overflow error
        while np.amax(I_buff) >= 1e+37:
            # increase range to next highest value and set it
            rang *= 10
            pam.RANG(rang)
            
            # take measuremenet and write to dataset file
            pam.BARM()
            I_buff = pam.BREAD()
        
        # get average and sem of current measurements from the buffer
        I_avg = np.mean(I_buff)
        I_sem = np.std(I_buff) / np.sqrt(Nrpt)
        
        # write the measurements into the current array
        I[i, :] = I_buff
        
        # no of readings taken so far
        Nnow = (i+1) * Nrpt
        # elapsed time in s
        Tnow = tm.time() - Ti
        if Tnow < 60:
            Tnow_print = f'{(Tnow / 60):.0f} min'
        else:
            Tnow_print = f'{Tnow:.0f} s'
        
        
        # time per one reading in s
        Tpor = Tnow / Nnow
        
        
        # estimated time till end of data acquisition
        Test = (Ntot - Nnow) * Tpor / 60
    
        print(f'{Nnow}/{Ntot}, ', \
              f'{(Nnow/Ntot)*100:.0f}%, ', \
              f'{smtim(Tnow)}, ', \
              f'{smtim(Test)}, ', \
              f'{Tpor*1e3:.0f} ms')
        print(f'I = {I_avg:.6e} ± {I_sem:.6e} A \t')
        print('\033[A\033[A\033[A')
        
        # loop over measurements from the buffer and write to the dataset
        for j in range(0, Nrpt):
            print(f'{V[i]:.3f}, {I_buff[j]}', file=ds)
            
            
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
I_avg = np.mean(I, axis=1)
I_sem = np.std(I, axis=1) / np.sqrt(Nrpt)

#print status
print('Plotting measured IV curve')

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