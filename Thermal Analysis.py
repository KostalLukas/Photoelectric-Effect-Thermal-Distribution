# -*- coding: utf-8 -*-
"""
Thermal Distribution Analysis v1.1

Lukas Kostal, 26.10.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.constants as sc
import scipy.odr as sr
import scipy.signal as ss


# linear function for fitting
def lin(x, m, c):
    y = m*x + c
    return y


# linear function for fitting with ODR
def srlin(q, x):
    m, c = q
    y = m*x + c
    return y


# function to apply low pass Butterworth filter
def butter_lp(arr, fs, fc, order):
    w = fc / fs * 2
    b, a = ss.butter(order, Wn=w, btype='low', analog=False)
    arr = ss.filtfilt(b, a, arr)
    return arr


# function of Fermi-Dirac distribution
def FermiDirac(x, Ef, T):
    y = 1 / (np.exp((x - Ef)/(sc.k * T)) + 1)
    return y


# function to get indices of rising edge above some threshold
def edge(arr, th):
    sgn = arr >= th
    idx = np.argwhere(np.convolve(sgn, [1, -1]) == 1)
    return idx


# specify dataset number to analyze
ds_num = 1

# approximate temp in K for theoretical distribution
T = 300

# fermi level of Potassium emitter in eV
Ef = 2.12

# voltage below which IV curve is ignored when finding rising edge in V
Vign = -8
 
# specify available filters
filt = ['UV', 'V', 'B', 'G', 'Y', 'R']

# constant error in power supply in V
psu_err_const = 0.005

# variable error in power supply as fraction of output voltage
psu_err_var = 5e-3

# specify voltage range on which to fit reverse current in V
Vi = -8
Vf = -4
 
# specify cutoff period in voltage increment for low pass filter in V
Vc = 0.5

# load wavelength and expected uncertainty in nm
wl, wl_err = np.loadtxt('Output/Peaks.csv', usecols=(1,2), delimiter=',', unpack=True, skiprows=1)

# calculate frequency in Hz
f = sc.c / wl * 1e9
f_err = f * wl_err / wl
 
# arrays to hold arrays of averaged voltage for each dataset
V_arr = np.zeros(len(filt), dtype=object)
V_err = np.zeros(len(filt), dtype=object)
I_arr = np.zeros(len(filt), dtype=object)
I_err = np.zeros(len(filt), dtype=object)
 
# load the data from background curve to subtract
V_bg, I_bg = np.loadtxt(f'Data/SerMes_{ds_num}_BG.csv', delimiter=',', unpack=True, skiprows=1)
 
# loop over all filters to subtract background and average over repeated measurements
for i in range(0, len(filt)):
     
    # load the data voltage V in V photocurrent I in nA and convert to A
    V, I = np.loadtxt(f'Data/SerMes_{ds_num}_{filt[i]}.csv', delimiter=',', unpack=True, skiprows=1)
 
    # subtract background photocurrent if specified
    I -= I_bg
 
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
        V_avg_err[j] = V_avg[j] * psu_err_var + psu_err_const
         
        # get average and sem of repeated photocurrent measurements
        I_avg[j] = np.mean(I_rpt)
        I_avg_err[j] = np.std(I_rpt) / np.sqrt(len(I_rpt))
         
    # write arrays with averages into arrays of objects
    V_arr[i] = V_avg
    V_err[i] = V_avg_err
    I_arr[i] = I_avg
    I_err[i] = I_avg_err

# loop over all filters to perform linear fit on reverse current and subtract
for i in range(len(filt)):
    
    # get increment in voltage per index
    Vinc = np.mean(np.diff(V_arr[i]))
    
    # get indices for Vi and Vf
    ni = int((Vi - V_arr[i][0]) / Vinc)
    nf = int((Vf - V_arr[i][0]) / Vinc)
    
    # slice voltage and current arrays to desired range for fitting
    V_lin = V_arr[i][ni:nf]
    I_lin = I_arr[i][ni:nf]
    V_sig = V_err[i][ni:nf]
    I_sig = I_err[i][ni:nf]
    
    # perform preliminary fit with usual curve fit
    lopt, lcov = so.curve_fit(lin, V_lin, I_lin, sigma=I_sig, absolute_sigma=True)

    # perform orthogonal distance regression with initial parameters from preliminary fit
    ldat = sr.RealData(V_lin, I_lin, sx=V_sig, sy=I_sig)
    lmod = sr.Model(srlin)
    lreg = sr.ODR(ldat, lmod, beta0=lopt)
    lout = lreg.run()
    lopt = np.array(lout.beta)
    lerr = np.array(lout.sd_beta)
    
    # subtract expected reverse current
    I_arr[i] -= lin(V_arr[i], *lopt)
    
    # modify expected error in current to account for error from reverse current subtraction
    I_err[i] = np.sqrt(I_err[i]**2 + V_arr[i]**2 * lerr[0]**2 + lerr[1]**2)


# arrays to hold differentiated photocurrent dI and theoretical distribution dI_fd
# note dV is not differentiated voltage just flipped voltage array corresponding to dI
dV_arr = np.zeros(len(filt), dtype=object)
dI_arr = np.zeros(len(filt), dtype=object)
dI_fd = np.zeros(len(filt), dtype=object)

# loop over all filters to analyze thermal distribution
for i in range(0, len(filt)):
    
    # get increment in voltage per index
    Vinc = np.mean(np.diff(V_arr[i]))
    
    # index below which voltage is to be ignored in finding rising edge 
    ign = int((Vign - V_arr[i][0]) / Vinc)
    
    # sliced arrays to exclude region to be ignored
    V = V_arr[i][ign:]
    I = I_arr[i][ign:]
    
    # apply low pass butterworth filter with specified critical period
    I = butter_lp(I, 1/0.01, 1/Vc, 1)
    
    # differentiate photocurrent wrt applied voltage and crop voltage array to match len
    dV = V[:-1]
    dI = np.diff(I) / np.diff(V)
    
    # normalize the differentiated photocurrent
    dI *= 1 / np.amax(dI)

    # get index of rising edge at height of 0.5
    idx = edge(dI, 0.5)[0,0]
    
    # shift the index to match rising edge to Fermi energy of emitter
    idx  += int(Ef/Vinc)
        
    # offset the voltage and differentiated photocurrent to match Fermi energy
    dV = dV[:idx] - dV[idx]
    dI = dI[:idx]
    
    # flip the curve in horizontal direction by taking -ve of applied voltage
    dV = -dV

    # get theoretical energy distribution from Fermi-Dirac with DOS function and normalize
    fd = FermiDirac(dV * sc.e, Ef * sc.e, T) * np.sqrt(dV)
    fd *= 1 / np.amax(fd)

    # write experimental and theoretical distribution arrays to array
    dV_arr[i] = dV
    dI_arr[i] = dI
    dI_fd[i] = fd
    
# calcualte energy of photon with expected error in eV
E_gam = sc.h * f / sc.e
E_err = E_gam * f_err / f

# parameters for plotting grid
plt.rc('grid', linestyle=':', color='black', alpha=0.8)

# prepare for plotting subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, constrained_layout=True)

# grid for each subplot
ax1.grid()
ax2.grid()
ax3.grid()

# title and axis labels for figure
fig.suptitle('Thermal Energy Distribution of Electrons', x=0.5, y=1.06, fontsize=14)
fig.supxlabel('electron energy $\epsilon$ (eV)')
fig.supylabel('probability measure $P(\epsilon)$ (a.u.)')

# plot distribution for green line
ax1.set_title(f'G Line $E_\gamma$ = {E_gam[3]:.3f}({(E_err[3]*1e3):.0f}) eV')
ax1.plot(dV_arr[3], dI_arr[3], c='blue', label='experimental')
ax1.plot(dV_arr[3], dI_fd[3], c='darkorange', label='theoretical')
ax1.axvline(Ef, ls=':',c='red', label='Fermi level $E_F$')

# plot distribution for yellow line
ax2.set_title(f'Y Line $E_\gamma$ = {E_gam[4]:.3f}({(E_err[4]*1e3):.0f}) eV')
ax2.plot(dV_arr[4], dI_arr[4], c='blue')
ax2.plot(dV_arr[4], dI_fd[4], c='darkorange')
ax2.axvline(Ef, ls=':',c='red')

# plot distribution for red line
ax3.set_title(f'R Line $E_\gamma$ = {E_gam[5]:.3f}({(E_err[5]*1e3):.0f}) eV')
ax3.plot(dV_arr[5], dI_arr[5], c='blue')
ax3.plot(dV_arr[5], dI_fd[5], c='darkorange')
ax3.axvline(Ef, ls=':',c='red')

# set x axis limit and legend
ax1.set_xlim(0, 6)
ax1.legend(loc=(0.1, 1.4), ncol=3)

# save the figure
plt.savefig('Output/Thermal.png', dpi=300, bbox_inches='tight')
plt.show()


    
    
    
    