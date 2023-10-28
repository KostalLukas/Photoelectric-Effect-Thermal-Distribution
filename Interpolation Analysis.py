# -*- coding: utf-8 -*-
"""
Interpolation Plancks Constant Analysis v2.1

Lukas Kostal, 25.10.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.constants as sc
import scipy.interpolate as si
import scipy.odr as sr


# linear function for fitting
def lin(x, m, c):
    y = m*x + c
    return y


# linear function for fitting with ODR
def srlin(q, x):
    m, c = q
    y = m*x + c
    return y


# specify dataset number to analyze
ds_num = 1
 
# specify available filters
filt = ['UV', 'V', 'B', 'G', 'Y', 'R']

# specify any filters to ignore in final linear fit
ign = ['R']

# constant error in power supply in V
psu_err_const = 0.005

# variable error in power supply as fraction of output voltage
psu_err_var = 5e-3
 
# specify voltage range on which to fit reverse current in V
Vi = -8
Vf = -4

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

# arrays to hold stopping potential and associated error
V0_arr = np.zeros(len(filt))
V0_err = np.zeros(len(filt))

# loop over all filters to fit cubic splines
for i in range(0, len(filt)):
    
    # fit cubic spline onto IV curve
    cs = si.CubicSpline(V_arr[i], I_arr[i])
    
    # numerically solve for zero of cubic spline to get stopping potential
    # -ve sign to account for -ve in electron charge in future
    V0_arr[i] = -np.amax(so.fsolve(cs, x0=0))
    
    # matrix of signs for combining errors for error propagation
    sg = np.array([[1,1], [1,-1], [-1, 1], [-1, -1]])
    
    # array to hold deviated cubic spline stopping potentials for erro propagation
    V0_dev = np.zeros(4)
    
    # loop over each way of combining the errors interpolate and get V0 for each
    for j in range(0, 4):
        cs_dev = si.CubicSpline(V_arr[i] + sg[j,0] * V_err[i], I_arr[i] + sg[j,1] * I_err[i])
        V0_dev[j] = np.amax(so.fsolve(cs_dev, x0=0))
        
    # take the greatest difference in the V0 from interpolations shifted by possible errors
    V0_err[i] = np.abs(np.ptp(V0_dev))

# create .csv file to store stopping potentials
with open('Output/Stopping.csv', 'w') as stp:
    
    # print headers to the file
    print('filter, stop pot V0 (V), error (V)', file=stp)
    
    # loop over each filter and print filter type, stopping potential and uncertainty in V
    for i in range(0, len(filt)):
        print(f'{filt[i]}, {V0_arr[i]}, {V0_err[i]}', file=stp)


# pick out filters to be ignored in final fit
ign_idx = filt.index(ign)
f_fit = np.delete(f, ign_idx)
f_fit_err = np.delete(f_err, ign_idx)
V0_fit = np.delete(V0_arr, ign_idx)
V0_fit_err = np.delete(V0_err, ign_idx)

# preliminary linear curve fit
popt, pcov = so.curve_fit(lin, f_fit, V0_fit, sigma=V0_fit_err, absolute_sigma=True)

# perform orthogonal distance regression with initial guess from preliminary fit
pdat = sr.RealData(f_fit, V0_fit, sx=f_fit_err, sy=V0_fit_err)
pmod = sr.Model(srlin)
preg = sr.ODR(pdat, pmod, beta0=popt)
pout = preg.run()
popt = np.array(pout.beta)
perr = np.array(pout.sd_beta)

# use optimised gradient to calcualte Planks constant and associated uncertainty
h_opt = popt[0] * sc.e
h_err = perr[0] * sc.e

# print the numerical value
print(f'h = {h_opt*1e34:.2f} ± {h_err*1e34:.2f} e-34 Js')
    
# parameters for plotting stopping potential against frequency
plt.figure(1)
plt.title('Stopping Potential against Frequency')
plt.xlabel(r'frequency $\nu$ (Hz)')
plt.ylabel(r'stopping potential $V_0$ (V)')
plt.rc('grid', linestyle=':', color='black', alpha=0.8)
plt.grid()

plt.errorbar(f_fit, V0_fit, xerr=f_fit_err, yerr=V0_fit_err, fmt='', ls='', capsize=4, c='blue')
plt.plot(f, lin(f, *popt), c='darkorange', zorder=0)
#plt.errorbar(f[-1], V0_arr[-1], xerr=f_err[-1], yerr=V0_err[-1], fmt='', ls='', capsize=4, c='red')

# save the plot
plt.savefig('Output/Stopping.png', dpi=300, bbox_inches='tight')
plt.show()    
 
 