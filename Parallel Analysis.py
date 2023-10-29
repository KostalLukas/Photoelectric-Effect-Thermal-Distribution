# -*- coding: utf-8 -*-
"""
Parallelized Optimization Plancks Constant Analysis v2.0

Lukas Kostal, 24.10.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.constants as sc
import multiprocessing as mp
import functools as ft
import scipy.odr as sr
import warnings


# linear function for fitting
def lin(x, m, c):
    y = m*x + c
    return y


# linear function for fitting with ODR
def srlin(q, x):
    m, c = q
    y = m*x + c
    return y


# quadratic function for fitting
def quad(x, a, b, c):
    y = a*x**2 + b*x + b**2 + c
    return y

# quadratic function for fitting with ODR
def srquad(q, x):
    a, b, c = q
    y = a*x**2 + b*x + c
    return y


# function to calculate reduced chi squared
def chi(e, o, n):
    chi = np.sqrt( np.sum((o - e)**2) / (len(o) - n) )
    return chi


# function to subtract background and average repeated readings
# for all filters 
def average(ds_num, filt, psu_err_const, psu_err_var):
    
    # arrays to hold arrays of averaged voltage for each dataset
    V_arr = np.zeros(len(filt), dtype=object)
    V_err = np.zeros(len(filt), dtype=object)
    I_arr = np.zeros(len(filt), dtype=object)
    I_err = np.zeros(len(filt), dtype=object)
     
    # load the data from background curve to subtract
    V_bg, I_bg = np.loadtxt(f'Data/SerMes/SerMes_{ds_num}_BG.csv', delimiter=',', unpack=True, skiprows=1)
     
    # loop over all filters to subtract background and average over repeated measurements
    for i in range(0, len(filt)):
         
        # load the data voltage V in V photocurrent I in nA and convert to A
        V, I = np.loadtxt(f'Data/SerMes/SerMes_{ds_num}_{filt[i]}.csv', delimiter=',', unpack=True, skiprows=1)
     
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

    return V_arr, V_err, I_arr, I_err

# function to preform reduced chi squared optimisation
# for specified filter for parellelising
def optimize(i, Vign, mpx, V_arr, V_err, I_arr, I_err):
    
    # supress warnings from numpy and scipy curve fit
    warnings.filterwarnings('ignore')
    
    # get increment in voltage per index
    Vinc = np.mean(np.diff(V_arr[i]))
    
    # index of voltage at which to start optimisation
    ign = int((Vign - V_arr[i][0]) / Vinc)
    
    # initial and final indices for linear fit
    Ni =  ign + 2
    Nf =  len(V_arr[i]) - 3

    # note mpx is increment usually 10 in which the loop increments

    # dummy loop to calculate required size of arrays indexed with s
    s = 0
    for n in range(Ni, Nf, mpx):
        for m in range(n + 3, len(V_arr[i]), mpx):
            s += 1
    
    # arrays to hold stopping potential and associated error for given filter
    V0f_arr = np.zeros(s)
    V0f_err = np.zeros(s)
    
    # array to hold reduced chi squared values for given filter
    chif_arr = np.zeros(s)
    
    # array to hold n, m indices for all combinations of given filter
    nmf_arr = np.zeros((s, 2), dtype=int)
    
    # linear and quadratic models for orthogonal distance regression
    lmod = sr.Model(srlin)
    qmod = sr.Model(srquad)

    # loop over all possible indices for linear fit
    s = 0
    for n in range(Ni, Nf, mpx):
        
        # slice voltage, current and associated error arrays
        V_lin = V_arr[i][ign:n] 
        I_lin = I_arr[i][ign:n]
        V_lin_err = V_err[i][ign:n]
        I_lin_err = I_err[i][ign:n]
        
        # preform preliminary linear fit
        lopt, lcov = so.curve_fit(lin, V_lin, I_lin, sigma=I_lin_err, absolute_sigma=True)
        
        # orthogonal distance regression with initial parameters from preliminary fit
        ldat = sr.RealData(V_lin, I_lin, sx=V_lin_err, sy=I_lin_err)
        lreg = sr.ODR(ldat, lmod, beta0=lopt, maxit=0)
        lout = lreg.run()
        lopt = lout.beta
        lerr = lout.sd_beta
        
        # initial and final indices for quadratic fit
        Mi = n + 3
        Mf = len(V_arr[i])
        
        # loop over all possible indices for quadratic fit
        for m in range(Mi, Mf, mpx):
            
            # slice voltage, current and associated error arrays
            V_quad = V_arr[i][n:m]
            I_quad = I_arr[i][n:m]
            V_quad_err = V_err[i][n:m]
            I_quad_err = I_err[i][n:m]
            
            # subtract linear fit and account for increase in uncertainty
            I_quad =  I_quad - lin(V_quad, *lopt)
            I_quad_err = np.sqrt(I_quad_err**2 + V_quad**2 * lerr[0]**2 + lerr[1]**2)
            
            # perform preliminary qudratic fit
            qopt, qcov = so.curve_fit(quad, V_quad, I_quad, sigma=I_quad_err, absolute_sigma=True)
            
            # orthogonal distance regression with initial parameters from preliminary fit
            qdat = sr.RealData(V_quad, I_quad, sx=V_quad_err, sy=I_quad_err)
            qreg = sr.ODR(qdat, qmod, beta0=qopt, maxit=0)
            qout = qreg.run()
            qopt = np.array(qout.beta)
            #qerr = np.array(qout.sd_beta)
            
            # unpack optimised quadratic fit parameters and calcualte determinant
            a, b, c = qopt
            det = np.sqrt(b**2 - 4 * a * c)
            
            # calcualte derivatives for error propagation
            dVda = (c / (a * det) - (b - det) / (2 * a**2))**2
            dVdb = (1 / (2 * a) - b / det)**2
            dVdc = (1 / det)**2
            
            # calculate stopping potential from optimised parameters and associated error
            V0f_arr[s] = (b - det) / (2 * a)
            V0f_err[s] = np.sqrt(dVda * qcov[0,0] + dVdb * qcov[1,1] + dVdc * qcov[2,2])
            
            # calculate expectation current from linear and quadratic fit
            I_exp = quad(V_quad, *qopt) + lin(V_quad, *lopt)
            
            # get reduced chi sqaured value
            chif_arr[s] = chi(I_exp, I_quad, 4)
            
            # write indices n and m to array
            nmf_arr[s, 0] = n
            nmf_arr[s, 1] = m
            
            # increment the overall index
            s += 1
            
    return V0f_arr, V0f_err, chif_arr, nmf_arr

# function to write np.nan in reduced chi sqaured array for np.inf and np.nan in other arrays
def inftonan(filt, V0_arr, V0_err, chi_arr):
    
    # loop over all filters
    for i in range(0, len(filt)):
        
        # get boolean arrays of np.nan and np.inf in V0_arr
        # for each boolean write np.nan in chi_arr
        Bnan = np.isnan(V0_arr[i])
        Binf = np.isinf(V0_arr[i])
        chi_arr[i][np.logical_or(Bnan, Binf)] = np.nan
        
        # get boolean arrays of np.nan and np.inf in V0_err
        # for each boolean write np.nan in chi_arr
        Bnan = np.isnan(V0_err[i])
        Binf = np.isinf(V0_err[i])
        chi_arr[i][np.logical_or(Bnan, Binf)] = np.nan
        
        # get boolean array of np.inf in chi_arr and instead write np.nan
        Binf = np.isinf(chi_arr[i])
        chi_arr[i][Binf] = np.nan
    
    return chi_arr
    
# function to choose optimum stopping potential from minimum reduced chi squared value
def chooseopt(filt, V0_arr, V0_err, chi_arr):
    
    # arrays of optimised stopping potential and associated error for each filter
    V0_opt = np.zeros(len(filt))
    V0_opt_err = np.zeros(len(filt))
    
    # loop over each filter and get optimised stopping potential and associated error
    for i in range(0, len(filt)):
        idx = np.nanargmin(chi_arr[i])
        
        V0_opt[i] = V0_arr[i][idx]
        V0_opt_err[i] = V0_err[i][idx]
        
    return V0_opt, V0_opt_err

# run main analysis
if __name__ == "__main__":
    
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
    
    # load wavelength and expected uncertainty in nm
    wl, wl_err = np.loadtxt('Output/Peaks.csv', usecols=(1,2), delimiter=',', unpack=True, skiprows=1)
    
    # calculate frequency in Hz
    f = sc.c / wl * 1e9
    f_err = f * wl_err / wl
    
    # average over repeated measurements for all filters
    V_arr, V_err, I_arr, I_err = average(ds_num, filt, psu_err_const, psu_err_var)
    
    # voltage below which chi squared optimisation is not performed
    Vign = -8
        
    # increment for loops in readuced chi squared optimisation
    mpx = 10
    
    # print status update
    print('starting parallel optimisation')
    print()
    
    # input arguments for function to be parallelised
    parallelize = ft.partial(optimize, Vign=Vign, mpx=mpx, V_arr=V_arr, V_err=V_err, I_arr=I_arr, I_err=I_err)

    # create pool of parallel processes and run optimisation in parallel
    pool = mp.Pool(processes=len(filt))
    output = pool.map(parallelize, range(0, len(filt)))
    pool.close()
    
    # arrays of objects to assign output from parallel processes
    V0_arr = np.zeros(len(filt), dtype=object)
    V0_err = np.zeros(len(filt), dtype=object)
    chi_arr = np.zeros(len(filt), dtype=object)
    nm_arr = np.zeros(len(filt), dtype=object)
    
    # loop over all filters and assign output from parallel processes
    for i in range(0, len(filt)):
        V0_arr[i] = output[i][0]
        V0_err[i] = output[i][1]
        chi_arr[i] = output[i][2]
        nm_arr[i] = output[i][3]

    # write np.nan in reduced chi squared for any np.inf or np.nan in any of input arrays
    chi_arr = inftonan(filt, V0_arr, V0_err, chi_arr)

    # get arrays of optimised stopping potential and associated error in V
    V0_opt, V0_opt_err = chooseopt(filt, V0_arr, V0_err, chi_arr)
    
    # pick out filters to be ignored in final fit
    ign_idx = filt.index(*ign)
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
    
    # show plot but dont save
    plt.show()    
    
