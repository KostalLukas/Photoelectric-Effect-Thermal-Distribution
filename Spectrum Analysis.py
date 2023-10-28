# -*- coding: utf-8 -*-
"""
Hg Lamp Spectrum Analysis v2.1

Lukas Kostal, 23.10.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss


# specify measured filters
filt = ['UV', 'V', 'B', 'G', 'Y', 'R']

# speciy expected unceertainty of spectrometer as fwhm value in nm
fwhm_spec = 2

# specify weather to plot subplots or individual plots
SPLT = True

# arrays to hold peak wavelength and expected uncertainty in nm
wl_avg = np.zeros(len(filt))
wl_err = np.zeros(len(filt))

# set figure parameters for plotting subplots
if SPLT == True:
    plt.rc('grid', linestyle=':', color='black', alpha=0.8)
    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    
    fig.suptitle('Hg discharge lamp line spectra')
    fig.supxlabel('wavelength $\lambda$ (nm)')
    fig.supylabel('normalised intensity $I$ (a.u.)')

# loop over spectra for all filters
for i in range(0, len(filt)):
    
    # load the data 
    wl, A = np.loadtxt(f'Data/Spectra/{filt[i]}_spec.csv', delimiter=',', unpack=True, skiprows=1)
    
    # normalise the intensity
    A /= np.amax(A)
    
    # find all peaks above 0.9 normalised intensity
    pks, _ = ss.find_peaks(A, height=0.9)
    
    # calcualte average wavelength of found peaks
    wl_avg[i] = np.mean(wl[pks])
    
    # get indices of the edges of the peaks at 0.5 normalised intensity
    idx = ss.peak_widths(A, pks, rel_height=0.5)[-2:]

    # get minimum and maximum wavelengths at 0.5 normalised intensity
    wl_min = wl[int(np.floor(np.amin(idx)))]
    wl_max = wl[int(np.ceil(np.amax(idx)))]
    
    # get FWHM in nm
    fwhm = wl_max - wl_min
    
    # calulcate Gaussian like error from FWHM
    # combine in quadrature with fwhm error of spectrometer
    wl_err[i] = np.sqrt(fwhm**2 + fwhm_spec**2) / np.sqrt(8 * np.log(2))
    
    # plot individual spectra
    if SPLT == False:
        # parameters for plotting
        plt.figure(1)
        plt.title(f'{filt[i]} Line Spectrum \n $\lambda$ = {wl_avg[i]:.1f} ± {wl_err[i]:.1f} nm')
        plt.xlabel('wavelength $\lambda$ (nm)')
        plt.ylabel('normalised intensity $I$ (a.u.)')
        plt.rc('grid', linestyle=':', color='black', alpha=0.8)
        plt.grid()
        
        # plot the spectrum, FWHM and peak wavelength with error
        plt.plot(wl, A, c='blue')
        plt.errorbar(wl_avg[i], 1, xerr=wl_err[i], fmt='x', capsize=5, c='green')
        plt.hlines(0.5, wl_min, wl_max, ls='-', color='red')
        plt.axvline(wl_avg[i], ls=':', color='red')
        plt.xlim(wl[pks[0]]-10, wl[pks[0]]+10)
        
        # save the plot
        plt.savefig(f'Output/Spec_{filt[i]}.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    # plot all peaks on a subplot
    if SPLT == True:
        # subplot title and grid
        plt.subplot(2,3, i+1)
        plt.title(f'{filt[i]} Line \n $\lambda$ = {wl_avg[i]:.1f} ± {wl_err[i]:.1f} nm')
        plt.rc('grid', linestyle=':', color='black', alpha=0.8)
        plt.grid()
    
        # plot peak, FWHM and error on wavelength
        plt.plot(wl, A, c='blue')
        plt.errorbar(wl_avg[i], 1, xerr=wl_err[i], fmt='x', capsize=5, c='green')
        plt.hlines(0.5, wl_min, wl_max, ls='-', color='red')
        plt.axvline(wl_avg[i], ls=':', color='red')
        plt.xlim(wl[pks[0]]-10, wl[pks[0]]+10)
        plt.ylim(0, 1.1)

# save the subplots
if SPLT == True:
    plt.savefig('Output/Spectra.png', dpi=300, bbox_inches='tight')
    plt.show()
    
# create .csv file to store peaks
with open('Output/Peaks.csv', 'w') as pks:
    
    # print headers to the file
    print('filter, wavelength (nm), error (nm)', file=pks)
    
    # loop over each filter and print filter type, wavelength and error in nm
    for i in range(0, len(filt)):
        print(f'{filt[i]}, {wl_avg[i]}, {wl_err[i]}', file=pks)