# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 09:25:11 2020

@author: Mark
"""
#%%
from model.base_model import (Spectrum, Gauss, Lorentz, VacuumExcitation, 
                        MeasuredSpectrum, ScatteringMedium, Calculation, 
                        Tougaard, Voigt, SyntheticSpectrum)

from model.algorithms.algorithm4 import Algorithm4
from model.model import Model

import matplotlib.pyplot as plt

filename = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\He\ARM22\Ag3d vacuum EX340 ARM22.TXT'
filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'

''' Instantiate Model object (needed for loading spectra and setting 
parameters)
'''
m = Model()
''' Load the scatterers from a JSON file.'''
m.loadScatterers(filename2)
''' Set the scatterer.'''
m.setCurrentScatterer('He')
''' Load a spectrum.'''
m.loadSpectrum(filename)
''' Set the unscattered spectrum to the loaded spectrum (this sets it to be the 
input spectrum of the algorithm).''' 
m.unscattered_spectrum = m.loaded_spectra[-1]
''' Set scattering medium parameters'''
m.scattering_medium.d_through_gas = 1000
m.scattering_medium.pressure = 4
m.scattering_medium.scatterer.inel_decay_factor = 1

''' Store the arguments to be used to initialize the Algorithm object as
local variables.'''
input_spect = m.unscattered_spectrum
scat_med = m.scattering_medium
params = {'option':''}

''' Instantiate Algorithm object and pass it the input spectrum and scattering 
medium as arguments.'''
a = Algorithm4(input_spect, scat_med, params)

''' Check to see if convolution is working correctly'''
a._convolve()
plt.plot(a.I.T)

''' Check to see that the Poisson factors are correctly calculated.'''
a._genPoisson()
plt.plot(a.poiss_inel)

''' Check to see if normalization factors are correctly calculated.'''
a._genNormFactors()
plt.plot(a.norm_factors)

''' Check to see that multiplication of Poisson factors and Norm factors is
working correctly.
'''
a._genScaleFactors()

''' Check the that thin film calculation is working correctly.'''
a._simulateFilm()
plt.plot(a.inel)
plt.plot(a.simulated)


''' Check the that thin film calculation is working correctly.'''
a._simulateBulk()
plt.plot(a.inel)
plt.plot(a.simulated)

