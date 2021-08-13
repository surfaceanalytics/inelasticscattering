# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:40:20 2021

@author: Mark
"""
from converters.data_converter import DataConverter
from model.base_model import MeasuredSpectrum, ScatteringMedium
from model.algorithms import algorithm5
from model.model import Model
import matplotlib.pyplot as plt
import os
import sys

abspath = os.path.abspath(__file__)
directory = os.path.dirname(abspath)

directory = os.getcwd()

filename = os.path.join(directory, r'data\test_spectrum.vms')

converter = DataConverter()
converter.load(filename)
data = converter.data[0]
input_spectrum = MeasuredSpectrum(data['data']['x'], data['data']['y0'])

model = Model()
model.loaded_spectra+=[input_spectrum]

plt.plot(model.loaded_spectra[0].x, model.loaded_spectra[0].lineshape)


filename2 = os.path.join(directory, r'data\scatterers2.json')
model.loadScatterers(filename2)

model.unscattered_spectrum = model.loaded_spectra[-1]

#%%
model.scattering_medium.d_through_gas = 10000
model.scattering_medium.pressure = 40

model.setCurrentScatterer('Shirley')
model.changeAlgorithm(1)
model.calculation.n_iter = 50
model.algorithm_option = 'bulk'
model.scattering_medium.scatterer.norm_factor = 0.6

shirley_params = {'edge':0, 'exponent':0.05, 'fermi_width':0.01}
model.scattering_medium.scatterer.loss_function.editComponent(0, shirley_params)

model.scatterSpectrum()

fig, ax = plt.subplots(1,2)
ax[0].plot(model.simulated_spectrum.x, model.simulated_spectrum.lineshape)
ax[1].plot(model.scattering_medium.scatterer.loss_function.x,model.scattering_medium.scatterer.loss_function.lineshape)
ax[1].set_xlim(0,200)
