# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:37:57 2020

@author: Mark
"""
from gasscattering.controller.controller import Controller
from gasscattering.model.model_classes import Spectrum, ScatteringMedium, MeasuredSpectrum, Scatterer, LossFunction, Peak, VacuumExcitation
import matplotlib.pyplot as plt
import numpy as np
#%%
spec = Controller()
spec.loadPrimarySpectrum('Au4f_0mbarO2.txt')
spec.loadScatterer('O2')
spec.loadSpectrumToFit('Au4f_2mbarO2.txt')

#%%

#spec.scattering_medium.scatterer.loss_function.components[-1].fermi_edge = 19
#spec.scattering_medium.scatterer.loss_function.components[-1].fermi_width = 2
spec.scattering_medium.scatterer.loss_function.components[-1].gauss_mean = 0
spec.scattering_medium.scatterer.loss_function.components[-1].gauss_fwhm = 65
spec.scattering_medium.scatterer.loss_function.components[-1].intensity = 450

spec.scattering_medium.scatterer.loss_function.reBuild()

#%%
spec.scattering_medium.d_through_gas = 800000
spec.scattering_medium.scatterer.setCrossSec(0.3)
spec.scattering_medium.setPressure(2)
spec.scatterSpectrum()

n = 400
m = 2800
x = [i-1400 for i in spec.spectrum_to_fit.x]
#plt.plot(spec.spectrum_to_fit.lineshape[n:m] / np.max(spec.spectrum_to_fit.lineshape), color = 'blue')
#plt.plot(spec.output_spectrum.lineshape[n:m] / np.max(spec.output_spectrum.lineshape), color = 'red')
plt.plot(spec.spectrum_to_fit.lineshape[n:m], color = 'blue')
plt.plot(spec.output_spectrum.lineshape[n:m], color = 'red')

plt.show()


loss = spec.scattering_medium.scatterer.loss_function.lineshape

plt.plot(loss)
#%%

spec.exportExcel('Co3O4 in O2')

#%%

