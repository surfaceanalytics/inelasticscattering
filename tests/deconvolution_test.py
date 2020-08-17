# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 11:48:22 2020

@author: Mark
"""

from analysis import Analysis
import matplotlib.pyplot as plt
import numpy as np
import math
import copy

from model import Model
#%%

loss_fns = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
filename = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\He\ARM22\Ag3d vacuum EX340 ARM22.TXT'
M = Model()
M.loadScatterers(loss_fns)
M.loadSpectrum(filename)
M.unscattered_spectrum = M.loaded_spectra[0]
M.algorithm_id = 0
M.scattering_medium.setPressure(4)
M.setCurrentScatterer('He')
M.scattering_medium.d_through_gas = 1000
M.scattering_medium.scatterer.norm_factor = 1
M.scatterSpectrum()
#%%
in_spec = M.simulation.P
loss = np.roll(np.flip(M.simulation.L) * 1, 0)

in_spec_padded = np.pad(in_spec, (len(loss)-len(in_spec),0), 'constant', constant_values = 0)

fft_in_spec = np.fft.fft(in_spec_padded)
fft_loss = np.fft.fft(loss)
fft_c = np.multiply(fft_in_spec, fft_loss)

ifft_c = np.fft.ifft(fft_c)
plt.plot(ifft_c[-900:])
plt.show()

decon = np.divide(fft_c, fft_loss)
ifft_dec = np.fft.ifft(decon)
plt.plot(ifft_dec[-900:])

fft_loss2 = np.multiply(fft_loss, fft_loss)
fft_loss2a = np.power(fft_loss, 2)

conv2 = np.multiply(fft_loss2a, fft_in_spec)
plt.plot(np.fft.ifft(conv2)[-900:])

comps = M.simulation.I.T

conv = M.simulation.simulated

fft_conv = np.fft.fft(conv)

total = np.multiply(fft_in_spec, np.exp(2*fft_loss))
ifft_total = np.fft.ifft(total)
plt.plot(ifft_total[-900:])



recover = np.divide(total, np.exp(2*fft_loss))
plt.plot(np.fft.ifft(recover)[-900:])
#%%

L = np.flip(M.simulation.L)

fft_loss = np.fft.fft(loss)

I = np.array([np.zeros(len(loss))])

n = 4
for i in range(0,n+1):
    new = np.power(fft_loss,i)
    new = np.array([new])
    I = np.concatenate((I,new),axis=0)
I = np.sum(I, axis=0)
I = np.fft.ifft(I)
plt.plot(I[-1000:])
