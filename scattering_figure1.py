# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:57:44 2020

@author: Mark
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json
import re
from model import *
import matplotlib.pyplot as plt

#%%

model = Model()
filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\He\ARM22\Ag3d vacuum EX340 ARM22.TXT'
model.loadSpectrum(filename1)
filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
model.loadScatterers(filename2)

x = model.loaded_spectra[0].x
y = model.loaded_spectra[0].lineshape

def all_scattered(label):
    model.setCurrentScatterer(label)
    model.unscattered_spectrum = model.loaded_spectra[0]
    model.algorithm2()
    
    return [model.intermediate_spectra1['inel_scattered']]


#%%
    
# this makes a figure where the 
#scatter_labels =  ['He','H2','N2','O2']
scatter_labels =  ['He','N2']

p = 4

model.scattering_medium.setPressure(p)

ncol = 2
nrow = int(len(scatter_labels)/2)

fig, axs = plt.subplots(nrow, ncol, figsize=(5,2), dpi=200)
fig.tight_layout(pad=0.05, w_pad=0.1, h_pad=0.7)
for i, ax in enumerate(fig.axes):
    label = scatter_labels[i]
    data = all_scattered(label)[0]
    data = np.divide(data,1000)
    for j in data:
        ax.plot(j)
    ax.set_xlabel('Energy [eV]', fontsize=8)
    ax.set_ylabel('Intensity [kcps]', fontsize=8)
    ax.tick_params(direction='out', length=2, width=1, colors='black',
           grid_color='black', labelsize=6, grid_alpha=0.5)
    text_y = 0.85 * (max(ax.get_ylim()))
    text_x = 0.05 * (max(ax.get_xlim()))
    spl = re.findall(r'(\w+?)(\d+)', label)
    if len(spl)==0:
        pass
    else:
        spl = re.findall(r'(\w+?)(\d+)', label)[0]
        label = spl[0] + '$' + '_' + spl[1] + '$'

    ax.text(text_x, text_y, label, fontsize=12)
plt.show()



#%%
'''    
trends = []
for i in range(50):
    f = i * 10000
    model.scattering_medium.d_through_gas = f
    model.algorithm2()
    trends += [model.inel_probs]
    plt.plot(trends[-1])
plt.show()

totals = []
for i in range(len(trends)):
    x = 0
    for j in range(len(trends[i])):
        x += trends[i][j]
    totals += [x]
        

plt.plot(totals)

s1 = np.sum(model.el_probs[1:])
s2 = np.sum(model.inel_probs[1:])
'''


    