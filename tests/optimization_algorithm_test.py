# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:52:00 2020

@author: Mark
"""
from model.model import Model
import matplotlib.pyplot as plt
import random
import numpy as np

#%%
filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\GasScattering\data\O2\Ag3d 4mbar O2 - EX321.txt'
filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
filename3 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\GasScattering\data\O2\Ag3d vacuum - EX321.txt'
''' First set-up the model.'''
m = Model()
m.loadScatterers(filename2) 
m.loadSpectrum(filename1)
m.loadSpectrum(filename3)
m.setCurrentScatterer('O2')
m.scattering_medium.setPressure(4)
m.scattering_medium.scatterer.inelastic_xsect = 0.003
m.scattering_medium.scatterer.norm_factor = 1
m.algorithm_id = 0
m.scattered_spectrum = m.loaded_spectra[0]
m.unscattered_spectrum = m.loaded_spectra[1]
m.scattering_medium.scatterer.loss_function.reBuild()
m.scatterSpectrum()

Scat = m.scattered_spectrum
Unscat = m.unscattered_spectrum
Sim = m.simulated_spectrum

plt.plot(Scat.lineshape)

#%%
def optimizePeakHeight(model):
    grad_rate = 0.01
    rate = 0.0001
    numbers = [-1,1]
    rand = random.choice(numbers)
    m = model
    inelastic_xsect = m.scattering_medium.scatterer.inelastic_xsect
    plus = inelastic_xsect + rand * rate
    
    m.scattering_medium.scatterer.inelastic_xsect = plus
    print(m.scattering_medium.scatterer.inelastic_xsect)
    m.scatterSpectrum()
    diff1 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms1 = np.sqrt(np.sum(np.square(diff1)))
    rms1 = rms1 / np.sum(m.simulated_spectrum.lineshape)
    print("rms1: " + str(rms1))
    
    minus = inelastic_xsect - rand * rate
    m.scattering_medium.scatterer.inelastic_xsect = minus
    print(m.scattering_medium.scatterer.inelastic_xsect)
    m.scatterSpectrum()
    diff2 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms2 = np.sqrt(np.sum(np.square(diff2)))
    rms2 = rms2 / np.sum(m.simulated_spectrum.lineshape)
    print("rms2: " + str(rms2))
    
    gradient = (rms2-rms1) / (minus-plus)
    print("gradient: " + str(gradient))
    
    new_val = inelastic_xsect - gradient * grad_rate
    m.scattering_medium.scatterer.inelastic_xsect = new_val
    print("new_val: " + str(new_val))
    m.scatterSpectrum()
    diff3 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms3 = np.sqrt(np.sum(np.square(diff3)))
    rms3 = rms3 / np.sum(m.simulated_spectrum.lineshape)
    
    print("rms3: " + str(rms3))
    return rms3

def iterateOnce(model):
    m = model
    learning_rate = 1
    Scat = m.scattered_spectrum
    L = m.scattering_medium.scatterer.loss_function
    component_variables = [comp.__dict__ for comp in L.components]
    numbers = [-1,1]
    step_sizes = {'width':0.005,'position':0.05,'intensity':1,
                  'B':100,'C':100,'D':10, 'Eg':0}
    limits = {'width':(0.01,10), 'position':(0,1000), 'intensity':(0,1000),
              'B':(0,1E+9), 'C':(0,1E+9), 'D':(0.01,1000), 'Eg':(0,0)}
    gradients = []
    for idx, comp in enumerate(component_variables):
        gradients += [{}]
        for key, value in comp.items():
            if step_sizes[key] != 0:
                
                rand = random.choice(numbers)
                plus = value + step_sizes[key] * rand
                L.editComponent(idx, {key: plus})
                m.scattering_medium.scatterer.loss_function = L
                m.scatterSpectrum()
                Sim = m.simulated_spectrum
                Diff1 = np.subtract(Scat.lineshape, Sim.lineshape)
                RMS1 = np.sqrt(np.sum(np.square(Diff1)))
                RMS1 = RMS1 / np.sum(Sim.lineshape)
                #print(RMS1)
                #print("plus: " + str(plus))
                
                minus = value - step_sizes[key] * rand
                L.editComponent(idx, {key: minus})

                m.scattering_medium.scatterer.loss_function = L
                m.scatterSpectrum()
                Sim = m.simulated_spectrum
                Diff2 = np.subtract(Scat.lineshape, Sim.lineshape)
                RMS2 = np.sqrt(np.sum(np.square(Diff2)))
                RMS2 = RMS2 / np.sum(Sim.lineshape)
                #print(RMS1)
                #print("minus: " + str(minus))
                    
                Gradient = ((RMS2-RMS1) / (minus-plus))
                #print("gradient: " + str(Gradient))
                
                
            else:
                Gradient = 0
            gradients[-1][key] = Gradient

            new_val = np.abs(value - Gradient * learning_rate)
            if (new_val > limits[key][0]) and (new_val < limits[key][1]):
                L.editComponent(idx, {key: new_val})
            else:
                old_val = component_variables[idx][key]
                L.editComponent(idx, {key: old_val})
            m.scattering_medium.scatterer.loss_function = L
            m.scatterSpectrum()
            
    diff3 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms3 = np.sqrt(np.sum(np.square(diff3))) #/ np.sum(Sim.lineshape)
    
    return rms3 

def optimizeFactor(model):
    grad_rate = 2
    rate = 0.1
    numbers = [-1,1]
    rand = random.choice(numbers)
    m = model
    norm_factor = m.scattering_medium.scatterer.norm_factor
    plus = norm_factor + rand * rate
    
    m.scattering_medium.scatterer.norm_factor = plus
    print(m.scattering_medium.scatterer.norm_factor)
    m.scatterSpectrum()
    diff1 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms1 = np.sqrt(np.sum(np.square(diff1)))
    rms1 = rms1 / np.sum(m.simulated_spectrum.lineshape)
    print("rms1: " + str(rms1))
    
    minus = norm_factor - rand * rate
    m.scattering_medium.scatterer.norm_factor = minus
    print(m.scattering_medium.scatterer.norm_factor)
    m.scatterSpectrum()
    diff2 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms2 = np.sqrt(np.sum(np.square(diff2)))
    rms2 = rms2 / np.sum(m.simulated_spectrum.lineshape)
    print("rms2: " + str(rms2))
    
    gradient = (rms2-rms1) / (minus-plus)
    print("gradient: " + str(gradient))
    
    new_val = norm_factor - gradient * grad_rate
    m.scattering_medium.scatterer.norm_factor = new_val
    print("new_val: " + str(new_val))
    m.scatterSpectrum()
    diff3 = np.subtract(m.scattered_spectrum.lineshape, m.simulated_spectrum.lineshape)
    rms3 = np.sqrt(np.sum(np.square(diff3)))
    rms3 = rms3 / np.sum(m.simulated_spectrum.lineshape)
    
    print("rms3: " + str(rms3))
    return rms3
#%%        

rms1 = []
for n in range(200):
    rms1 += [optimizePeakHeight(m)]

plt.scatter(list(range(len(rms1))),rms1)

#%%

rms2 = []
for n in range(50):
    rms2 += [optimizeFactor(m)]

plt.scatter(list(range(len(rms2))),rms2)

#%%
rms = []
#%%
plt.plot(m.scattered_spectrum.lineshape[:-100])
I = 0
for i in range(I):
    rms += [iterateOnce(m)]

plt.plot(m.simulated_spectrum.lineshape[:-100])

plt.show()

#%%
plt.scatter(list(range(len(rms))),rms)

#%%
plt.plot(m.scattering_medium.scatterer.loss_function.lineshape[:-11500])


#%%
L = m.scattering_medium.scatterer.loss_function
component_variables = [comp.__dict__ for comp in L.components]


m.updateScatterersDict()
m.saveScatterers('optimizedHe')
