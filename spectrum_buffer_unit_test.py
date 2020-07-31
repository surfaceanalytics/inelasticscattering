# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 09:30:47 2020

@author: Mark
"""

import numpy as np
import matplotlib.pyplot as plt
from base_model import (Spectrum, Gauss, Lorentz, MeasuredSpectrum, Voigt,
                        SyntheticSpectrum, LineShapeBuffer)

#%% Test adding and deleting components        
synth = SyntheticSpectrum(0,10,0.1)
synth.addComponent(Gauss(2,0.1,5))
synth.addComponent(Gauss(4,0.1,5))
synth.addComponent(Gauss(6,0.1,5))
synth.addComponent(Gauss(8,0.1,5))


plt.plot(synth.lineshape)
synth.removeComponent(2)
plt.plot(synth.lineshape)
synth.removeComponent(1)
plt.plot(synth.lineshape)
synth.removeComponent(1)
plt.plot(synth.lineshape)
synth.removeComponent(0)
plt.plot(synth.lineshape)
synth.addComponent(Gauss(2,0.1,5))
synth.addComponent(Gauss(4,0.1,5))
synth.addComponent(Gauss(6,0.1,5))
synth.addComponent(Gauss(8,0.1,5))
plt.plot(synth.lineshape)
synth.removeComponent(0)
synth.removeComponent(0)
synth.removeComponent(0)
plt.plot(synth.lineshape)


#%% Test modifying components

synth = SyntheticSpectrum(0,10,0.1)
synth.addComponent(Gauss(2,0.1,5))
comps = synth.components
params = {'width':1,'position':5,'intensity':10}
synth.editComponent(0, params)

arr = synth.buffer.array
l = synth.lineshape
