# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:27:57 2020

@author: Mark
"""
from gasscattering.model.model_classes import Scatterer
import pickle 
import os

datapath = os.path.dirname(os.path.abspath(__file__))

def read_scatterers(label):
    file = datapath+"\\scatterers"
    infile = open(file,'rb')
    scatterer = pickle.load(infile)[label]
    infile.close()
    return scatterer

default = Scatterer()
default.label = 'default'
default.cross_sec = 0.01
default.gas_diameter = 0.2
default.loss_function.addPeak(5,0.15,8.0)
default.loss_function.addVacuumExcitation(5,12,2,2,50)

He = Scatterer()
He.label = 'He'
He.cross_sec = 0.018
He.gas_diameter = 0.28
He.loss_function.addPeak(21.2,0.15,8.0)
He.loss_function.addPeak(23.03,0.15,2.2) 
He.loss_function.addPeak(23.67,0.15,0.9)
He.loss_function.addPeak(24.01,0.25,0.43)
He.loss_function.addPeak(24.2,0.25,0.33)
He.loss_function.addPeak(24.35,0.25,0.33)
He.loss_function.addVacuumExcitation(19,12,24.5,2,50)

N2 = Scatterer()
N2.label = 'N2'
N2.cross_sec = 0.04
N2.gas_diameter = 0.364
N2.loss_function.addPeak(9.1,1,1)
N2.loss_function.addPeak(12.95,0.35,5.8) # energyloss, broadening, probability
N2.loss_function.addPeak(14.15,0.35,3.8)
N2.loss_function.addPeak(15.95,0.50,3)
N2.loss_function.addPeak(17.6,0.45,1.5)
N2.loss_function.addPeak(18.6,0.45,1)
N2.loss_function.addVacuumExcitation(19,15,19,2,30)

O2 = Scatterer()
O2.label = 'O2'
O2.cross_sec = 0.148
O2.diameter = 0.346
O2.loss_function.addPeak(8.5,0.4,10.3)
O2.loss_function.addPeak(15.3,0.4,5.26)
O2.loss_function.addPeak(13.1,0.4,5)
O2.loss_function.addPeak(14.1,0.4,3)
O2.loss_function.addPeak(15.3,0.4,3.5)
O2.loss_function.addPeak(17,0.8,12)
O2.loss_function.addPeak(18.8,0.8,5.7)
O2.loss_function.addPeak(21.2,0.8,3)
O2.loss_function.addPeak(10,1.5,3)
O2.loss_function.addVacuumExcitation(10,20,19,2,230)

def write_scatterers():
    scatterers = {'default':default, 'He':He, 'N2':N2, 'O2':O2}
    file = datapath+"\\scatterers"
    outfile = open(file,'wb')
    pickle.dump(scatterers,outfile)
    outfile.close()

#%%
'''
spec.scattering_medium.scatterer.loss_function.components[0].mean = 8.5
spec.scattering_medium.scatterer.loss_function.components[0].stdev = 0.4
spec.scattering_medium.scatterer.loss_function.components[0].intensity = 10.3

spec.scattering_medium.scatterer.loss_function.components[1].mean = 15.3
spec.scattering_medium.scatterer.loss_function.components[1].stdev = 0.4
spec.scattering_medium.scatterer.loss_function.components[1].intensity = 5.26

spec.scattering_medium.scatterer.loss_function.components[2].mean = 13.1
spec.scattering_medium.scatterer.loss_function.components[2].stdev = 0.4
spec.scattering_medium.scatterer.loss_function.components[2].intensity = 5

spec.scattering_medium.scatterer.loss_function.components[3].mean = 14.1
spec.scattering_medium.scatterer.loss_function.components[3].stdev = 0.4
spec.scattering_medium.scatterer.loss_function.components[3].intensity = 3

spec.scattering_medium.scatterer.loss_function.components[4].mean = 15.3
spec.scattering_medium.scatterer.loss_function.components[4].stdev = 0.4
spec.scattering_medium.scatterer.loss_function.components[4].intensity = 3.5

spec.scattering_medium.scatterer.loss_function.components[5].mean = 17
spec.scattering_medium.scatterer.loss_function.components[5].stdev = 0.8
spec.scattering_medium.scatterer.loss_function.components[5].intensity = 12

spec.scattering_medium.scatterer.loss_function.components[7].mean = 18.8
spec.scattering_medium.scatterer.loss_function.components[7].stdev = 0.8
spec.scattering_medium.scatterer.loss_function.components[7].intensity = 5.7

spec.scattering_medium.scatterer.loss_function.components[6].fermi_edge = 19
spec.scattering_medium.scatterer.loss_function.components[6].fermi_width = 2
spec.scattering_medium.scatterer.loss_function.components[6].gauss_mean = 10
spec.scattering_medium.scatterer.loss_function.components[6].gauss_fwhm = 20
spec.scattering_medium.scatterer.loss_function.components[6].intensity = 230

spec.scattering_medium.scatterer.loss_function.components[8].mean = 21.2
spec.scattering_medium.scatterer.loss_function.components[8].stdev = 0.8
spec.scattering_medium.scatterer.loss_function.components[8].intensity = 3.0

spec.scattering_medium.scatterer.loss_function.components[9].mean = 10
spec.scattering_medium.scatterer.loss_function.components[9].stdev = 1.5
spec.scattering_medium.scatterer.loss_function.components[9].intensity = 3.0

spec.scattering_medium.scatterer.loss_function.reBuild()
'''
