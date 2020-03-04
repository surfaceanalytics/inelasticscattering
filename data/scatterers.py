# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:27:57 2020

@author: Mark
"""
from model import Scatterer
import os
import json

#datapath = os.path.dirname(os.path.abspath(__file__))
datapath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data'

default = Scatterer()
default.label = 'default'
default.cross_sec = 0.01
default.gas_diameter = 0.2
default.loss_function.addPeak(5,0.15,8.0)
default.loss_function.addVacuumExcitation(1,2,2,50) #exponent, edge, fermi_width, intensity

He = Scatterer()
He.label = 'He'
He.cross_sec = 0.038
He.gas_diameter = 0.28
He.loss_function.addPeak(21.2,0.05,12.0)
He.loss_function.addPeak(23.03,0.05,5.0) 
He.loss_function.addPeak(23.67,0.05,1)
He.loss_function.addPeak(24.01,0.1,0.1)
He.loss_function.addPeak(24.2,0.1,0.3)
He.loss_function.addPeak(27,6,0.24)
He.loss_function.addVacuumExcitation(0.006,25,10,0.1) #exponent, edge, fermi_width, intensity

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
N2.loss_function.addVacuumExcitation(0.1,19,2,30)

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
O2.loss_function.addVacuumExcitation(0.1,22,2,230)

def write_scatterers():
    scatterers = {'default':default, 'He':He, 'N2':N2, 'O2':O2}
    file = datapath+"\\scatterers"
    outfile = open(file,'wb')
    pickle.dump(scatterers,outfile)
    outfile.close()
    
write_scatterers()


scatterers = [default, He, N2, O2]

def build_dict(scatterer):
    d = {'cross_section':scatterer.cross_sec,
         'gas_diameter':scatterer.gas_diameter, 
         'loss_function':[(lambda x: {'id':x, 
                                      'type':scatterer.loss_function.components[x].__class__.__name__, 
                                      'params':scatterer.loss_function.components[x].__dict__})(i) 
                                        for i in range(len(scatterer.loss_function.components))]}
    return d

scatterers_dict = {}
for scatterer in scatterers:
    scatterers_dict[scatterer.label] = build_dict(scatterer)
    
datapath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data'
filename = 'scatterers.json'

with open(datapath+filename, 'w') as json_file:
  json.dump(scatterers_dict, json_file, indent=4)
  
  
