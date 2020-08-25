# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:32:48 2020

@author: Mark
"""
import numpy as np
from converters.data_converter import DataConverter
# Usage example        
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\XPS data analysis manuscript with Neal\C1s scans.xy'
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Data Science\xps_data_conversion_tools\EX337 - test.vms'
filepath = r'C:\Users\Mark\Desktop\Ir75Ru25.vms'
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\H2\Ag3p vacuum EX322.txt'
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\EX340_CEC356_Au O2 ARM22.vms'
filepath = r'C:\Users\Mark\Desktop\EX434 Ag in He.xy'
converter  = DataConverter()

converter.load(filepath)
d = converter.data
converter.write('test', out_format = 'Vamas')
converter.write('test', out_format = 'JSON')
converter.write('test', out_format = 'Excel')

for d in converter.data:
    if d['spectrum_type'] == 'Ag3d (2)':
        d['spectrum_type'] = 'Ag3d'

converter.data = converter.data[:-1]

types = list(set([d['spectrum_type'] for d in converter.data]))
lengths = {}
for t in types:
    for d in converter.data:
        if (d['spectrum_type'] not in lengths.keys()) & (d['spectrum_type'] == t):
            lengths[d['spectrum_type']] = list(np.zeros(len(d['data']['y0'])))
        else:
            continue
        
 
for t in types:
    temp = lengths[t]
    temp_time = 0
    for dat in converter.data:
        if dat['spectrum_type'] == t:
            if len(dat['data']['y0']) == len(temp):
                dat['data']['y0'] = list(np.add(np.array(dat['data']['y0']), np.array(temp)))
                dat['settings']['dwell_time'] = float(dat['settings']['dwell_time']) + temp_time
                temp = dat['data']['y0']
                temp_time = dat['settings']['dwell_time']
            

