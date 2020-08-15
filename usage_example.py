# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:32:48 2020

@author: Mark
"""

from data_converter import DataConverter
# Usage example        
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\XPS data analysis manuscript with Neal\C1s scans.xy'
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Data Science\xps_data_conversion_tools\EX337 - test.vms'
filepath = r'C:\Users\Mark\Desktop\Ir75Ru25.vms'
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\H2\Ag3p vacuum EX322.txt'
#filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\EX340_CEC356_Au O2 ARM22.vms'
data  = DataConverter()

data.load(filepath)
d = data.data
data.write('test', out_format = 'Vamas')
data.write('test', out_format = 'JSON')
data.write('test', out_format = 'Excel')

x = data.data[0]['data']['x']
