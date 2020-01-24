# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:10:26 2020

@author: Mark
"""
import sys
import os

folder_name = "gasscattering"

path_name = os.path.dirname(
        os.path.abspath(__file__)).partition(folder_name)[0] + folder_name
        
if path_name not in sys.path:
    sys.path.insert(0, path_name)
    
else:
    pass