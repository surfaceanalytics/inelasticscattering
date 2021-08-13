# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 09:50:10 2020

@author: Mark
"""
import numpy as np

class PathlengthApproximator():
    def __init__(self):
        self.formula = Formula()
        
    def function(self,x):
        return self.formula.function(x)
    

class Formula():
    def __init__(self):
        self.params = {'y0':1.84636,
                       'A1':-0.49188,
                       't1':0.19019,
                       'A2':-0.34152,
                       't2':0.57686}
        
    def function(self, x):
        y0 = self.params['y0']
        A1 = self.params['A1']
        t1 = self.params['t1']
        A2 = self.params['A2']
        t2 = self.params['t2']
        
        y = x * (y0 + A1*(1-np.exp(-x/t1)) + A2*(1-np.exp(-x/t2)))
        return y
            