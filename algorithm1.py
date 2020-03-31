# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 11:57:39 2020

@author: Mark
"""

from angular_spread_lorentzian import AngularSpreadCalc
import numpy as np

'''
This algorithm: 
1) iteratively convolves an imput lineshape with a loss function
2) re-scales the spectra by the inelastic and elastic probabilities for each 
    iteration 
    
Parameters
----------
n: the number of iterations to use
P: an array that represents the inpout spectrum
L: an array that represents the loss function
I: an array to hold the intermediate inelastically scattered spectra
    It has the shape (m, n), where m is the number of spectra and n is the number
    of data points per spectrum.
elastic_prob: a float to represent how the convolution thouls be re-scaled each
iteration
inelastic_xsect: a float to represent how the non-convolved portion of the 
spectrum should be scaled each iteration

inel_decay_factor: a float used to make the inelastic intensity decay with each
iteration
el_decay_factor: a float used to make the inelastic intensity decay with each
iteration

Returns
-------



Options
-------
The algorithm can return either the bulk or the film simulation

'''

    
class Algorithm1:
    def __init__(self, params):
        self.n = int(params['n'])
        self.P = params['P'] # primary input spectrum
        self.L = params['L'] # the loss function
        self.I = params['I'] # the inelastically scattered spectra
        self.elastic_prob = params['elastic_prob']
        self.inelastic_prob = params['inelastic_prob']
        self.inel_decay_factor = params['inel_decay_factor']
        self.el_decay_factor = params['el_decay_factor']
        self.option = params['option']
        self.option = params['option']

    def run(self):
        # This rescales the loss function by the total probability per elastic collision
        L = self.inelastic_prob * self.L 
        for i in range(self.n):
            # This condition is so that the first iteration uses the primary spectrum as input
            # All subsequent iterations use the convolved spectrum from the previous iteration
            if i == 0:
                self.I[-1] = self.P
            new = np.convolve(self.I[-1],np.flip(L))
            # Trim convolved spectrum so that it has the same lenth as the primary
            # spectrum
            l = len(self.P)
            new = new[-l:]
            # scale scattered spectrum by inelastic angle factor
            new = new * self.inel_decay_factor
            # Note: self.I has the shape of an (m,n) array, where m is the number
            # of spectra, and n is the number pf points per spectrum
            # therefore, new is made into a (1,n) array to match the shape
            # of self.I
            new = np.array([new])
            # factor to re-scale th input spectrum
            f = (1 - self.elastic_prob 
                 - self.inelastic_prob 
                 + self.elastic_prob 
                 * self.el_decay_factor) 
            # sum together the scattered and non-scattered spectra
            new = np.add(new, self.I[-1] * f)
            # append the new spectrum to the array of intermediate spectra
            self.I = np.concatenate((self.I,new),axis=0)
            
        if (self.option == 'film') | (len(self.option) == 0):
            simulated = self.I[-1]
            return simulated
        elif self.option == 'bulk':
            simulated = np.sum(self.I, axis=0)
            return simulated