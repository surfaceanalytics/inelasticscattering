# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 11:57:39 2020

@author: Mark
"""
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

norm_factor: a float used to make the inelastic intensity decay with each
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
    def __init__(self, inputSpec, scattering_medium, params):
        self.inputSpec = inputSpec
        self.scattering_medium = scattering_medium
        self.P = self.inputSpec.lineshape
        self.L = (self.scattering_medium.scatterer.loss_function.lineshape * 
              self.scattering_medium.scatterer.loss_function.step) 
        self.I = np.array([np.zeros(len(self.inputSpec.lineshape))])     
        self.n = int(params['n'])
        self.inelastic_prob = self.scattering_medium.scatterer.inelastic_prob
        self.norm_factor = self.scattering_medium.scatterer.norm_factor
        self.option = params['option']
        
    def _removeMin(self):
        """
        This function subtracts the minimum value from the XPS spectrum
        before it convolves it.

        Returns
        -------
        None.

        """
        self.min_value = np.min(self.P)
        self.P = self.P - self.min_value
        
    def _addMin(self):
        """
        This function adds back the minimum value after the convolution is
        complete.

        Returns
        -------
        None.

        """
        self.P = self.P + self.min_value

    def run(self):
        """This rescales the loss function by the total probability per elastic 
        collision.
        """
        L = self.inelastic_prob * self.L 
        self._removeMin()
        for i in range(self.n):
            ''' This condition is so that the first iteration uses the primary 
            spectrum as input. All subsequent iterations use the convolved 
            spectrum from the previous iteration.'''
            if i == 0:
                self.I[-1] = self.P
            new = np.convolve(self.I[-1],np.flip(L))
            ''' Trim convolved spectrum so that it has the same lenth as the 
            primary spectrum.'''
            l = len(self.P)
            new = new[-l:]
            ''' Scale scattered spectrum by the normalization factor.'''
            new = new * self.norm_factor
            ''' Note: self.I has the shape of an (m,n) array, where m is the 
            number of spectra, and n is the number pf points per spectrum
            therefore, new is made into a (1,n) array to match the shape
            of self.I'''
            new = np.array([new])
            '''Factor to re-scale th input spectrum'''
            f = (1 - self.inelastic_prob) 
            ''' Sum together the scattered and non-scattered spectra.'''
            new = np.add(new, self.I[-1] * f)
            ''' Append the new spectrum to the array of intermediate spectra.'''
            self.I = np.concatenate((self.I,new),axis=0)
            
        self._addMin()
            
        if (self.option == 'film') | (len(self.option) == 0):
            simulated = self.I[-1] + self.min_value
            return simulated
        elif self.option == 'bulk':
            simulated = np.sum(self.I, axis=0) + self.min_value
            return simulated