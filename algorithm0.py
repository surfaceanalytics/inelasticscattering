# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 09:50:50 2020

@author: Mark
"""
from angular_spread_lorentzian import AngularSpreadCalc
import numpy as np

'''
This algorithm: 
1) iteratively convolves an input lineshape with a loss function
2) determines the number of iterations to use, based on the pressure
3) scales the spectra by an angular dispersion parameter
    
Parameters
----------
n: the number of scattering events to simulate
P: an array that represents the inpout spectrum
L: an array that represents the loss function
I: an array to hold the intermediate inelastically scattered spectra
    It has the shape (m, n), where m is the number of spectra and n is the number
    of data points per spectrum.
elastic_xsect: a float to represent the elastic scattering cross section of 
    the scatterer
inelastic_xsect: a float to represent the inelastic cross section
    density: a float to represent the density of the scatterer
distance: a float to represent the distance through the scattering medium the 
    electrons must travel
inel_decay_factor: a float used in the anglular spread algorithm
el_decay_factor: a float used in the angular spread algorithm
acceptance_angle: a float used in the angular spread algorithm
nr_iter_per_mfp: number of iterations to use per MFP

Returns
-------
An array that represents the simulated spectrum

Options
-------
The algorithm can return either the bulk or the film simulation

'''

class Algorithm0:
    def __init__(self, params):
        self.P = params['P'] # primary input spectrum
        self.L = params['L'] # the loss function
        self.I = params['I'] # the inelastically scattered spectra
        self.elastic_xsect = params['elastic_xsect']
        self.inelastic_xsect = params['inelastic_xsect']
        self.density = params['density']
        self.distance = params['distance']
        self.inel_decay_factor = params['inel_decay_factor']
        self.el_decay_factor = params['el_decay_factor']
        self.acceptance_angle = params['acceptance_angle']
        self.nr_iter_per_mfp = params['nr_iter_per_mfp']
        self.option = params['option']
        self._convertDist()
        self._calcMFP()
        self._calcProbs()
        
        
    def _convertDist(self):
        '''Distance is received in mm but is used in nm for calculations'''
        self.distance_nm = self.distance * 1000000

    def _calcMFP(self):
        if (self.elastic_xsect !=0) & (self.inelastic_xsect !=0):
            self.mfp = 1/(self.elastic_xsect * self.density) # the elastic mean free path in nm
            self.imfp = 1/(self.inelastic_xsect * self.density) # the inelastic mean free path in nm
            self.d_mfp = self.distance_nm / min(self.mfp,self.imfp) # this is the distance from sample to spectrometer in units of mfp
        elif self.elastic_xsect == 0:
            self.imfp = 1/(self.inelastic_xsect * self.density)    
            self.d_mfp = self.distance_nm / abs(self.imfp) # this is the distance from sample to spectrometer in units of mfp
        elif self.inelastic_xsect == 0:
            self.imfp = 1/(self.elastic_xsect * self.density)    
            self.d_mfp = self.distance_nm / abs(self.mfp)
        print('imfp: ' + str(self.imfp))
        print('d imfp: ' + str(self.d_mfp))

    def _calcProbs(self):
        if self.d_mfp < 1:
            self.n = 10
        else:
            # this is the number of iterations to use in the calculation
            self.n = int(self.d_mfp * self.nr_iter_per_mfp) 
            # this is the distance for one iteration in nm
        self.d_iter = self.distance_nm / self.n 
        self.inelastic_prob = 1 - np.exp(-1 * (self.distance_nm / self.n) * self.inelastic_xsect * self.density)
        self.elastic_prob = 1 - np.exp(-1 * (self.distance_nm / self.n) * self.elastic_xsect * self.density)
        
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