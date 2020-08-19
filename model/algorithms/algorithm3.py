# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:07:22 2020

@author: Mark
"""

from angular_spread_lorentzian import AngularSpreadCalc
import numpy as np

'''
This algorithm uses: 
1) convolution to determine the inelastically scattered
lineshape, 
2) Poisson statistics to determine the scattering probabilities
3) an angular spread algorithm to determin the intensity loss factors
    
Parameters
----------
n: the number of scattering events to simulate
P: an array that represents the input spectrum
L: an array that represents the loss function
I: an array to hold the intermediate inelastically scattered spectra
    It has the shape (s, p), where s is the number of spectra and p is the number
    of data points per spectrum.
elastic_xsect: a float to represent the elastic scattering cross section of 
    the scatterer
inelastic_xsect: a float to represent the inelastic cross section
    density: a float to represent the density of the scatterer
distance: a float to represent the distance through the scattering medium the 
    electrons must travel

Returns
-------

inel: an n-d array, representing one lineshape per n events. The data in the 
arrays can be interpreted as the spectra of the inelastically scattered 
spectra, after all intensity scaling has been applied where n is the number of
inelastic scattering events
el: an array representing the elastically scattered spectrum
non: an array representing the un-scattered spectrum
simulated: an array representing the scaled sum of inel, el and non. It can be 
interpreted as the measured spectrum

Options
-------
The algorithm can return either the bulk or the film simulation

'''

class Algorithm3:
    algorithm_type = 'convolution'
    
    def __init__(self,params):
        self.n = params['n'] # the number of scattering events to calculate
        self.P = params['P'] # primary input spectrum
        self.L = params['L'] # the loss function
        self.I = params['I'] # the inelastically scattered spectra
        self.elastic_xsect = params['elastic_xsect']
        self.inelastic_xsect = params['inelastic_xsect']
        self.density = params['density']
        self.distance = params['distance']
        self.option = params['option']
        self._convertDist()
        
    def _convertDist(self):
        '''Distance is received in mm but is needed in nm for calculations'''
        self.distance_nm = self.distance * 1000000

    def run(self):
        self._convolve()
        self._genPoisson() # generates vectors that contain Poisson dist for inelastic and elastic scattering
        self._makePoissonArray()
        self._genFactors()
        if (self.option == 'film') | (len(self.option) == 0):
            return self._simulateFilm()
        elif self.option == 'bulk':
            return self._simulateBulk()

    def _convolve(self):
        ''' This convolves the input spectrum n times with the loss function.
        it returns an array of shape (n,p), where n is the number of iterations
        minus 1, and p is the number of data points per spectrum
        '''
        for i in range(1,self.n):
            '''This condition is so that the first iteration uses the primary 
            spectrum as input all subsequent iterations use the convolved 
            spectrum from the previous iteration
            '''
            if i == 1:
                new = np.convolve(self.P,np.flip(self.L))
            else:
                new = np.convolve(self.I[-1],np.flip(self.L))
            l = len(self.P)
            # Trim the unneeded data in the convolved spectrum
            new = new[-l:]  
            new = np.array([new])
            self.I = np.concatenate((self.I,new),axis=0)
                
    def Poisson(self, n, distance, sigma, density):
        ''' This function generates a point in a Poisson distribution'''
        p = ((1/np.math.factorial(n)) 
            * ((distance * density * sigma)**n) 
            * np.exp(-1 * distance * density * sigma))
        return p
        
    def _genPoisson(self):
        ''' This function generates the Poisson distributions for inelastic 
        and elastic scattering events. It returns an two arras of length n, 
        where n is the number of scattering events.
        '''
        self.poiss_inel = np.array([self.Poisson(i,
                                        self.distance_nm,
                                        self.inelastic_xsect,
                                        self.density) 
                                        for i in range(self.n)])
        self.poiss_el = np.array([self.Poisson(i,
                                          self.distance_nm,
                                          self.elastic_xsect,
                                          self.density) 
                                            for i in range(self.n)])
     
    def _makePoissonArray(self):
        '''Calculate the dot product of the two Poisson distributions
        Rows represent number of times electron is inelastically scattered
        Columns represent number of times an electron was elastically scattered
        '''
        self.T = np.dot(np.array([self.poiss_inel]).T, np.array([self.poiss_el]))
        self._genPFactors()
        
    def _genPFactors(self):
        '''This function calculates the scattering probabilities for elastic
        (p_el), inelastic (p_inel), and non-scattering (p_non).
        p_non is a float, p_el is a 1D array of length n, and p_inel is a 2D
        array of shape (n-1, n)
        '''
        self.p_non = self.T[0,0] # this is the probability of not being scattered over the distance d
        self.p_el = self.T[0,1:]
        self.p_inel = np.sum(self.T[1:,:], axis=1)
       
    def _genFactors(self):
        self.inel_factor = np.sum(self.T[1:,:], axis=1)
        self.el_factor = np.sum(self.T[0,1:])
        self.non_factor = self.T[0,0]

    def _simulateFilm(self):
        ''' This function is for the case of scattering though a film.
        first the convolved spectra are dotted with the factors vector to
        scale all the inelastic scatterd spectra by the Poisson and angular 
        factors. Then the elastically scattered and non-scattered spectra
        are scaled by their respective factors.
        '''
        self.inel = np.array(np.dot(self.I[1:].T,
                               np.array([self.inel_factor]).T))[:,0] 
        self.el = self.P * self.el_factor 
        self.non = self.P * self.non_factor 
        self.simulated = np.sum([np.array(self.inel),np.array(self.el), 
                                 np.array(self.non)], axis=0)
        
        return self.inel, self.el, self.non, self.simulated
    
    def _simulateBulk(self):
        # For the case of scattering though bulk
        self.imfp = 1/(self.density * self.inelastic_xsect)
        self.mfp = 1/(self.density * self.elastic_xsect)
        inel_factor = self.angle_inel * self.imfp #np.ones(self.n) * self.imfp
        el_factor = 1 #np.sum(self.angle_el * self.mfp)
        self.inel = np.array(np.dot(self.I[1:].T,np.array([inel_factor]).T[1:]))[:,0] # scale all inelastic scattered spectra by their probabilities
        self.el = self.P * el_factor # scale the primary spectrum by the elastic probability for the elastically scattered signal
        self.non = self.P * self.imfp # scale the primary spectrum by the probability for no scattering 
        self.simulated = np.sum([np.array(self.inel),np.array(self.el), 
                                 np.array(self.non)], axis=0)
        
        return self.inel, self.el, self.non, self.simulated
