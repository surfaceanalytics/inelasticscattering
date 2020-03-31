# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:51:27 2020

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
inel_angle_factor: a float used in the anglular spread algorithm
el_angle_factor: a float used in the angular spread algorithm
acceptance_angle: a float used in the angular spread algorithm

Returns
-------

inel: an n-d array, representing one lineshape per n events. The data in the 
arrays can be interpreted as the spectra of the inelastically scattered spectra,
where n is the number of inelastic scattering events
el: an array representing the elastically scattered spectrum
non: an array representing the un-scattered spectrum
simulated: an array representing the scaled sum of inel, el and non. It can be 
interpreted as the measured spectrum

Options
-------
The algorithm can return either the bulk or the film simulation

'''

class Algorithm2:
    def __init__(self,params):
        self.n = params['n'] # the number of scattering events to calculate
        self.P = params['P'] # primary input spectrum
        self.L = params['L'] # the loss function
        self.I = params['I'] # the inelastically scattered spectra
        self.elastic_xsect = params['elastic_xsect']
        self.inelastic_xsect = params['inelastic_xsect']
        self.density = params['density']
        self.distance = params['distance']
        self.inel_angle_factor = params['inel_angle_factor']
        self.el_angle_factor = params['el_angle_factor']
        self.acceptance_angle = params['acceptance_angle']
        self.option = params['option']
        self._convertDist()
        
    def _convertDist(self):
        '''Distance is received in mm but is used in nm for calculations'''
        self.distance_nm = self.distance * 1000000

    def run(self):
        self._convolve()
        self._genPoisson() # generates vectors that contain Poisson dist for inelastic and elastic scattering
        self._makePoissonArray()
        self._genAngle()
        self._genAngleArray()
        self._genFactorArray()
        self._genFactors()
        if (self.option == 'film') | (len(self.option) == 0):
            return self._simulateFilm()
        elif self.option == 'bulk':
            return self._simulateBulk()

    def _convolve(self):
        for i in range(1,self.n):
            # This condition is so that the first iteration uses the primary spectrum as input
            # All subsequent iterations use the convolved spectrum from the previous iteration
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
        p = ((1/np.math.factorial(n)) 
            * ((distance * density * sigma)**n) 
            * np.exp(-1 * distance * density * sigma))
        return p
        

    def _genPoisson(self):
        # These are the Poisson distributions for inelastic and elastic scattering
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
        # Calculate the dot product of the two Poisson distributions
        # Rows represent number of times electron is inelastically scattered
        # Columns represent number of time an electron was elastically scattered
        self.T = np.dot(np.array([self.poiss_inel]).T, np.array([self.poiss_el]))
        self._genPFactors()
        
    def _genPFactors(self):
        # Calculate the scattering probabilities for elastic, inelastic and non-scattering
        self.p_non = self.T[0,0] # this is the probability of not being scattered over the distance d
        self.p_el = self.T[0,1:]
        self.p_inel = np.sum(self.T[1:,:], axis=1)

    def _calcAngleDist(self, width):
        angle_dist = AngularSpreadCalc(iterations=self.n,
                                acceptance_angle=self.acceptance_angle,
                                energy=500,
                                width=width)
        # Generate lorenztian
        angle_dist.gen_lorentzian_cross_section()
        # Run convolution
        angle_dist.run_convolution()
        # Limit by acceptance angle
        angle_dist.limit_by_acceptance_angle()
        # Calculate accepted area under curve
        area_sum = angle_dist.calc_area_under_curve()
        return area_sum
        
    def _genAngle(self):
        # Calculate the scaling factors due to angular spread for inelastic and elastic scattering
        self.angle_inel = self._calcAngleDist(width = self.inel_angle_factor)[0:-1]
        self.angle_el = self._calcAngleDist(width = self.el_angle_factor)[0:-1]

    def _genAngleArray(self):
        # An array of the angular spread factors for inelastic and elastic scattering
        self.A = np.dot(np.array([self.angle_inel]).T,np.array([self.angle_el]))

    def _genFactorArray(self):
        # Multipliy element wise the angle factors and the Poisson factors
        self.F = self.T * self.A
        
    def _genFactors(self):
        self.inel_factor = np.sum(self.F[1:,:], axis=1)
        self.el_factor = np.sum(self.F[0,1:])
        self.non_factor = self.F[0,0]

    def _simulateFilm(self):
        # For the case of scattering though film
        inel = np.array(np.dot(self.I[1:].T,np.array([self.inel_factor]).T))[:,0] 
        el = self.P * self.el_factor 
        non = self.P * self.non_factor 
        simulated = np.sum([np.array(inel),np.array(el), np.array(non)], axis=0)
        
        return inel, el, non, simulated
    
    def _simulateBulk(self):
        # For the case of scattering though bulk
        self.imfp = self.density * self.inelastic_xsect
        self.mfp = self.density * self.elastic_xsect
        inel_factor = self.angle_inel * self.imfp
        el_factor = np.sum(self.angle_el * self.mfp)
        inel = np.array(np.dot(self.I[1:].T,np.array([inel_factor]).T[1:]))[:,0] # scale all inelastic scattered spectra by their probabilities
        el = self.P * el_factor # scale the primary spectrum by the elastic probability for the elastically scattered signal
        non = self.P * (self.imfp + self.mfp) / 2 # scale the primary spectrum by the probability for no scattering 
        simulated = np.sum([np.array(inel),np.array(el), np.array(non)], axis=0)
        
        return inel, el, non, simulated
