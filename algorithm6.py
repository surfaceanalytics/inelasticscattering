# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:55:44 2020

@author: Mark
"""
import numpy as np



class Algorithm6:
    """
    This algorithm uses: 
    1) convolution to determine the inelastically scattered
    lineshape, 
    2) Poisson statistics to determine the scattering probabilities
    3) an angular spread algorithm to determine the intensity loss factors
        
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
    norm_factor: a float used in the anglular spread algorithm
    acceptance_angle: a float used in the angular spread algorithm
    
    Returns
    -------
    
    inel: an n-d array, representing one lineshape per n events. The data in 
    the arrays can be interpreted as the spectra of the inelastically 
    scattered spectra, where n is the number of inelastic scattering events
    el: an array representing the elastically scattered spectrum
    non: an array representing the un-scattered spectrum
    simulated: an array representing the scaled sum of inel, el and non. 
    It can be interpreted as the measured spectrum.
    
    Options
    -------
    The algorithm can return either the bulk or the film simulation
    """
        
    def __init__(self, inputSpec, scattering_medium, params):
        self.inputSpec = inputSpec
        self.scattering_medium = scattering_medium
        self.P = self.inputSpec.lineshape
        self.L = (self.scattering_medium.scatterer.loss_function.lineshape * 
              self.scattering_medium.scatterer.loss_function.step) 
        self.I = np.array([np.zeros(len(self.L))]) 
        self.inelastic_xsect = self.scattering_medium.scatterer.inelastic_xsect
        self.scattering_medium.calcDensity()
        self.density = self.scattering_medium.density
        self.distance = self.scattering_medium.distance
        self.norm_factor = self.scattering_medium.scatterer.norm_factor
        self.n = params['n_events'] 
        self.option = params['option']
        self._convertDist()
        
        self.user_inputs = [{'label':'Inel.\nX-sect.',
                             'variable':self.inelastic_xsect},
                            {'label':'Norm.\nfactor',
                             'variable':self.norm_factor},]
        
    def _convertDist(self):
        '''Distance is received in mm but is needed in nm for calculations'''
        self.distance_nm = self.distance * 1000000

    def run(self):
        self._genPoisson()
        self._genNormFactors()
        self._genScaleFactors()
        
        if (self.option == 'film') | (len(self.option) == 0):
            return self._simulateFilm()
        elif self.option == 'bulk':
            return self._simulateBulk()
        
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
        
                
    def Poisson(self, n, distance, sigma, density):
        """ This function generates a point in a Poisson distribution"""
        p = ((1/np.math.factorial(n)) 
            * ((distance * density * sigma)**n) 
            * np.exp(-1 * distance * density * sigma))
        return p
        
    def _genPoisson(self):
        """ This function generates the Poisson distributions for inelastic 
        scattering events. It returns an n-by-1 array, where n is the 
        number of scattering events.
        """
        self.poiss_inel = np.array([self.Poisson(i,
                                        self.distance_nm,
                                        self.inelastic_xsect,
                                        self.density) 
                                        for i in range(self.n+1)])
 
    def _genNormFactors(self):
        """ This generates an array of normalization factors. The array
        has shape (n,), where n is the number of scattering events.
        """
        self.norm_factors = np.array([self.norm_factor ** (i+1) 
                                        for i in range(self.n+1)])
        self.norm_factors[0] = 1

    def _genScaleFactors(self):
        """ This multiplies element-wise the Poisson factors with the 
        normalization factors. It results in an (n,) array, where n is the
        number of scattering events.
        """
        self.scale_factors = np.multiply(self.poiss_inel, self.norm_factors)

    def _simulateFilm(self):
        """ This function is for the case of scattering though a film.
        first the convolved spectra are dotted with the factors vector to
        scale all the inelastic scatterd spectra by the Poisson and angular 
        factors. Then the elastically scattered and non-scattered spectra
        are scaled by their respective factors.
        """
        
        self._removeMin()
        fft_loss = np.fft.fft(np.flip(self.L))
        for i in range(0,self.n+1):
            new = np.power(fft_loss,i) * self.scale_factors[i]
            new = np.array([new])
            self.I = np.concatenate((self.I,new),axis=0)
        self.I = np.sum(self.I, axis=0)

        
        input_spec_padded = np.pad(self.P, (len(self.L)-len(self.P),0), 
                                'constant', constant_values = 0)
        fft_input_spec = np.fft.fft(input_spec_padded)
        
        fft_deconv = np.divide(fft_input_spec, self.I)
        deconv = np.real(np.fft.ifft(fft_deconv))
        deconv = deconv[-len(self.P):]
        
        self.I = np.real(np.fft.ifft(self.I))
        
        self._addMin()
     
        self.inel = np.subtract(self.P, deconv)
        self.simulated = deconv + self.min_value  
        return self.simulated
    
    def _simulateBulk(self):
        """ This function is run in the case of simulating scattering though 
        bulk. It does not require the Poisson factors to be calculated, and
        therefore does not require the pressure or distance parameters.
        """
        self.inel = np.multiply(self.I.T, self.norm_factors)
        self.simulated = np.sum(self.inel, axis=1) + self.min_value       
        return self.simulated
