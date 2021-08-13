# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:55:44 2020

@author: Mark
"""
import numpy as np
from model.algorithms.pathlengthapproximator import PathlengthApproximator

class Algorithm5:
    """
    Gas-phase scattering algorithm
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
    inelastic_xsect: a float to represent the inelastic cross section
        density: a float to represent the density of the scatterer
    distance: a float to represent the distance through the scattering medium the 
        electrons must travel
    norm_factor: a float used to rescale the convolved spectra.
    
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
    algorithm_type = 'convolution'
    
    def __init__(self, inputSpec, scattering_medium, params):
        self.inputSpec = inputSpec
        self.scattering_medium = scattering_medium
        self.P = self.inputSpec.lineshape
        self.L = (self.scattering_medium.scatterer.loss_function.lineshape * 
              self.scattering_medium.scatterer.loss_function.step) 
        self.I = np.array([np.zeros(len(self.inputSpec.lineshape))]) 
        self.inelastic_xsect = self.scattering_medium.scatterer.inelastic_xsect
        self.scattering_medium.calcDensity()
        self.density = self.scattering_medium.density
        self.distance = self.scattering_medium.distance
        self.norm_factor = self.scattering_medium.scatterer.norm_factor
        if 'n_events' in params.keys():
            self.n = params['n_events']
        else:
            self.n = 50
        self.option = params['option']
        self._convertDist()

        self.user_inputs = [{'label':'Inel.\nX-sect.',
                             'variable':self.inelastic_xsect},
                            {'label':'Norm.\nfactor',
                             'variable':self.norm_factor},]
        
    def _convertDist(self):
        """
        Converts the distance (received in mm) to nm, which is needed for the
        calculations.

        Returns
        -------
        None.

        """
        self.distance_nm = self.distance * 1000000

    def run(self):
        """
        This function runs the simulation. It first prepares the needed 
        factors, then calls a simulation method, depending on the kind of 
        simulation that should be run.    

        Returns
        -------
        Returns a numpy array, representing the simulated spectrum.

        """
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
        """
        This function generates a point in a Poisson distribution.

        Parameters
        ----------
        n : INT
            The Poisson random variable for a Poisson point proccess.
            In the case of electron scattering, it represents the number of
            times the electron was scattered.
        distance : FLOAT or INT
            Used to calculate lambda in the Poisson point process.
            In the context of electron scattering, it represents the distance
            through the scattering medium that the electron travels.
        sigma : FLOAT
            Used to calculate lambda in the Poisson point process.
            In the context of electron scattering, it represents the 
            scattering cross section of the scatterer.
        density : FLOAT
            Used to calculate lambda in the Poisson point process.
            It represents the density of the scatterer.

        Returns
        -------
        p : FLOAT
            The Poisson probability.
            In the context of electron scattering, it represents the 
            probability that an electron is scattered n times when travelling
            a distance d through the scattering medium, having a scatterer
            with cross section and density.

        """
        p = ((1/np.math.factorial(n)) 
            * ((distance * density * sigma)**n) 
            * np.exp(-1 * distance * density * sigma))
        #p = distance * density * sigma
        return p

    def _simulateFilm(self):
        """ This function is for the case of scattering though a film.
        first the convolved spectra are dotted with the factors vector to
        scale all the inelastic scatterd spectra by the Poisson and angular 
        factors. Then the elastically scattered and non-scattered spectra
        are scaled by their respective factors.
        """
        self._removeMin()
        loss = np.flip(self.L)
        ''' Pad the imput spectrum with zeros so that it has the same
        dimensions as the loss function.
        '''
        input_spec_padded = np.pad(self.P, (len(loss)-len(self.P),0), 
                                'constant', constant_values = 0)
        ''' Take Fourier transform of the input spectrum.
        '''
        fft_input_spec = np.fft.fft(input_spec_padded)
        
        ''' Approximate the average path length'''
        d = self._getPathLength(self.distance_nm)
        ''' Take the Fourier transform of the loss function.
        '''
        fft_loss = np.fft.fft(loss)
        poisson_factor = d * self.inelastic_xsect * self.density
        norm_factor = self.norm_factor
        total_factor = poisson_factor * norm_factor
        
        exp_factor = np.exp(-1 * poisson_factor)
        
        #print('exp_factor: ' + str(exp_factor))
        #print('total_factor: ' + str(total_factor))
        #print('sum fft input: ' + str(np.sum(fft_input_spec)))
        #print('sum fft input: ' + str(np.sum(fft_loss)))
              
        fft_total = exp_factor * np.multiply(fft_input_spec, 
                                             np.exp(total_factor*fft_loss))
        
        ''' Take the inverse Fourier transform of the convolved spectrum.
        '''
        total = np.real(np.fft.ifft(fft_total)[-len(self.P):])
        
        self.inel = total - (self.P * np.exp(-self.distance_nm / self._getLambda()))
        self.simulated = total + self.min_value  
        self._addMin()

        return self.simulated
    
    def _getLambda(self):
        return 1/(self.inelastic_xsect * self.density)
    
    def _simulateBulk(self):
        """ This function is run in the case of simulating scattering though 
        bulk. It does not require the Poisson factors to be calculated, and
        therefore does not require the pressure or distance parameters.
        """
        self._removeMin() 
        n_events = 50 #self.n_events
        
        loss = np.flip(self.L)

        input_spec_padded = np.pad(self.P, (len(loss)-len(self.P),0), 
                                'constant', constant_values = 0)

        fft_input_spec = np.fft.fft(input_spec_padded)
        fft_loss = np.fft.fft(loss)

        fft_total = np.zeros((n_events, len(fft_input_spec)))
        
        for i in range(n_events):
            fft_total[i] = np.multiply(fft_input_spec, np.power(fft_loss,i)).copy()
        
        fft_total = np.sum(fft_total, axis=0)

        
        total = np.real(np.fft.ifft(fft_total)[-len(self.P):]) * self._getLambda()
        
        self.inel = total - self.P
        self.simulated = total + self.min_value  
        self._addMin()
        return self.simulated
    
    def _getFactor(self):
        d = self._getPathLength(self.distance_nm)
        poisson_factor = d * self.inelastic_xsect * self.density
        norm_factor = self.norm_factor
        total_factor = poisson_factor * norm_factor
        
        return total_factor
    
    def _getPathLength(self,x):
        """
        This function uses a fit function to estimate the average path-length
        electrons travel on their way from the sample to the nozzle.

        Parameters
        ----------
        x : FLOAT
           The distance between the sample and the nozzle.

        Returns
        -------
        None.

        """
        pathlength = PathlengthApproximator()
        self._checkValidity()
        p = pathlength.function(x)
        return p
    
    def _checkValidity(self):
        """
        This function checks if the approximation used for calculating path-
        length is valid under the used parameters.

        Returns
        -------
        None.

        """
        cutoff_factor = 6
        imfp = self._getLambda()
        print(imfp)
        d = self.distance_nm
        if (d/imfp > cutoff_factor):
            self.warning ='''The approximation made for pathlength may not be valid under these conditions,
because the sample-nozzle distance is more than {} times greater than the sample-nozzle distance.'''.format(cutoff_factor)