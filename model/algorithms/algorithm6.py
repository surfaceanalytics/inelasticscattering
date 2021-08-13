# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:55:44 2020

@author: Mark
"""
import numpy as np

class Algorithm6:
    """
    This it eh 'Unscatter' algorithm.
    This algorithm deconvolves spectra. It takes a lineshape that represents
    the convolution of two spectra. Then it takes the known 'loss function'
    spectrum and deconvolves it from the input spectrum.
        
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
    norm_factor: a float used in the anglular spread algorithm
    
    Returns
    -------
    
    inel: An n-d array, representing one lineshape per n events. The data in 
    the arrays can be interpreted as the spectra of the inelastically 
    scattered spectra, where n is the number of inelastic scattering events.
    
    simulated: An array representing the decovolved spectrum.
    
    Options
    -------
    The algorithm can return either the bulk or the film simulation
    """
    algorithm_type = 'deconvolution'
            
    def __init__(self, inputSpec, scattering_medium, params):
        """
        Parameters
        ----------
        inputSpec : A Spectrum object, as defined in the base_model module.
            It must have an attribute called lineshape, which is a numpy array
            representing the y intensities of an input spectrum.
        scattering_medium : ScatteringMedium, as defined in the base_model
            module. It must have attributes for distance, density, a method
            called calcDensity, as well as an attribute called scatterer, 
            which is a Scatterer object as defined in the base_model mudule. 
            The Scatterer must have attributes inelastic_xsect [FLOAT], 
            norm_factor [FLOAT] and lineshape [NUMPY ARRAY].
            DESCRIPTION.
        params : DICTIONARY

        Returns
        -------
        None.

        """

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
        return self._simulateFilm()
        
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
             
    def _simulateFilm(self):
        """ This function is for the case of scattering though a film.
        first the convolved spectra are dotted with the factors vector to
        scale all the inelastic scatterd spectra by the Poisson and angular 
        factors. Then the elastically scattered and non-scattered spectra
        are scaled by their respective factors.
        """
        
        self._removeMin()
        '''Calculate the Fourier transform of the loss function.'''
        fft_loss = np.fft.fft(np.flip(self.L))
        
        poisson_factor = self.distance_nm * self.inelastic_xsect * self.density
        norm_factor = self.norm_factor
        total_factor = poisson_factor * norm_factor
        fft_multi_loss = np.exp(total_factor*fft_loss)
        
        ''' Pad the input spectrum with zeros so that it has the same 
        dimensions as the loss function.'''
        input_spec_padded = np.pad(self.P, (len(self.L)-len(self.P),0), 
                                'constant', constant_values = 0)
        fft_input_spec = np.fft.fft(input_spec_padded)
        exp_factor = np.exp(-1 * poisson_factor)
        
        fft_deconv = np.divide(fft_input_spec, fft_multi_loss) / exp_factor
        deconv = np.real(np.fft.ifft(fft_deconv))
        deconv = deconv[-len(self.P):]
        
        self.I = deconv
        
        self._addMin()
     
        self.inel = np.subtract(self.P, deconv)
        self.simulated = deconv + self.min_value  
        return self.simulated