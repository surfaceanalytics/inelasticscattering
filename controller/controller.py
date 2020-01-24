# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:21:27 2020

@author: Mark
"""
from gasscattering.model.model_classes import Spectrum, ScatteringMedium, MeasuredSpectrum
from  gasscattering.data.scatterers import read_scatterers
import numpy as np
import pandas as pd
import os

datapath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + 'data\\'

class Controller():
    def __init__(self):
        self.start = 1200
        self.stop = 1400
        self.step = 0.1
        
        self.input_spectrum = Spectrum(self.start,self.stop,self.step)
        self.spectrum_to_fit = Spectrum(self.start,self.stop,self.step)
        self.output_spectrum = Spectrum(self.start,self.stop,self.step)
        self.scattering_medium = ScatteringMedium()
        self.loadScatterer('default')
        
        self.intermediate_spectra =[]
        self.scattered_spectra = []
        self.bulk_spectrum = []
        
    def loadScatterer(self, label):
        self.scattering_medium.scatterer = read_scatterers(label)
        self.scattering_medium.scatterer.loss_function.step = self.step
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()

    def loadPrimarySpectrum(self, filename):
        self.input_spectrum = MeasuredSpectrum(datapath + filename)
        self.start = self.input_spectrum.start    # The step width must be defined by the measured spectrum 
        self.stop = self.input_spectrum.stop      # All synthetic pectra need to have their step widths redefined
        self.step = self.input_spectrum.step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        #self.scattering_medium.scatterer.loss_function.clearLineshape()
        #self.scattering_medium.scatterer.loss_function.buildLine() # Re-build loss function lineshape
        #self.scattering_medium.scatterer.loss_function.normalize()
        self.output_spectrum = Spectrum(self.start,self.stop,self.step) # Overwrite old output spectrum with new settings
        
    def loadSpectrumToFit(self, filename):
        self.spectrum_to_fit = MeasuredSpectrum(datapath + filename)
        
    def scatterSpectrum(self):
        n = self.scattering_medium.n_iter
        a = self.input_spectrum.lineshape
        p = self.scattering_medium.elast_prob * self.scattering_medium.scatterer.cross_sec
        b = p * self.scattering_medium.scatterer.loss_function.lineshape
        self.intermediate_spectra = [a]
        for i in range(n):
            c = np.convolve(a,np.flip(b))
            l = len(a)
            c = c[-l:]
            a = (1-p) * a
            a = np.add(c,a)
            self.intermediate_spectra += [a]
        self.output_spectrum.lineshape = a
        self.intermediate_spectra = np.stack(self.intermediate_spectra)
        self.bulk_spectrum = np.sum(self.intermediate_spectra, axis=0)
        
    def exportExcel(self, filename):
        input_x = self.input_spectrum.x
        input_spec = self.input_spectrum.lineshape / np.max(self.input_spectrum.lineshape)
        output_spec = self.output_spectrum.lineshape / np.max(self.output_spectrum.lineshape)
        fit = self.spectrum_to_fit.lineshape / np.max(self.spectrum_to_fit.lineshape)
        df1 = pd.DataFrame(np.stack([input_x,input_spec,output_spec,fit]).T, columns=['kinetic energy (eV)','unscattered spectrum', 'fit spectrum','scattered spectrum'])
        loss_x = self.scattering_medium.scatterer.loss_function.x
        loss = self.scattering_medium.scatterer.loss_function.lineshape
        df2 =  pd.DataFrame(np.stack([loss_x,loss]).T,columns=['energy loss (eV)','proability'])
        file = filename + '.xlsx'
        with pd.ExcelWriter(file) as writer:
            df1.to_excel(writer,sheet_name='spectra',startcol=0)
            df2.to_excel(writer,sheet_name='spectra',startcol=5)
