# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np

class Spectrum:
    def __init__(self,start,stop,step):
        self.start = start
        self.stop = stop
        self.step = step
        self.x = np.arange(self.start,self.stop,self.step)
        self.clearLineshape()
        self.visibility = 'visible'
        self.kind = 'none'
        
    def clearLineshape(self):
        self.lineshape = np.zeros(int((self.stop-self.start)/self.step))
        self.x = np.arange(self.start,self.stop,self.step)

class Peak:
    def __init__(self, mean, stdev, intensity):
        self.mean = mean
        self.stdev = stdev
        self.intensity = intensity

    def Gauss(self, x):
        g = self.intensity / (self.stdev * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x-self.mean)/self.stdev)**2)
        return g
    
    def Lorentzian(self,x):
        l = self.intensity * 1 / (1 + ((self.mean-x)/(self.stdev/2))**2)
        return l
        
    def function(self, x):
        f = self.Lorentzian(x)
        return f
    
class SyntheticSpectrum(Spectrum):
    def __init__(self,start,stop,step):
        Spectrum.__init__(self,start,stop,step)
        self.components = []
    
    def buildLine(self):
        self.clearLineshape()
        for component in self.components:
            y = np.array([component.function(x) for x in self.x])
            self.lineshape = np.add(self.lineshape,y)
            
    def normalize(self):
        self.lineshape = self.lineshape / np.sum(self.lineshape) 
    
class MeasuredSpectrum(Spectrum):
    def __init__(self, filename):
        self.filename = filename
        self.data = self.convert(self.filename)
        x = self.data[:,0]
        x1 = np.roll(x,-1)
        diff = np.abs(np.subtract(x,x1))
        self.step = round(np.min(diff[diff!=0]),2)
        x = x[diff !=0]
        self.start = np.min(x)
        self.stop = np.max(x)
        Spectrum.__init__(self, self.start,self.stop,self.step)
        self.lineshape = self.data[:,2][diff != 0]
        self.x = x
        if (self.stop-self.start)/self.step > len(self.x):
            self.interpolate()
        self.lineshape = self.lineshape - np.min(self.lineshape)

    def convert(self, filename):
        file = open(filename,'r')
        lines = []
        for line in file.readlines():
            lines += [line]
        lines = lines[4:]
        lines = [[float(i) for i in line.split()] for line in lines]
        data = np.array(lines)
        return data
    
    def interpolate(self):
        new_x = []
        new_y = []
        for i in range(len(self.x)-1):
            diff = np.abs(np.around(self.x[i+1]-self.x[i],2))
            if (diff > self.step) & (diff < 10):
                for j in range(int(np.round(diff/self.step))):
                    new_x += [self.x[i] + j*self.step]
                    k = j / int(diff/self.step)
                    new_y += [self.lineshape[i]*(1-k) + self.lineshape[i+1]*k]
            else:
                new_x += [self.x[i]]
                new_y += [self.lineshape[i]]
                
        new_x += [self.x[-1]]
        new_y += [self.lineshape[-1]]
        self.x = new_x
        self.lineshape = np.array(new_y)
    
class LossFunction(SyntheticSpectrum):
    def __init__(self,start,stop,step):
        SyntheticSpectrum.__init__(self,start,stop,step)
        
    def addPeak(self,mean,stdev,intensity):
        self.components += [Peak(mean,stdev,intensity)]
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1     
            
    def addVacuumExcitation(self, power_exponent, fermi_edge, fermi_width, intensity):
        self.components += [VacuumExcitation(power_exponent, fermi_edge, fermi_width, intensity)]
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1  
        
    def reBuild(self):
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1  

class VacuumExcitation():
    def __init__(self, exponent, edge, fermi_width, intensity):

        self.edge = edge
        self.fermi_width = fermi_width
        self.intensity = intensity
        self.exponent = exponent
        
    def Fermi(self, x):
        k = 0.1
        f = 1/(np.exp((x-self.edge)/(k*self.fermi_width))+1)
        return f
    
    def Power(self, x):
        p = np.exp(-1 * (x+self.edge) * self.exponent)
        return p
    
    def function(self,x):
        f = (1-self.Fermi(x)) * self.Power(x) * self.intensity
        return f 

class Scatterer():
    def __init__(self):
        self.label = 'default'
        self.loss_function = LossFunction(0,200,0.1)
        self.cross_sec = 0.01
        self.gas_diameter = 0.2 #In nanometers
        self.gas_cross_section = np.pi * (self.gas_diameter / 2)**2
        
    def setCrossSec(self, cross_sec):
        self.cross_sec = cross_sec
        self.loss_function.buildLine()
        self.loss_function.normalize()
        
class ScatteringMedium():
    def __init__(self):
        self.scatterer = Scatterer()
        self.d_through_gas = 800000 # In nanometers
        self.pressure = 1 # In mbar
        self.const = 9323 # calculated from the gas constant at T = 300 K
        self.n_iter = int(100)
        self.calcParams()
        
    def calcParams(self):
        self.mean_free_path = self.const/(self.scatterer.gas_diameter**2 * self.pressure)
        self.d_mfp = self.d_through_gas / self.mean_free_path # this is the distance in units of mfp
        self.elast_prob = 0.5 * (self.d_through_gas / self.n_iter) / self.mean_free_path # This is the elastic scattering probability per iteration (0.5 is the elast_scatter prob for 1 * MFP)
        
    def setPressure(self, pressure):
        self.pressure = pressure
        self.calcParams()
        
    def setDistance(self,distance):
        self.d_through_gas = distance
        self.calcParams()    