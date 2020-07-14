# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 11:00:22 2020

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
        self.x = np.arange(self.start,self.stop,self.step)
        self.lineshape = np.zeros(len(self.x))

class Peak:
    def __init__(self, position, width, intensity):
        self.position = position
        self.width = width
        self.intensity = intensity
        
class Gauss(Peak):
    def __init__(self,position,width,intensity):
        Peak.__init__(self, position, width, intensity)

    def function(self, x):
        if self.width != 0:
            g = self.intensity / (self.width * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x-self.position)/self.width)**2)
            return g

class Lorentz(Peak):
    def __init__(self,position,width,intensity):
        Peak.__init__(self, position, width, intensity)    
    
    def function(self, x):
        if self.width != 0:
            l = self.intensity * 1 / (1 + ((self.position-x)/(self.width/2))**2)
            return l

class Voigt(Peak):
    def __init__(self, position, width, intensity, fraction_gauss):
        Peak.__init__(self, position, width, intensity)
        self.gauss = Gauss(position, width, intensity)
        self.lorentz = Lorentz(position, width, intensity)
        self.fraction_gauss = fraction_gauss
        
    def function(self,x):
        v = (self.fraction_gauss * self.gauss.function(x) 
        + (1-self.fraction_gauss) * self.lorentz.function(x))
        return v
        

class VacuumExcitation():
    def __init__(self, edge, fermi_width, intensity, exponent):
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
        if self.fermi_width !=0:
            f = (1-self.Fermi(x)) * self.Power(x) * self.intensity
            return f 
    
class Tougaard():
    def __init__(self,B, C, D, Eg):
        self.B = B
        self.C = C
        self.D = D
        self.Eg = Eg
        #self.t = 300 # Temperature in Kelvin
        #self.kb = 0.000086 # Boltzman constant
        
    def function(self, x):
        kb = 0.000086
        t = 300
        C = self.C * 20
        f = ((self.B * x) / ((C-x**2)**2 + self.D*x**2)
        * 1/(np.exp((self.Eg - x)/(t * kb)) + 1))
        return f
    
class SyntheticSpectrum(Spectrum):
    def __init__(self,start,stop,step):
        Spectrum.__init__(self,start,stop,step)
        self.components = []
        self.buildLine()
    
    def buildLine(self):
        self.clearLineshape()
        if len(self.components)==0:
            y = np.zeros(len(self.x))
            self.lineshape = y
        else:
            for component in self.components:
                y = np.array([component.function(x) for x in self.x])
                self.lineshape = np.add(self.lineshape,y)
            
    def normalize(self):
        if np.sum(self.lineshape) != 0:
            self.lineshape = self.lineshape / np.sum(self.lineshape)
            
    def addComponent(self,component, rebuild=True):
        self.components += [component]
        if rebuild == True:
            self.reBuild()   
        
    def removeComponent(self, comp_idx):
        del self.components[comp_idx]
        self.reBuild()

    def reBuild(self):
        self.updateRange()
        self.buildLine()
        
    def updateRange(self):
        self.x = np.arange(self.start,self.stop,self.step)
        
class SimulatedSpectrum(Spectrum):
    def __init__(self,start,stop,step):
        Spectrum.__init__(self,start,stop,step)
        self.convolved = np.array([]) # array to hold all the convolved spectrsa
        self.primary = np.array([]) # array to hold the non-scattered spectra
        self.poiss_inel = []
        self.poiss_el = []
        self.inel_angle = []
        self.el_angle = []
        self.inel_factor = []
        self.el_factor = []
        self.p_el = '' # total probability of elastic scattering over a distance d
        self.p_inel = '' # total probability of inelastic scattering over a distance d
        self.p_non = '' # total probability of not being scattered over a distance d

        self.bulk = {}
        self.film = {}
    
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
        self.lineshape = self.lineshape[:-1] - np.min(self.lineshape[:-1])
        self.x = self.x[:-1]

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
        
    def normalize(self):
        if np.sum(self.lineshape) != 0:
            self.lineshape = self.lineshape / (np.sum(self.lineshape)*self.step)
        
    def reBuild(self): # redefine the rebuild method for loss function (polymorphism)
        self.updateRange()
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1   
        
    def addVacuumExcitation(self, edge, fermi_width, intensity, exponent):
        self.components += [VacuumExcitation(edge, fermi_width, intensity, exponent)]
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1  
        
class Scatterer():
    def __init__(self):
        self.label = 'default'      
        self.loss_function = LossFunction(0,1200,0.1)
        self.cross_section = 0.01
        self.gas_diameter = 0.2 #In nanometers
        self.gas_cross_section = np.pi * (self.gas_diameter / 2)**2
        self.elastic_xsect = 0.01 # in units of nm^3
        self.inelastic_xsect = 0.01 # in units of nm^3
        self.inel_angle_factor = 10 
        self.el_angle_factor = 10
        self.inel_decay_factor = 1 
        self.el_decay_factor = 1
        self.inelastic_prob = 0.2
        self.elastic_prob = 0.2
        
class ScatteringMedium():
    def __init__(self):
        self.scatterer = Scatterer()
        self.R = 8.314463E+25 # gas constant in units of nm^3.mbar.K^-1.mol^-1
        self.avagadro = 6.022141E+23 # Avagadro's contant
        self.T = 300 # temperature in Kelvin
        self.pressure = 1 # In mbar
        self.density = self.pressure / (self.R * self.T) * self.avagadro # molecular density in units of particles per nm^3
        self.distance = 0.80 # In millimeters
        
    def calcDensity(self):
        self.density = self.pressure / (self.R * self.T) * self.avagadro # units are particles per nm^3
        
    def setPressure(self, pressure):
        self.pressure = pressure
        self.calcDensity()
        
    def setDistance(self,distance):
        self.distance = distance

class Calculation():
    def __init__(self):
        self.n_iter = 1
        self.nr_iter_per_mfp = 10
        self.n_events = 50