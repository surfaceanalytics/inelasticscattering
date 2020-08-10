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
    def __init__(self, position, width, intensity, fraction_gauss=0.5):
        Peak.__init__(self, position, width, intensity)
        self.fraction_gauss = fraction_gauss
        
    def function(self,x):
        if self.width != 0:
            v = ((self.fraction_gauss 
                 * Gauss(self.position,self.width,self.intensity).function(x))
                + ((1- self.fraction_gauss)
                * Lorentz(self.position,self.width,self.intensity).function(x)))
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
        self.buffer = LineShapeBuffer(self)
    
    def buildLine(self):
        """ This function removes the previous lineshape and re-constructs
        the lineshape from the present components in the componetns list.
        """
        self.clearLineshape()
        if len(self.components)==0:
            y = np.zeros(len(self.x))
            self.lineshape = y
        else:
            '''for component in self.components:
                y = np.array([component.function(x) for x in self.x])
                self.lineshape = np.add(self.lineshape,y)'''
            self.buffer._sum()
            
    def normalize(self):
        if np.sum(self.lineshape) != 0:
            self.lineshape = self.lineshape / np.sum(self.lineshape)
            
    def addComponent(self,component, rebuild=True):
        self.components += [component]
        self.buffer._addComponent()
        if rebuild == True:
            self.reBuild()
        
    def removeComponent(self, comp_idx):
        del self.components[comp_idx]
        self.buffer._delComponent(comp_idx)
        self.reBuild()
        
    def editComponent(self, comp_idx, params):
        """ params should be a dictionary where they keys are the attribute
        names and the values are the attribute values.
        """
        for i in params:
            try:
                params[i] = float(params[i])
            except:
                continue
            setattr(self.components[comp_idx], i, params[i])
        self.buffer._updateComponent(comp_idx) 
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
    def __init__(self, x, y):          
        self.start = x[0]
        self.step = x[1] - x[0]
        self.stop = x[-1]
        Spectrum.__init__(self, self.start, self.stop, self.step)
        self.x = x
        self.lineshape = y
    
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
        
class Scatterer():
    def __init__(self):
        self.label = 'default'      
        self.loss_function = LossFunction(0,1200,0.1)
        self.cross_section = 0.01
        self.gas_diameter = 0.2 #In nanometers
        self.gas_cross_section = np.pi * (self.gas_diameter / 2)**2
        self.inelastic_xsect = 0.01 # in units of nm^3
        self.norm_factor = 1 
        self.inelastic_prob = 0.2
        
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

class Test():
    def __init__(self):
        pass

class LineShapeBuffer():
    def __init__(self, synthetic_spectrum):
        """ The argument must be of type SyntheticSpectrum. It must have an
        attribute called components, consisting of a list of objects of type 
        Peak (such as Gauss, Lorentz or Voigt). The Peak object must have
        a method called function()."""
        self.synth_spec = synthetic_spectrum
        self._buildArray()
        
    def _buildArray(self):
        """ This function builds an array, where of shape (n,m) where n is the
        number of components in the synthetic spectrum and m is the number of
        data points per component (i.e. x values).
        """
        if np.sum(self.synth_spec.lineshape) == 0:
            self.array = np.array(self.synth_spec.lineshape)
        else:
            y = []
            for comp in self.synth_spec.components:
                y += [[comp.function(x) for x in self.synth_spec.x]]
            self.array = np.array(y)
        
    def _delComponent(self, comp_idx):
        self.array = np.delete(self.array, comp_idx, 0)
        self._sum()
        
    def _addComponent(self):
        new = np.array([self.synth_spec.components[-1].function(x) 
                        for x in self.synth_spec.x])
        if len(self.synth_spec.components) == 1:
            ''' This condition is run if this is the first component in
            the components list. This means the buffer array contains only 
            zeros. Then the zeros are replaced by the first component.
            '''
            self.array = new
        elif len(np.shape(self.array)) == 1:
            '''This condition is run when there is only one component. In this
            case, the shape of the buffer.array is (n,), where n is the
            number of points in the spectrum. In this case np.append works 
            differently than when the array has shape (n,m).'''
            self.array = np.append([self.array], [new], axis=0)
        else:
            '''This condition is run when there is more than one component.
            '''
            self.array = np.append(self.array, [new], axis=0)
        self._sum()
        
    def _updateComponent(self, comp_idx):
        new = np.array([self.synth_spec.components[comp_idx].function(x) 
                        for x in self.synth_spec.x])
        if len(np.shape(self.array)) == 1:
            self.array = new
        else:
            self.array[comp_idx,:] = new
    
    def _sum(self):
        if len(np.shape(self.array)) == 1:
            self.synth_spec.lineshape = self.array
        else:   
            self.synth_spec.lineshape = np.sum(self.array, axis=0)