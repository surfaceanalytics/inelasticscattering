# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json

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
        g = self.intensity / (self.width * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x-self.position)/self.width)**2)
        return g

class Lorentz(Peak):
    def __init__(self,position,width,intensity):
        Peak.__init__(self, position, width, intensity)    
    
    def function(self, x):
        l = self.intensity * 1 / (1 + ((self.position-x)/(self.width/2))**2)
        return l
    
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
        f = (1-self.Fermi(x)) * self.Power(x) * self.intensity
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
            
    def addComponent(self,component):
        self.components += [component]
        self.reBuild()   
        
    def removeComponent(self, comp_idx):
        del self.components[comp_idx]
        self.reBuild()

    def reBuild(self):
        self.updateRange()
        self.buildLine()
        
    def updateRange(self):
        self.x = np.arange(self.start,self.stop,self.step)
    
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
        self.loss_function = LossFunction(0,200,0.1)
        self.cross_section = 0.01
        self.gas_diameter = 0.2 #In nanometers
        self.gas_cross_section = np.pi * (self.gas_diameter / 2)**2
        self.angle_factor = 1
        
    def setCrossSec(self, cross_section):
        self.cross_section= cross_section
        self.loss_function.buildLine()
        self.loss_function.normalize()
        
class ScatteringMedium():
    def __init__(self):
        self.scatterer = Scatterer()
        self.d_through_gas = 800000 # In nanometers
        self.pressure = 1 # In mbar
        self.const = 9323 # calculated from the gas constant at T = 300 K
        self.n_iter = int(100)
        self.nr_iter_per_mfp = 10

        self.calcParams()
        
    def calcParams(self):
        self.mean_free_path = self.const/(self.scatterer.gas_diameter**2 * self.pressure)
        self.d_mfp = self.d_through_gas / self.mean_free_path # this is the distance in units of mfp
        self.n_iter = int(self.d_mfp * self.nr_iter_per_mfp)
        self.d_iter = self.d_through_gas / self.n_iter
        #if d_iter > self.mean_free_path:
        self.collis_prob = 0.5 * (self.d_through_gas / self.n_iter) / self.mean_free_path # This is the elastic scattering probability per iteration (0.5 is the elast_scatter prob for 1 * MFP)
           
    def setPressure(self, pressure):
        self.pressure = pressure
        self.calcParams()
        
    def setDistance(self,distance):
        self.d_through_gas = distance
        self.calcParams()    
              
class Model():
    def __init__(self, controller):
        self.controller = controller
        self.start = 1200
        self.stop = 1400
        self.step = 0.1
        self.loaded_spectra = []
        self.unscattered_spectrum = Spectrum(self.start,self.stop,self.step)
        self.scattered_spectrum = Spectrum(self.start,self.stop,self.step)
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step)
        self.scattering_medium = ScatteringMedium()
        self.intermediate_spectra =[]
        self.scattered_spectra = []
        self.bulk_spectrum = []
        self.readScatterers()
        self.loss_component_kinds = ['Gauss', 'Lorentz', 'VacuumExcitation']
        self.peak_kinds = ['Gauss', 'Lorentz']

    def loadSpectrum(self, filename):
        self.loaded_spectra += [MeasuredSpectrum(filename)]
        self.start = self.loaded_spectra[-1].start    # The step width must be defined by the measured spectrum 
        self.stop = self.loaded_spectra[-1].stop      # All synthetic pectra need to have their step widths redefined
        self.step = self.loaded_spectra[-1].step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step) # Overwrite old output spectrum with new settings

    def scatterSpectrum(self):
        n = self.scattering_medium.n_iter
        a = self.unscattered_spectrum.lineshape # this is the initial input spectrum
        p_coll = self.scattering_medium.collis_prob # this is the collision probability per unit distance
        p_elast = self.scattering_medium.collis_prob * (1 - self.scattering_medium.scatterer.cross_section) # this is the probability of elastic scattering per unit distance
        p_inelast= self.scattering_medium.collis_prob * self.scattering_medium.scatterer.cross_section # this gives the total probability of inelastic scattering per unit distance
        #print("elatic probability: " + str(self.scattering_medium.collis_prob))
        #print("d_mfp: " + str(self.scattering_medium.d_mfp))
        #print("mfp: " + str(self.scattering_medium.mean_free_path))
        #print("p: " + str(p))
        #print("d_iter: " + str(self.scattering_medium.d_iter))
        self.scattering_medium.scatterer.step = self.unscattered_spectrum.step
        self.scattering_medium.scatterer.loss_function.reBuild()
        b = self.scattering_medium.scatterer.angle_factor * p_inelast* self.scattering_medium.scatterer.loss_function.lineshape # This rescales the loss function by the total probability per elastic collision
        self.intermediate_spectra = [a]
        for i in range(n):
            c = np.convolve(a,np.flip(b)) # this convolves the input spectrum with the scaled loss function
            l = len(a)
            c = c[-l:] # this trims the unneeded data in the convolved spectrum
            a = (1-p_coll + p_elast * self.scattering_medium.scatterer.angle_factor) * a # this rescales the non-scattered portion of the spectrum
            a = np.add(c,a) # this sums the non-scattered and scattered spectra. it will be the input spectrum for next iteration
            self.intermediate_spectra += [a]
        self.simulated_spectrum.lineshape = a
        self.simulated_spectrum.x = self.unscattered_spectrum.x 
        self.simulated_spectrum.kind = 'Simulated'
        self.bulk_spectrum = np.sum(self.intermediate_spectra, axis=0)

    def setCurrentScatterer(self, label):
        self.scattering_medium.scatterer.label = label
        self.scattering_medium.scatterer.cross_section = self.scatterers[label]['cross_section']
        self.scattering_medium.scatterer.gas_diameter = self.scatterers[label]['gas_diameter']
        self.scattering_medium.scatterer.angle_factor = self.scatterers[label]['angle_factor']
        self.scattering_medium.scatterer.loss_function.components = []
        for i in self.scatterers[label]['loss_function']:
            if i['type'] == 'Gauss':
                self.scattering_medium.scatterer.loss_function.addComponent(
                        Gauss(i['params']['position'], i['params']['width'], 
                        i['params']['intensity']))
            elif i['type'] == 'Lorentz':
                self.scattering_medium.scatterer.loss_function.addComponent(
                        Lorentz(i['params']['position'], i['params']['width'], 
                        i['params']['intensity']))
            elif i['type'] == 'VacuumExcitation':
                self.scattering_medium.scatterer.loss_function.addComponent(
                        VacuumExcitation(
                        i['params']['edge'], i['params']['fermi_width'], 
                        i['params']['intensity'], i['params']['exponent']))
        self.scattering_medium.scatterer.loss_function.step = self.step
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()
        
    def changeLossFunction(self, comp_nr, params):
        for i in params:
            setattr(self.scattering_medium.scatterer.loss_function.components[comp_nr], i, params[i])
        self.scattering_medium.scatterer.loss_function.reBuild()
        
    def dictFromScatterer(self, scatterer):
        d = {'cross_section':scatterer.cross_section ,
             'gas_diameter':scatterer.gas_diameter,
             'angle_factor':scatterer.angle_factor,
             'loss_function':[(lambda x: {'id':x, 
                                          'type':scatterer.loss_function.components[x].__class__.__name__, 
                                          'params':scatterer.loss_function.components[x].__dict__})(i) 
                                            for i in range(len(scatterer.loss_function.components))]}
        return d
    
    def updateScatterersDict(self):
        self.scatterers[self.scattering_medium.scatterer.label] = self.dictFromScatterer(self.scattering_medium.scatterer)

    def saveScatterers(self):
        self.updateScatterersDict()
        file = self.controller.scatt_datapath+"\\scatterers.json"
        with open(file, 'w') as json_file:
            json.dump(self.scatterers, json_file, indent=4)
            
    def readScatterers(self):
        file = self.controller.scatt_datapath+"\\scatterers.json"
        with open(file) as json_file: 
            scatterers = json.load(json_file)
        self.scatterers = scatterers
        self.controller.scatterer_choices = [i for i in self.scatterers]
        
    def addComponent(self, comp_kind):
        if comp_kind == 'Gauss':
            self.scattering_medium.scatterer.loss_function.addComponent(Gauss(1,1,0))
        elif comp_kind == 'Lorentz':
            self.scattering_medium.scatterer.loss_function.addComponent(Lorentz(1,1,0))
        elif comp_kind == 'VacuumExcitation': 
            self.scattering_medium.scatterer.loss_function.addComponent(VacuumExcitation(0,1,0,0.1))
        
    def addPeak(self, spec_idx, peak_kind):
        if peak_kind == 'Gauss':
            self.loaded_spectra[spec_idx].addComponent(Gauss(1,1,0))
        elif peak_kind == 'Lorentz':
            self.loaded_spectra[spec_idx].addComponent(Lorentz(1,1,0))
            
    def updatePeak(self, new_values):
        spec_idx = int(new_values['spec_idx'])
        peak_idx = int(new_values['peak_idx'])
        self.loaded_spectra[spec_idx].components[peak_idx].position = new_values['position']
        self.loaded_spectra[spec_idx].components[peak_idx].width = new_values['width']
        self.loaded_spectra[spec_idx].components[peak_idx].intensity = new_values['intensity']
        self.loaded_spectra[spec_idx].reBuild()

'''       
    def unScatterSpectrum(self):
            n = self.scattering_medium.n_iter
            a = self.scattered_spectrum.lineshape # this is the initial input spectrum
            p_inelast= self.scattering_medium.collis_prob * self.scattering_medium.scatterer.cross_section # this gives the total probability of inelastic scattering per unit distance
            self.scattering_medium.scatterer.loss_function.buildLine()
            self.scattering_medium.scatterer.loss_function.normalize()
            b = p_inelast* self.scattering_medium.scatterer.loss_function.lineshape # This rescales the loss function by the total probability per elastic collision
            self.intermediate_spectra = [a]
            for i in range(n):
                c = deconvolve(a,np.flip(b)) 
                l = len(a)
                #c = c[-l:] # this trims the unneeded data in the convolved spectrum
                a = (1+p) * a # this rescales the non-scattered portion of the spectrum
                a = np.subtract(a,c) # this sums the non-scattered and scattered spectra. it will be the input spectrum for next iteration
                self.intermediate_spectra += [a]
            self.simulated_spectrum.lineshape = a
            self.simulated_spectrum.x = self.unscattered_spectrum.x 
            self.simulated_spectrum.kind = 'Simulated'
'''