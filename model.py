# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json
import re
from angular_spread_lorentzian import AngularSpreadCalc

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
        self.loss_function = LossFunction(0,200,0.1)
        self.cross_section = 0.01
        self.gas_diameter = 0.2 #In nanometers
        self.gas_cross_section = np.pi * (self.gas_diameter / 2)**2
        self.elastic_xsect = 0.01 # in units of nm^3
        self.inelastic_xsect = 0.01 # in units of nm^3
        self.inel_angle_factor = 1 
        self.el_angle_factor = 1
        
class ScatteringMedium():
    def __init__(self):
        self.scatterer = Scatterer()
        self.R = 8.314463E+25 # gas constant in units of nm^3.mbar.K^-1.mol^-1
        self.avagadro = 6.022141E+23 # Avagadro's contant
        self.T = 300 # temperature in Kelvin
        self.pressure = 1 # In mbar
        self.density = self.pressure / (self.R * self.T) * self.avagadro # molecular density in units of particles per nm^3
        self.d_through_gas = 800000 # In nanometers
        
    def calcDensity(self):
        self.density = self.pressure / (self.R * self.T) * self.avagadro # units are particles per nm^3

    def calcMFP(self):
        if (self.scatterer.elastic_xsect !=0) & (self.scatterer.inelastic_xsect !=0):
            self.mfp = 1/(self.scatterer.elastic_xsect * self.density) # the elastic mean free path in nm
            self.imfp = 1/(self.scatterer.inelastic_xsect * self.density) # the inelastic mean free path in nm
            self.d_mfp = self.d_through_gas / min(self.mfp,self.imfp) # this is the distance from sample to spectrometer in units of mfp
        elif self.scatterer.elastic_xsect == 0:
            self.imfp = 1/(self.scatterer.inelastic_xsect * self.density)    
            self.d_mfp = self.d_through_gas / abs(self.imfp) # this is the distance from sample to spectrometer in units of mfp
        elif self.scatterer.inelastic_xsect == 0:
            self.imfp = 1/(self.scatterer.elastic_xsect * self.density)    
            self.d_mfp = self.d_through_gas / abs(self.mfp)
            
        #print('mfp: '+str(self.mfp))
        #print('imfp: '+str(self.imfp))
        #print('d_mfp: '+str(self.d_mfp))
        #print('n_iter: '+str(self.n_iter))
        #print('d_iter:'+str(self.d_iter))
        #print('inelast prob: '+str(self.inelastic_prob))
        #print('elast prob: '+str(self.elastic_prob))
        
    def setPressure(self, pressure):
        self.pressure = pressure
        self.calcDensity()
        
    def setDistance(self,distance):
        self.d_through_gas = distance
              
class Model():
    def __init__(self):
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
        self.scatterers = {}
        self.loss_component_kinds = ['Gauss', 'Lorentz', 'VacuumExcitation']
        self.peak_kinds = ['Gauss', 'Lorentz']
        self.calc_variant = 0
        self.n_events = 50
        
        self.acceptance_angle = 10
        self.n_iter = int(100)
        self.nr_iter_per_mfp = 10

    def calcParams(self):
        if self.scattering_medium.d_mfp < 1:
            self.n_iter = 1
        else:
            self.n_iter = int(self.scattering_medium.d_mfp * self.nr_iter_per_mfp)# this is the number of iterations to use            
        self.d_iter = self.scattering_medium.d_through_gas / self.n_iter # this is the distance for one iteration in nm       
        self.inelastic_prob = 1 - np.exp(-1 * (self.scattering_medium.d_through_gas / self.n_iter) * self.scattering_medium.scatterer.inelastic_xsect * self.scattering_medium.density)
        self.elastic_prob = 1 - np.exp(-1 * (self.scattering_medium.d_through_gas / self.n_iter) * self.scattering_medium.scatterer.elastic_xsect * self.scattering_medium.density)

    def loadSpectrum(self, filename):
        self.loaded_spectra += [MeasuredSpectrum(filename)]
        self.start = self.loaded_spectra[-1].start    # The step width must be defined by the measured spectrum 
        self.stop = self.loaded_spectra[-1].stop      # All synthetic pectra need to have their step widths redefined
        self.step = self.loaded_spectra[-1].step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step) # Overwrite old output spectrum with new settings

    def scatterSpectrum(self, version):
        self.prepSpectra()
        if version == 0:
            self.scattering_medium.calcDensity()
            self.scattering_medium.calcMFP()
            self.calcParams()
            self.algorithm1()
        elif version == 1:
            self.algorithm1()
        elif version == 2:
            self.scattering_medium.calcDensity()
            self.algorithm2()
            
    def prepSpectra(self):
        self.scattering_medium.scatterer.step = self.unscattered_spectrum.step
        self.scattering_medium.scatterer.loss_function.reBuild()
        
    def algorithm1(self):
        n = self.n_iter
        a = self.unscattered_spectrum.lineshape # this is the initial input spectrum
        p_elast = self.elastic_prob # this is the probability of elastic scattering per iteration
        p_inelast= self.inelastic_prob # this gives the total probability of inelastic scattering per unit distance
        b = p_inelast * self.scattering_medium.scatterer.loss_function.lineshape * self.scattering_medium.scatterer.loss_function.step # This rescales the loss function by the total probability per elastic collision
        self.intermediate_spectra = [a]
        for i in range(n):
            c = np.convolve(a,np.flip(b)) # this convolves the input spectrum with the scaled loss function
            l = len(a)
            c = c[-l:]  # this trims the unneeded data in the convolved spectrum
            c = c * self.scattering_medium.scatterer.inel_angle_factor # this rescales the scattered spectrum by the angle factor
            #print(self.scattering_medium.scatterer.inel_angle_factor)
            #print('area of inelastic spectrum: ' + str(np.sum(c) * self.scattering_medium.scatterer.loss_function.step))
            a = (1 - p_elast - p_inelast + p_elast * self.scattering_medium.scatterer.el_angle_factor) * a # this rescales the non-scattered portion of the spectrum
            #print('area of unscattered spectrum: ' + str(np.sum(a)*self.unscattered_spectrum.step))
            a = np.add(c,a) # this sums the non-scattered and scattered spectra. it will be the input spectrum for next iteration
            #print('area of output spectrum: ' + str(np.sum(a)*self.unscattered_spectrum.step))
            self.intermediate_spectra += [a]
        self.simulated_spectrum.lineshape = a
        self.simulated_spectrum.x = self.unscattered_spectrum.x 
        self.simulated_spectrum.kind = 'Simulated'
        self.bulk_spectrum = np.sum(self.intermediate_spectra, axis=0)
        
    def algorithm2(self):
        self.simulated_spectrum = SimulatedSpectrum(self.start, self.stop, self.step)
        n = self.n_events # number of scattering events to simulate
        P = self.unscattered_spectrum.lineshape # primary input spectrum
        L = self.scattering_medium.scatterer.loss_function.lineshape * self.scattering_medium.scatterer.loss_function.step # This rescales the loss function by the total probability per elastic collision
        self.intermediate_spectra1 = {'unscattered':P, 'el_scattered':np.zeros(len(P)),'inel_scattered':np.array([np.zeros(len(P))])}
        I = np.array([np.zeros(len(P))]) # the matrix of inelastically scattered spectra
        for i in range(1,n):
            # this condition is so that the first iteration uses the primary spectrum as input
            # all subsequent iterations use the convolved spectrum from the previous iteration
            if i == 1:
                new = np.convolve(P,np.flip(L))
            else:
                new = np.convolve(I[-1],np.flip(L))
            l = len(P)
            new = new[-l:]  # this trims the unneeded data in the convolved spectrum
            new = np.array([new])
            I = np.concatenate((I,new),axis=0)
        
        self.simulated_spectrum.convolved = I
        
        # these are the needed parameters for the Poisson distribution
        elastic_xsect = self.scattering_medium.scatterer.elastic_xsect
        inelastic_xsect = self.scattering_medium.scatterer.inelastic_xsect
        density = self.scattering_medium.density
        distance = self.scattering_medium.d_through_gas 
           
        def Poisson(n, distance, sigma, density):
            p = (1/np.math.factorial(n)) * ((distance * density * sigma)**n) * np.exp(-1 * distance * density * sigma)
            return p

        # These are the Poisson distributions for inelastic and elastic scattering
        self.simulated_spectrum.poiss_inel = np.array([Poisson(i,distance,inelastic_xsect,density) for i in range(n)])
        self.simulated_spectrum.poiss_el = np.array([Poisson(i,distance,elastic_xsect,density) for i in range(n)])
        
        # calculate the dot product of the two Poisson distributions
        # rows represent number of times electron is inelastically scattered
        # columns represent number of time an electron was elastically scattered
        T = np.dot(np.array([self.simulated_spectrum.poiss_inel]).T, np.array([self.simulated_spectrum.poiss_el]))
        
        # calculate the scattering probabilities for elastic, inelastic and non-scattering
        self.simulated_spectrum.p_non = T[0,0] # this is the probability of not being scattered over the distance d
        self.simulated_spectrum.p_el = T[0,1:]
        self.simulated_spectrum.p_inel = np.sum(T[1:,:], axis=1)

        # calculate the scaling factors due to angular spread for inelastic and elastic scattering
        self.simulated_spectrum.inel_angle = self.calcAngleDist(width = self.scattering_medium.scatterer.inel_angle_factor)[0:-1]
        self.simulated_spectrum.el_angle = self.calcAngleDist(width = self.scattering_medium.scatterer.el_angle_factor)[0:-1]

        # an array of the angular spread factors for inelastic and elastic scattering
        A = np.dot(np.array([self.simulated_spectrum.inel_angle]).T,np.array([self.simulated_spectrum.el_angle]))
        
        # Multipliy element wise the angle factors and the Poisson factors
        F = T * A
        
        # for the case of scattering though film
        inel_factor = np.sum(F[1:,:], axis=1)
        el_factor = np.sum(F[0,1:])
        non_factor = F[0,0]
        inel = np.array(np.dot(I[1:].T,np.array([inel_factor]).T))[:,0] # scale all inelastic scattered spectra by their probabilities
        el = P * el_factor # scale the primary spectrum by the elastic probability for the elastically scattered signal
        non = P * non_factor # scale the primary spectrum by the probability for no scattering 
        simulated = np.sum([np.array(inel),np.array(el), np.array(non)], axis=0)
        self.simulated_spectrum.film['inel'] = inel
        self.simulated_spectrum.film['el'] = el
        self.simulated_spectrum.film['non'] = non
        self.simulated_spectrum.film['sum'] = simulated
        
        self.simulated_spectrum.lineshape = simulated
        self.simulated_spectrum.x = self.unscattered_spectrum.x 
        self.simulated_spectrum.kind = 'Simulated'

        # for the case of scattering though bulk
        imfp = density * inelastic_xsect
        mfp = density * elastic_xsect
        inel_factor = self.simulated_spectrum.inel_angle * imfp
        el_factor = np.sum(self.simulated_spectrum.el_angle * mfp)
        inel = np.array(np.dot(I[1:].T,np.array([inel_factor]).T[1:]))[:,0] # scale all inelastic scattered spectra by their probabilities
        el = P * el_factor # scale the primary spectrum by the elastic probability for the elastically scattered signal
        non = P * (imfp + mfp) / 2 # scale the primary spectrum by the probability for no scattering 
        simulated = np.sum([np.array(inel),np.array(el), np.array(non)], axis=0)
        self.simulated_spectrum.bulk['inel'] = inel
        self.simulated_spectrum.bulk['el'] = el
        self.simulated_spectrum.bulk['non'] = non
        self.simulated_spectrum.bulk['sum'] = simulated

        self.inel = inel
        self.el = el
        self.non = non
        self.I = I
        self.P = P
        
        self.bulk_spectrum = np.add(np.sum(I, axis=0), P)
        
    def calcAngleDist(self, width):
        angle_dist = AngularSpreadCalc(iterations=self.n_events,
                                acceptance_angle=self.acceptance_angle,
                                energy=500,
                                width=width)

        # generate lorenztian
        angle_dist.gen_lorentzian_cross_section()
        
        # run convolution
        angle_dist.run_convolution()
        # limit by acceptance angle
        angle_dist.limit_by_acceptance_angle()
        # calculate accepted area under curve
        area_sum = angle_dist.calc_area_under_curve()
        
        return area_sum
        
    def setCurrentScatterer(self, label):
        self.scattering_medium.scatterer.label = label
        self.scattering_medium.scatterer.inelastic_xsect = self.scatterers[label]['inelastic_xsect']
        self.scattering_medium.scatterer.elastic_xsect = self.scatterers[label]['elastic_xsect']
        self.scattering_medium.scatterer.inel_angle_factor = self.scatterers[label]['inel_angle_factor']
        self.scattering_medium.scatterer.el_angle_factor = self.scatterers[label]['el_angle_factor']
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
        d = {'inelastic_xsect':scatterer.inelastic_xsect ,
             'elastic_xsect':scatterer.elastic_xsect,
             'inel_angle_factor':scatterer.inel_angle_factor,
             'el_angle_factor':scatterer.el_angle_factor,
             'loss_function':[(lambda x: {'id':x, 
                                          'type':scatterer.loss_function.components[x].__class__.__name__, 
                                          'params':scatterer.loss_function.components[x].__dict__})(i) 
                                            for i in range(len(scatterer.loss_function.components))]}
        return d
    
    def updateScatterersDict(self):
        self.scatterers[self.scattering_medium.scatterer.label] = self.dictFromScatterer(self.scattering_medium.scatterer)
        
    def newScatterer(self, name):
        self.scatterers[name] = {'inelastic_xsect':0.01,'elastic_xsect':0.01,'inel_angle_factor':1,'el_angle_factor':1,'loss_function':[]}
        return self.updateScattererChoices() # this will return a list of scatterer chouces

    def saveScatterers(self, file):
        self.updateScatterersDict()
        if file[-5:] != '.json':
            filename = file + '.json'
        else:
            filename = file
        with open(filename, 'w') as json_file:
            json.dump(self.scatterers, json_file, indent=4)

    def loadScatterers(self, loss_fn_file):
        """ take user selected file and load the scatters contained within that
        file """
        with open(loss_fn_file) as json_file:
            scatterers = json.load(json_file)
        self.scatterers = scatterers
        return self.updateScattererChoices() # this will return a list of scatterer choices

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
        
    def setInelAngle(self, new_value):
        self.scattering_medium.scatterer.inel_angle_factor = new_value
        
    def setElAngle(self, new_value):
        self.scattering_medium.scatterer.el_angle_factor = new_value
        
    def setInelasticXSect(self, inelastic_xsect):
        self.scattering_medium.scatterer.inelastic_xsect = inelastic_xsect
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()
        
    def setElasticXSect(self, elastic_xsect):
        self.scattering_medium.scatterer.elastic_xsect = elastic_xsect
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()

    def updateScattererChoices(self):
        choices = [i for i in self.scatterers]
        return choices
        

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

#%%
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    model = Model()
    filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\He\ARM22\Ag3d vacuum EX340 ARM22.TXT'
    model.loadSpectrum(filename1)
    filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
    model.loadScatterers(filename2)
    
    x = model.loaded_spectra[0].x
    y = model.loaded_spectra[0].lineshape
    
    model.setCurrentScatterer('He')
    #plt.plot(x,y)
    model.unscattered_spectrum = model.loaded_spectra[0]
    model.scattering_medium.setPressure(4)
    #model.scattering_medium.d_through_gas = 50000
    #model.scattering_medium.density = 0.005
    model.scattering_medium.scatterer.inelastic_xsect = 0.01
    model.scattering_medium.scatterer.inelastic_xsect = 0.01
    model.scattering_medium.scatterer.inel_angle_factor = 30
    model.scattering_medium.scatterer.el_angle_factor = 360
    model.n_events = 50
    model.algorithm2()
        
    limit = 6
    #plt.plot(model.poiss_inel[:limit])
    #plt.plot(model.poiss_el[:limit])

    plt.plot(model.simulated_spectrum.lineshape)
    
    plt.plot(model.simulated_spectrum.film['sum'])

    
    '''plt.plot(model.simulated_spectrum.x,model.simulated_spectrum.lineshape)

    for i in range(6):
        plt.plot(model.I[i,:])
    plt.show()

    plt.plot(model.poiss_inel[1:11])
    plt.plot(model.inel_angle[:10])
    plt.plot(model.inel_factor[:10])
    plt.show()
    
    print(model.poiss_inel[1])
    print(model.inel_factor[0])'''


#%%
    
