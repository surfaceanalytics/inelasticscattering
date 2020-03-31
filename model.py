# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json
from algorithm0 import Algorithm0
from algorithm1 import Algorithm1
from algorithm2 import Algorithm2

from base_model import Spectrum, Gauss, Lorentz, VacuumExcitation, MeasuredSpectrum, ScatteringMedium, Calculation
              
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
        self.scatterers = {}
        self.loss_component_kinds = ['Gauss', 'Lorentz', 'VacuumExcitation']
        self.peak_kinds = ['Gauss', 'Lorentz']
        self.acceptance_angle = 10 # This is used in the AngularSpreadCalc class 
        
        self.calculation = Calculation()
        self.calculation.n_iter = int(100)
        
        self.algorithm_option = ''
        self.algorithm_id = 2 # this is the default algorithm to use
        self.algorithmInputs()
        
        ''' The variable called 'variable_mapping' is used to store a collection of objects and their
        attribute names for the parameters that are exchanged with the controller
        for the algorithm parameter inputs'''
        self.variable_mapping = {
            'inelastic_xsect': self.scattering_medium.scatterer,
            'elastic_xsect': self.scattering_medium.scatterer,
            'distance': self.scattering_medium,
            'pressure': self.scattering_medium,
            'inel_decay_factor': self.scattering_medium.scatterer,
            'el_decay_factor': self.scattering_medium.scatterer,
            'n_iter': self.calculation,
            'inel_angle_factor': self.scattering_medium.scatterer,
            'el_angle_factor': self.scattering_medium.scatterer,
            'inelastic_prob': self.scattering_medium.scatterer,
            'elastic_prob': self.scattering_medium.scatterer
            }
        self._populateInputs()

    def loadSpectrum(self, filename):
        self.loaded_spectra += [MeasuredSpectrum(filename)]
        self.start = self.loaded_spectra[-1].start    # The step width must be defined by the measured spectrum 
        self.stop = self.loaded_spectra[-1].stop      # All synthetic pectra need to have their step widths redefined
        self.step = self.loaded_spectra[-1].step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step) # Overwrite old output spectrum with new settings

    def scatterSpectrum(self):
        ''' This function runs one of the calculations. The choice of calculations
        is determined from the self.algorithm_id. The variable 'params' stores
        all the parameters needed for the given algorithm. The algorithm is run
        and the results are stored in a simulated_spectrum object, then added
        to the list of loaded spectra.
        '''
        algorithm_id = self.algorithm_id
        self._prepSpectra()
        params = self._getAlgorithmParams(algorithm_id)
        if algorithm_id == 0:
            self.scattering_medium.calcDensity()
            simulation = Algorithm0(params)
            simulated = simulation.run()
        elif algorithm_id == 1:
            simulation = Algorithm1(params)
            simulated = simulation.run()
        elif algorithm_id == 2:
            self.scattering_medium.calcDensity()
            simulation = Algorithm2(params)
            inel, el, non, simulated = simulation.run()
        
        self.simulated_spectrum.lineshape = simulated
        self.simulated_spectrum.x = self.unscattered_spectrum.x 
        self.simulated_spectrum.kind = 'Simulated'
        self.intermediate_spectra = simulation.I
        
        # this ensures that there is only one simulated spectrum
        if not ('Simulated' in [i.kind for i in self.loaded_spectra]):
            self.loaded_spectra += [self.simulated_spectrum]
        else:
            idx = [i.kind for i in self.loaded_spectra].index('Simulated')
            self.loaded_spectra[idx] = self.simulated_spectrum
        
    def _prepSpectra(self):
        self.scattering_medium.scatterer.step = self.unscattered_spectrum.step
        self.scattering_medium.scatterer.loss_function.reBuild()
    
    def algorithmInputs(self):
        '''inputs_dict holds all the inputs needed for each algorithm, as well
        as their values, the labels to be used in the View, and the variable
        name used in the objects in Model'''
        self.inputs_dict = {
            0:[
                {'name': 'P [mbar]', 'value':'', 'variable':'pressure'},
                {'name': 'D [mm]', 'value':'','variable':'distance'},
                {'name':'Inelastic X-sect', 'value':'', 'variable':'inelastic_xsect'},
                {'name': 'f(Decay)', 'value':'', 'variable':'inel_decay_factor'},
                {'name': 'Elastic X-sect', 'value':'', 'variable':'elastic_xsect'},
                {'name':'f(Decay)', 'value':'','variable':'el_decay_factor'}],
            1:[
                {'name':'Inelastic Prob.', 'value':'', 'variable':'inelastic_prob'},
                {'name': 'f(Decay)', 'value':'', 'variable':'inel_decay_factor'},
                {'name': 'Elastic Prob', 'value':'', 'variable':'elastic_prob'},
                {'name':'f(Decay)', 'value':'','variable':'el_decay_factor'},
                {'name':'Nr. Iter.', 'value':'','variable':'n_iter'}],
            2:[
                {'name': 'P [mbar]', 'value':'', 'variable':'pressure'},
                {'name': 'D [mm]', 'value':'','variable':'distance'},
                {'name':'Inelastic X-sect', 'value':'', 'variable':'inelastic_xsect'},
                {'name': 'f(Angle)', 'value':'', 'variable':'inel_angle_factor'},
                {'name': 'Elastic X-sect', 'value':'', 'variable':'elastic_xsect'},
                {'name':'f(Angle)', 'value':'','variable':'el_angle_factor'}]
            }
    def changeAlgorithm(self,new_id):
        self.algorithm_id = int(new_id)
        self._populateInputs()

    def _populateInputs(self):
        ''' inputs is the list of dictionaries that is exchanged between the
        model and the controller to hold the algorithm input parameters.
        Inputs_dict stores a dictionary of dictionaries, where one finds the
        dictionary of input parameters needed for each algorithm.
        The loop iterates through the list of dicts of self.inputs, then sets
        the value of the item in self.inputs to the corresponding value from 
        the object attribute, stored in Model. self.inputs will eventually be 
        sent to the controller, to be sent to the view.
        '''
        self.inputs = self.inputs_dict[self.algorithm_id] # this is a list of dicts
        for i in self.inputs:
            var = i['variable']
            obj = self.variable_mapping[var]
            i['value'] = getattr(obj,var)
  
    def updateAlgorithmParams(self, params):
        ''' params comes from controller. It is a list of dictionaries
        Each dictionary containing names, and values algorithm parameters the 
        user has inputted. 
        This iterates over the items in the list, getss the name of the variable.
        The variable name maps 1:1 to the attribute name of the model object.
        The variable_mapping contains a dictionary with the variable names and
        object references for the items that need to be updated.'''
        for i in params:
            var = i['variable']
            obj = self.variable_mapping[var]
            new_val = i['value']
            setattr(obj,var,new_val)

    def _getAlgorithmParams(self, algorithm_id):
        '''This function prepares all of the inputs for the scattering algorithm
        calculations. Each algorithm uses a slightly different set of parameters.
        The dictionary of parameters is then returned, and eventually passed to
        the respective Algorithm class.
        '''
        if algorithm_id == 0:
            params = {'P':self.unscattered_spectrum.lineshape, # primary input spectrum
              'L': self.scattering_medium.scatterer.loss_function.lineshape * 
              self.scattering_medium.scatterer.loss_function.step,
              'I': np.array([np.zeros(len(self.unscattered_spectrum.lineshape))]),
              'elastic_xsect':self.scattering_medium.scatterer.elastic_xsect,
              'inelastic_xsect': self.scattering_medium.scatterer.inelastic_xsect,
              'density': self.scattering_medium.density,
              'distance': self.scattering_medium.distance,
              'inel_decay_factor': self.scattering_medium.scatterer.inel_decay_factor,
              'el_decay_factor': self.scattering_medium.scatterer.el_decay_factor,
              'acceptance_angle': self.acceptance_angle,
              'nr_iter_per_mfp' : self.calculation.nr_iter_per_mfp,
              'option': self.algorithm_option
              }
        elif algorithm_id == 1:
            params = {'n':self.calculation.n_iter,
              'P':self.unscattered_spectrum.lineshape,
              'L': self.scattering_medium.scatterer.loss_function.lineshape * 
              self.scattering_medium.scatterer.loss_function.step,
              'I': np.array([np.zeros(len(self.unscattered_spectrum.lineshape))]),
              'elastic_prob':self.scattering_medium.scatterer.elastic_prob,
              'inelastic_prob': self.scattering_medium.scatterer.inelastic_prob,
              'inel_decay_factor': self.scattering_medium.scatterer.inel_decay_factor,
              'el_decay_factor': self.scattering_medium.scatterer.el_decay_factor,
              'option': self.algorithm_option
              }
        elif algorithm_id == 2:
            params = {'n':self.calculation.n_events,
              'P':self.unscattered_spectrum.lineshape, # primary input spectrum
              'L': self.scattering_medium.scatterer.loss_function.lineshape * 
              self.scattering_medium.scatterer.loss_function.step,
              'I': np.array([np.zeros(len(self.unscattered_spectrum.lineshape))]),
              'elastic_xsect':self.scattering_medium.scatterer.elastic_xsect,
              'inelastic_xsect': self.scattering_medium.scatterer.inelastic_xsect,
              'density': self.scattering_medium.density,
              'distance': self.scattering_medium.distance,
              'inel_angle_factor': self.scattering_medium.scatterer.inel_angle_factor,
              'el_angle_factor': self.scattering_medium.scatterer.el_angle_factor,
              'acceptance_angle': self.acceptance_angle,
              'option': self.algorithm_option
              }
        return params
          
    def setCurrentScatterer(self, label):
        ''' This function sets the current scatterer to whatever the variable
        'label' tells it to do. Then it looks up the parameters from the
        corresponding scartterer in the scatterers dictionary, and sets the 
        attributes of the loaded scatterer to the selected scatterer's values.'''
        self.scattering_medium.scatterer.label = label
        for k, v in self.scatterers[label].items():
            if k in list(self.variable_mapping.keys()):
                var = k
                obj = self.variable_mapping[k]
                new_val = v
                setattr(obj,var,new_val)
        self._populateInputs()
        
        '''This part builds the spectrum of the loss function
        '''
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
             'inel_decay_factor':scatterer.inel_decay_factor,
             'el_decay_factor':scatterer.el_decay_factor,
             'inelastic_prob':scatterer.inelastic_prob,
             'elastic_prob':scatterer.elastic_prob,
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

    def loadScatterers(self, file):
        """ take user selected file and load the scatters contained within that
        file """
        with open(file) as json_file:
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
    model.scattering_medium.distance = 0.8
    #model.scattering_medium.density = 0.005
    model.scattering_medium.scatterer.inelastic_xsect = 0.01
    model.scattering_medium.scatterer.inelastic_xsect = 0.01
    model.scattering_medium.scatterer.inel_angle_factor = 30
    model.scattering_medium.scatterer.el_angle_factor = 360
    #model.n_events = 50
    model.algorithm_id = 0
    model.scatterSpectrum()
    
    I = model.intermediate_spectra
#%%
        
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
    
