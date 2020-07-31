# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json
import xlsxwriter
from algorithm1 import Algorithm1
from algorithm4 import Algorithm4
import vamas as vamas
import datetime
from base_model import (Spectrum, Gauss, Lorentz, VacuumExcitation, 
                        MeasuredSpectrum, ScatteringMedium, Calculation, 
                        Tougaard, Voigt, SyntheticSpectrum)

import time
              
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
        self.component_spectra = []
        self.scatterers = {}
        self.loss_component_kinds = ['Gauss', 'Lorentz', 'VacuumExcitation',
                                     'Tougaard']
        self.peak_kinds = ['Gauss', 'Lorentz', 'Voigt']
        self.acceptance_angle = 10 # This is used in the AngularSpreadCalc class 
        
        self.calculation = Calculation()
        self.calculation.n_iter = int(100)
        
        ''' The next two attributes are for the synthetic spectrum builder
        '''
        self.default_spectrum_fields = [{'name':'position',
                                         'label':'Position',
                                         'value':0.0},
                                        {'name':'width',
                                         'label':'Width',
                                         'value':0.0},
                                        {'name':'Intensity',
                                         'label':'Intensity',
                                         'value':0.0}]
        
        self.spec_builder_columns = [{'name':'Nr.','width':30},
                                     {'name':'Type', 'width':50},
                                     {'name':'Position','width':60},
                                     {'name':'Width', 'width':60},
                                     {'name':'Intensity', 'width':60}]
       
        self.spec_builder_column_attributes = ['position', 'width', 
                                               'intensity']
        
        self.algorithm_option = ''
        self.algorithm_id = 0 # this is the default algorithm to use
        self.algorithmInputs()
        
        ''' The variable called 'variable_mapping' is used to store a 
        collection of objects and their attribute names for the parameters that 
        are exchanged with the controller for the algorithm parameter inputs'''
        self.variable_mapping = {
            'inelastic_xsect': self.scattering_medium.scatterer,
            'distance': self.scattering_medium,
            'pressure': self.scattering_medium,
            'norm_factor': self.scattering_medium.scatterer,
            'n_iter': self.calculation,
            'inelastic_prob': self.scattering_medium.scatterer,
            }

    def loadSpectrum(self, filename):
        self.loaded_spectra += [MeasuredSpectrum(filename)]
        self.start = self.loaded_spectra[-1].start    # The step width must be defined by the measured spectrum 
        self.stop = self.loaded_spectra[-1].stop      # All synthetic pectra need to have their step widths redefined
        self.step = self.loaded_spectra[-1].step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step) # Overwrite old output spectrum with new settings

    def scatterSpectrum(self):
        """ This function runs one of the calculations. The choice of calculations
        is determined from the algorithm_id. The variable 'params' stores
        all the parameters needed for the given algorithm. The algorithm is run
        and the results are stored in a simulated_spectrum object, then added
        to the list of loaded spectra.
        """
        algorithm_id = self.algorithm_id
        params = self._getAlgorithmParams(algorithm_id)
        self._prepSpectra()
        if algorithm_id == 0:
            self.simulation = Algorithm4(self.unscattered_spectrum, 
                                         self.scattering_medium, params)
            simulated = self.simulation.run()
        elif algorithm_id == 1:
            self.simulation = Algorithm1(self.unscattered_spectrum,
                                         self.scattering_medium,params)
            simulated = self.simulation.run()
                
        self.simulated_spectrum.lineshape = simulated
        self.simulated_spectrum.x = self.unscattered_spectrum.x 
        self.simulated_spectrum.kind = 'Simulated'
        self.intermediate_spectra = self.simulation.I
        
        ''' This condition ensures that there is only one simulated spectrum.'''
        if not ('Simulated' in [i.kind for i in self.loaded_spectra]):
            self.loaded_spectra += [self.simulated_spectrum]
        else:
            idx = [i.kind for i in self.loaded_spectra].index('Simulated')
            self.loaded_spectra[idx] = self.simulated_spectrum
        
    def _prepSpectra(self):
        self.scattering_medium.scatterer.step = self.unscattered_spectrum.step
        self.scattering_medium.scatterer.loss_function.reBuild()
    
    def algorithmInputs(self):
        """inputs_dict holds all the inputs needed for each algorithm, as well
        as their values, the labels to be used in the View, and the variable
        name used in the objects in Model"""
        self.inputs_dict = {
            0:[
                {'label': 'P [mbar]', 'value':'', 'variable':'pressure'},
                {'label': 'D [mm]', 'value':'','variable':'distance'},
                {'label':'Inelastic X-sect', 'value':'', 'variable':'inelastic_xsect'},
                {'label': 'f(Norm.)', 'value':'', 'variable':'norm_factor'}],
            1:[
                {'label':'Inelastic Prob.', 'value':'', 'variable':'inelastic_prob'},
                {'label': 'f(Norm.)', 'value':'', 'variable':'norm_factor'},
                {'label':'Nr. Iter.', 'value':'','variable':'n_iter'}]
            }

    def changeAlgorithm(self,new_id):
        self.algorithm_id = int(new_id)
        return self._returnInputFields()

    def _returnInputFields(self):
        """ Inputs is the list of dictionaries that is exchanged between the
        model and the controller to hold the algorithm input parameters.
        Inputs_dict stores a dictionary of dictionaries, where one finds the
        dictionary of input parameters needed for each algorithm.
        The loop iterates through the list of dicts of self.inputs, then sets
        the value of the item in self.inputs to the corresponding value from 
        the object attribute, stored in Model. self.inputs will eventually be 
        sent to the controller, to be sent to the view.
        """
        inputs = self.inputs_dict[self.algorithm_id] # this is a list of dicts
        for i in inputs:
            var = i['variable']
            obj = self.variable_mapping[var]
            i['value'] = getattr(obj,var)
        return inputs

    def updateAlgorithmParams(self, params):
        """ params comes from controller. It is a list of dictionaries
        Each dictionary containing names, and values algorithm parameters the 
        user has inputted. 
        This iterates over the items in the list, gets the name of the variable.
        The variable name maps 1:1 to the attribute name of the model object.
        The variable_mapping contains a dictionary with the variable names and
        object references for the items that need to be updated."""
        for i in params:
            var = i['variable']
            obj = self.variable_mapping[var]
            new_val = i['value']
            setattr(obj,var,new_val)

    def _getAlgorithmParams(self, algorithm_id):
        """This function prepares all of the inputs for the scattering algorithm
        calculations. Each algorithm uses a slightly different set of parameters.
        The dictionary of parameters is then returned, and eventually passed to
        the respective Algorithm class.
        """
        if algorithm_id == 0:
            params = {'n_events':self.calculation.n_events, 
                      'option': self.algorithm_option}
        elif algorithm_id == 1:
            params = {'n':self.calculation.n_iter,
                      'option': self.algorithm_option}
        return params
          
    def setCurrentScatterer(self, label):
        """ This function sets the current scatterer to whatever the variable
        'label' tells it to do. Then it looks up the parameters from the
        corresponding scatterer in the scatterers dictionary, and sets the 
        attributes of the loaded scatterer to the selected scatterer's values.
        """
        
        ''' This updates the algorithm parameters input frame.
        '''
        self.scattering_medium.scatterer.label = label
        for k, v in self.scatterers[label].items():
            if k in list(self.variable_mapping.keys()):
                var = k
                obj = self.variable_mapping[k]
                new_val = v
                setattr(obj,var,new_val)
         
        self._buildLossFromJSON(label)
        self.scattering_medium.scatterer.loss_function.step = self.step
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()
        
        return self._returnInputFields()
    
    def _buildLossFromJSON(self, label):
        """ This function builds the spectrum of the loss function from the
        components in a scatterrer loaded from JSON
        """     
        self.scattering_medium.scatterer.loss_function.components = []
        loss_fn = self.scattering_medium.scatterer.loss_function
        
        '''This is a list of dictionaries'''
        components = self.scatterers[label]['loss_function'] 
        
        for i in components:
            if i['type'] == 'Gauss':
                loss_fn.addComponent(
                        Gauss(i['params']['position'], i['params']['width'], 
                        i['params']['intensity']), rebuild = False)
            elif i['type'] == 'Lorentz':
                loss_fn.addComponent(
                        Lorentz(i['params']['position'], i['params']['width'], 
                        i['params']['intensity']), rebuild = False)
            elif i['type'] == 'VacuumExcitation':
                loss_fn.addComponent(
                        VacuumExcitation(
                        i['params']['edge'], i['params']['fermi_width'], 
                        i['params']['intensity'], i['params']['exponent']), 
                        rebuild = False)
            elif i['type'] == 'Tougaard':
                loss_fn.addComponent(
                        Tougaard(
                        i['params']['B'], i['params']['C'], 
                        i['params']['D'], i['params']['Eg']), rebuild = False)
        
    def modifyLossLineshape(self, comp_nr, params):
        """ This is run any time the user selects edits a component of the 
        loss function."""
        loss_function = self.scattering_medium.scatterer.loss_function
        loss_function.editComponent(comp_nr, params)
        
    def dictFromScatterer(self, scatterer):
        """ This creates a dictionary, to eventually be exported as JSON"""
        d = {'inelastic_xsect':scatterer.inelastic_xsect ,
             'elastic_xsect':scatterer.elastic_xsect,
             'norm_factor':scatterer.norm_factor,
             'inelastic_prob':scatterer.inelastic_prob,
             'loss_function':[(lambda x: {'id':x, 
                                          'type':scatterer.loss_function.components[x].__class__.__name__, 
                                          'params':scatterer.loss_function.components[x].__dict__})(i) 
                                            for i in range(len(scatterer.loss_function.components))]}
        return d
    
    def updateScatterersDict(self):
        """This updates the dictionary, to be exported as JSON, with any 
        changes made by the user"""
        label = self.scattering_medium.scatterer.label
        scatterer = self.scattering_medium.scatterer
        self.scatterers[label] = self.dictFromScatterer(scatterer)
        
    def newScatterer(self, name):
        self.scatterers[name] = {'inelastic_xsect':0.01,
                       'norm_factor':1,
                       'loss_function':[]}
        return self.updateScattererChoices() # this returns a list of scatterer names

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

    def addLossComponent(self, comp_kind):
        """This method is for adding a new synthetic component to the loss 
        function"""
        loss_function = self.scattering_medium.scatterer.loss_function
        if comp_kind == 'Gauss':
            loss_function.addComponent(Gauss(1,1,0))
        elif comp_kind == 'Lorentz':
            loss_function.addComponent(Lorentz(1,1,0))
        elif comp_kind == 'VacuumExcitation': 
            loss_function.addComponent(VacuumExcitation(0,1,0,0.1))
        elif comp_kind == 'Tougaard':
            loss_function.addComponent(Tougaard(200,100,500,0))
            
    def _returnComponentValues(self, spec_idx, comp_idx):
        """ This returns a list of dictionaries containing all the user-
        editable fields for a component.
        """
        comp = self.loaded_spectra[spec_idx].components[comp_idx]
        values = []
        for key in comp.__dict__.keys():
            values += [{'name':key,
                        'label':key.capitalize(),
                        'value':comp.__dict__[key]}]
        return values
            
    def addSynthSpec(self, start, stop, step):
        """ This function adds a synthetic spectrum to the list of loaded 
        spectra and returns a dictionary.
        """
        self.loaded_spectra += [SyntheticSpectrum(start, stop, step)]
        return_params = {'spec_idx':len(self.loaded_spectra)-1,
                         'columns': self.spec_builder_columns,
                         'entry_fields': self.default_spectrum_fields}
        return return_params
        
    def returnSpecComps(self, spec_idx):
        """ This function provides the component parameters for a selected 
        synthetic spectrum. It is called when a synthetic spectrum is edited.
        """
        entry_fields = self._returnComponentValues(spec_idx,0)
        return_params = {'spec_idx':spec_idx,
                         'columns': self.spec_builder_columns,
                         'entry_fields': entry_fields}
        return return_params
        
        
    def addSpecComponent(self, spec_idx, peak_kind):
        """This method adds a new synthetic component to the spectrum lineshape
        of the currently selected synthertic spectrum
        """
        spectrum = self.loaded_spectra[spec_idx]
        if peak_kind == 'Gauss':
            spectrum.addComponent(Gauss(1.0,1.0,0.0))
        elif peak_kind == 'Lorentz':
            spectrum.addComponent(Lorentz(1.0,1.0,0.0))
        elif peak_kind == 'Voigt':
            spectrum.addComponent(Voigt(1.0,1.0,0.0,0.5))
        
    def specBuilderTableData(self, spec_idx):
        table_data = {}
        for idx, comp in enumerate(self.loaded_spectra[spec_idx].components):
            table_data[idx] = ([idx, comp.__class__.__name__] + 
                           [attr for attr in comp.__dict__.values()])
        return table_data
                
    def updateComponent(self, new_values):
        """ This function updates the parameters of a synthetic spectrum
        component. It requires the spectrum id (spec_idx) and component id 
        (comp_idx) to know which component to change.
        """
        spec_idx = int(new_values['spec_idx'])
        new_values.pop('spec_idx')
        comp_idx = int(new_values['comp_idx'])
        new_values.pop('comp_idx')
        params = new_values
        spectrum = self.loaded_spectra[spec_idx]
        spectrum.editComponent(comp_idx, params)
        
    def removeSpecComp(self, spec_idx, comp_idx):
        """ This function removes a component with idx = comp_idx from the 
        synthetic spectrum with idx = spec_idx.
        """
        self.loaded_spectra[int(spec_idx)].removeComponent(int(comp_idx))
        
    def setInelasticXSect(self, inelastic_xsect):
        self.scattering_medium.scatterer.inelastic_xsect = inelastic_xsect
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()
        
    def updateScattererChoices(self):
        choices = [i for i in self.scatterers]
        return choices

    def exportToExcel(self, file):
        file = file + '.xlsx'
        workbook = xlsxwriter.Workbook(file)
        worksheet = workbook.add_worksheet()
        cols_per_spec = 2
        for i, d in enumerate(self.component_spectra):
            start_col = i * cols_per_spec
            worksheet.write(0,start_col, i)
            worksheet.write(1,start_col, 'Energy [eV]')
            worksheet.write(1,start_col+1, 'Intensity')
            x = d.x
            y = d.lineshape
            for i, val in enumerate(x):
                worksheet.write(2+i,start_col, val)
            for i, val in enumerate(y):
                worksheet.write(2+i, start_col+1, val)
        workbook.close()

    def exportToVamas(self, file):
        now = datetime.datetime.now()
        file = file + '.vms'
        data = []
        for d in self.component_spectra:
            data += [{'x':d.x,
                      'y':d.lineshape,
                      'blockID':'spectrum',
                      'sampleID':d.label,
                      'month':now.month,
                      'day':now.day,
                      'year':now.year,
                      'hour':now.hour,
                      'minute':now.minute,
                      'second':now.second,
                      'noCommentLines':0,
                      'techniqe':'XPS',
                      'sourceLabel':'Al',
                      'sourceEnergy':1486.7,
                      'sourceAnalyzerAngle':56.5,
                      'analyzerMode':'FAT',
                      'passEnergy':0,
                      'workFunction':0
                      }]
        dataset = vamas.Dataset()
        dataset.data_to_blocks(data)
        dataset.writeVamas(file)
            

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

