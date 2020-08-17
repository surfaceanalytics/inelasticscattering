# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json
import xlsxwriter
from model.algorithms.algorithm1 import Algorithm1
from model.algorithms.algorithm4 import Algorithm4
from model.algorithms.algorithm5 import Algorithm5
from model.algorithms.algorithm6 import Algorithm6
import datetime
from model.base_model import (Spectrum, Gauss, Lorentz, VacuumExcitation, 
                        MeasuredSpectrum, ScatteringMedium, Calculation, 
                        Tougaard, Voigt, SyntheticSpectrum)

from converters.data_converter import DataConverter
from converters.vamas import Vamas

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
        
        self.selector_table_params = {'table_name':'selector', 
                            'height':14,
                            'on_double_click':'visibility',
                            'selectmode':'extended',
                            'colour_keys':False,
                            'columns':[{'name':'Nr.', 
                                        'width':30},
                                       {'name':'Type',
                                        'width':75},
                                       {'name':'Name',
                                        'width':150}
                                       ]
                            }
        self.selector_fig_params = {'xlabel':'Energy [eV]', 
                          'ylabel': 'Intensity [cts./sec.]',
                          'axis_label_fontsize':8,
                          'title':'', 
                          'size':(4,3)}
                
        
        
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

    def loadFile(self, filename):
        """
        This function instantiates a Converter object to parse a file.
        The parsed data is stored inside the Converter's data attribute.

        Parameters
        ----------
        filename : STRING
            The name of the file to be parsed

        Returns
        -------
        INT
            The number of items in the converter's data list (i.e. the number
            of spectra in the file.

        """
        self.converter = DataConverter()
        self.converter.load(filename)
        if filename.rsplit('.')[-1] == 'vms':
            for data in self.converter.data:
                if data['settings']['y_units'] == 'counts':
                    y = np.array(data['data']['y0'])
                    dwell = data['settings']['dwell_time']
                    scans = data['scans']
                    y = y / dwell / scans
                    data['data']['y0'] = list(y)
                    data['settings']['y_units'] = 'counts_per_second'         
        return len(self.converter.data)
    
    def returnSelectorParams(self):
        return self.selector_table_params, self.selector_fig_params
    
    def returnFileContents(self):
        data = self.converter.data
        contents = [[idx, d['spectrum_type'], d['group_name']] 
                    for idx, d in enumerate(data)]
        return contents
    
    def loadSpectra(self, selection):
        """
        This function loads the selected spectra from the file into the 
        loaded_spectra attribute

        Parameters
        ----------
        selection : LIST of integers
            The indices of the Converter object's data attribute (a list), 
            representing the items in the data attribute that should be 
            stored in the Models loaded_spectra attribute.

        Returns
        -------
        None.

        """
        for idx in selection:
            x = self.converter.data[idx]['data']['x']
            y = self.converter.data[idx]['data']['y0']
            self.loaded_spectra += [MeasuredSpectrum(x,y)]
        self.converter = None
        self.start = self.loaded_spectra[-1].start    # The step width must be defined by the measured spectrum 
        self.stop = self.loaded_spectra[-1].stop       # All synthetic spectra need to have their step widths redefined
        self.step = self.loaded_spectra[-1].step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        
    def loadSpectrum(self, filename):
        selection = [self.loadFile(filename) - 1]
        self.loadSpectra(selection)
        

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
            self.simulation = Algorithm5(self.unscattered_spectrum, 
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
        
        self._onlyOneSimulated()
            
    def unScatterSpectrum(self):
        algorithm_id = 2
        params = self._getAlgorithmParams(algorithm_id)
        self._prepSpectra()
        self.simulation = Algorithm6(self.scattered_spectrum, 
                                         self.scattering_medium, params)
        simulated = self.simulation.run()
        
        self.simulated_spectrum.lineshape = simulated
        self.simulated_spectrum.x = self.scattered_spectrum.x 
        self.simulated_spectrum.kind = 'Simulated'
        self.intermediate_spectra = self.simulation.I
        
        self._onlyOneSimulated()
            
    def _onlyOneSimulated(self):
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
                {'label': 'P [mbar]', 'value':'', 'variable':'pressure', 
                 'tip':'''P represents the pressure in mbar, of the scattering
medium'''},
                {'label': 'D [mm]', 'value':'','variable':'distance',
                 'tip':'''D represents the distance from the sample to the
spectrometer nozzle, the is, the distance 
electrons travel through the scattering medium.'''},
                {'label':'Inelastic X-sect', 'value':'', 
                 'variable':'inelastic_xsect',
                 'tip':'''Inelastic X-sect is the inelastic scattering cross-
section of the scattering medium, in units of nm^2'''},
                {'label': 'f(Norm.)', 'value':'', 'variable':'norm_factor',
                 'tip': '''f(Norm.) is a normalization factor. It is used to 
account for loss of electron signal, due to 
energy loss processes not unaccounted for in 
the loss function, such as core-level excitations.'''
                 }],
            1:[
                {'label':'Inelastic Prob.', 'value':'', 
                 'variable':'inelastic_prob',
                 'tip':'''Inelastic Prob. is the probability of an inelastic
scattering event for one iteration.'''},
                {'label': 'f(Norm.)', 'value':'', 'variable':'norm_factor',
                 'tip':'''f(Norm.) is a normalization factor. It is used to 
account for loss of electron signal, due to energy loss processes not 
unaccounted for in the loss function, such as core-level excitations.'''},
                {'label':'Nr. Iter.', 'value':'','variable':'n_iter',
                 'tip':'''Nr. Iter. represents the number of iterations to 
calculate the convolution. It can be though of as a 
small delta-distance an electron travels through in 
the scattering medium.'''}]
            }

    def changeAlgorithm(self,new_id):
        self.algorithm_id = int(new_id)
        return self._returnInputFields()

    def _returnInputFields(self):
        """
        Inputs is the list of dictionaries that is exchanged between the
        model and the controller to hold the algorithm input parameters.
        Inputs_dict stores a dictionary of dictionaries, where one finds the
        dictionary of input parameters needed for each algorithm.
        The loop iterates through the list of dicts of self.inputs, then sets
        the value of the item in self.inputs to the corresponding value from 
        the object attribute, stored in Model. self.inputs will eventually be 
        sent to the controller, to be sent to the view.        

        Returns
        -------
        inputs : LIST of DICT's'
            Keys: 
                'label' : The label to be displayed in the view.
                'value' : The value to be set in the input field of the view.
                'variable' : The name of the variable used in the model.
                'tip' : The text to be displayed in the tooltip.
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
        elif algorithm_id == 2:
            params = {'n_events':self.calculation.n_events, 
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
        
    def _dictFromScatterer(self, scatterer):
        """
        This creates a dictionary of the scatterer's attributes,
        to eventually be exported as JSON.

        Parameters
        ----------
        scatterer : STRING
            The name label of the scatterer.

        Returns
        -------
        d : DICTIONARY
            Keys: 
                'inelastic_xsect' : FLOAT
                'norm_factor' : FLOAT
                'inelastic_prob' : FLOAT
                'loss_function' : DICT
                        Keys: 
                            'id' : INT
                            'type' : STRING
                            'params' : DICT
                                    keys depend on component type.

        """
        d = {'inelastic_xsect':scatterer.inelastic_xsect ,
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
        self.scatterers[label] = self._dictFromScatterer(scatterer)
        
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
        """
        Takes user selected file and load the scatters contained within that
        file. 

        Parameters
        ----------
        file : STRING
            A filename of a file in JSON format.

        Returns
        -------
        Call to the updateScattererChoices function
        which returns a list of scatterer names..

        """
        with open(file) as json_file:
            scatterers = json.load(json_file)
        self.scatterers = scatterers
        return self.updateScattererChoices() # this will return a list of scatterer choices

    def addLossComponent(self, comp_kind):
        """
        This method is for adding a new synthetic component to the loss 
        function.

        Parameters
        ----------
        comp_kind : STRING
            The name of the component class to be added as a synthetic
            component.

        Returns
        -------
        None.

        """
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
            
    def getSpecBuilderParams(self, spec_idx):
        """
        Provides all the parameters needed to build the SpecBuilder pop-up.

        Parameters
        ----------
        spec_idx : INT
            The index of the chosen spectrum

        Returns
        -------
        return_params : DICT
            A dictionary containing keys 'spec_idx', 'columns', 'entry_fields'
        """
        return_params = {'spec_idx':spec_idx,
                 'columns': self.spec_builder_columns,
                 'entry_fields': self.default_spectrum_fields}
        return return_params
    
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
        """
        A function called when a scatterers file is opened.

        Returns
        -------
        choices : LIST
            A list of strings, containing the names of the scatterers 
            available.
        """
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
        if file.rsplit('.',1)[-1] == 'vms':
            file = file.rsplit('.',1)[0]
        now = datetime.datetime.now()
        date = now.strftime("%Y-%m-%d %H:%M:%S")
        spectra = []
        for d in self.loaded_spectra:
            data = {'date': date,
                    'spectrum_type':d.kind,
                    'group_name':d.__class__.__name__,
                    'scans':1,
                    'settings':{'analysis_method':'XPS',
                                'dwell_time': 1,
                                'excitation_energy':1486.61,
                                'pass_energy':0,
                                'scan_mode':'FAT',
                                'source_label':'Al',
                                'workfunction':0,
                                'x_units':'kinetic energy',
                                'y_units':'counts',
                                },
                    'data':{'x':np.around(d.x,3),
                            'y0':np.around(d.lineshape,6)}
                    
                    }
            spectra += [data]
        self.converter = DataConverter()
        self.converter.data = spectra
        self.converter.write(file, 'Vamas')
        self.converter = None
          
                
#%%
if __name__ == "__main__":
    model = Model()
    filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Data Science\xps_data_conversion_tools\EX337 - test.vms'
    d = model.openSpectra(filepath)
    


