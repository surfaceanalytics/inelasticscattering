# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:21:27 2020

@author: Mark
"""
from base_model import SyntheticSpectrum
from model import Model
from view import View, LossEditor, SpecBuilder
import numpy as np
import os
import tkinter
import tkinter.ttk as tk
from tkinter import PhotoImage, Menu
import threading
from concurrent import futures

class Controller():
    def __init__(self):
        self.root = tkinter.Tk()

        self.fig2_list = {'loss function':[]}
        
        self.scatterer_choices = []
        self.default_algorithm = 0
        
        self.thread_pool_executer = futures.ThreadPoolExecutor(max_workers=1)
        
        self.press = None
        self.moved_while_pressed = 0
        
        self.spectra_table_choices = {1:['none','Scattered','Unscattered'],
                                      2:['visible','hidden']}
        self.table1_entries = []
        self.datapath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\data'
        self.resourcepath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\resources'
        self.scatt_datapath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\data'
        
        self.model = Model()
        
        self.component_choices = self.model.loss_component_kinds
        self.peak_choices = self.model.peak_kinds
        
        self.view = View(self, self.root)
        
        self._sendVariantChoices()

        self.img_collection = []

        self.colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        self.table_names = ['spectra','loss_function', 'selector']
                                
    def mainloop(self):
        self.root.mainloop()
        
    def _sendVariantChoices(self):
        """ This send the choices of algorithm variant to the view."""
        params = list(self.model.inputs_dict.keys())
        self.view.addVariantChoices(params)
        '''set variant default'''
        self.view.variant.set(self.default_algorithm)
        self.toggleVariant(self.default_algorithm)
        
    def toggleVariant(self, new_var):
        """ This function sete the new algorithm_id in model, then gets the
        Algorithm input fields for the chosen algorithm."""
        params = self.model.changeAlgorithm(new_var)
        self.sendInputsToView(params)
        
    def sendInputsToView(self, params):
        """
        This sends the user inputs for the algorithm choice to the view.
        Parameters
        ----------
        params: LIST of DICT's'
            Keys: 
                'label' : The label to be displayed in the view.
                'value' : The value to be set in the input field of the view.
                'variable' : The name of the variable used in the model.
                'tip' : The text to be displayed in the tooltip.'''

        Returns
        -------
        None.
        """
        self.view.buildAlgorithmFrame(params)
                
    def updateValuesInModel(self):
        ''' Get new values from view'''
        params = self.view.inputs_frame.sendValues()
        ''' Send new values to model'''
        self.model.updateAlgorithmParams(params)
        
    def getInputsFromModel(self):
        """ This function gets the input fields for the currently selected 
        algorithm.
        """
        params = self.model.inputs
        self.sendInputsToView(params)
                
    def scatterSpectrum(self):
        self.updateValuesInModel()
        ''' This makes sure there it is defined which spectrum should be used as 
        the input spectrum'''
        if (('Unscattered' not in [x.kind for x in self.model.loaded_spectra]) 
            | (len(self.model.scattering_medium.scatterer.loss_function.components)==0)):
            self.view.noScatterer()
        else:
            ''' This checks if the user wants to display the 'bulk' spectrum'''
            if self.view.bulk.get() == 1:
                self.model.algorithm_option = 'bulk'
            else:
                self.model.algorithm_option = 'film'
            
            self.model.scatterSpectrum()
            self.view.tables[0].fillTable()
            self.rePlotFig(1, rescale=False)
            
    def unScatterSpectrum(self):
        self.updateValuesInModel()
        if (('Scattered' not in [x.kind for x in self.model.loaded_spectra]) 
            | (len(self.model.scattering_medium.scatterer.loss_function.components)==0)):
            self.view.noScatterer()
        else:
            ''' This checks if the user wants to display the 'bulk' spectrum'''
            if self.view.bulk.get() == 1:
                self.model.algorithm_option = 'bulk'
            else:
                self.model.algorithm_option = 'film'
            
            self.model.unScatterSpectrum()
            self.view.tables[0].fillTable()
            self.rePlotFig(1, rescale=False)
        
                
    def setScatteredSpectrum(self, spectrum):
        self.model.scattered_spectrum = spectrum
        
    def setUnscatteredSpectrum(self, spectrum):
        self.model.unscattered_spectrum = spectrum
        self.model.scattering_medium.scatterer.loss_function.step = spectrum.step
        self.model.scattering_medium.scatterer.loss_function.reBuild()
        self.rePlotFig(2)
        
    def loadFile(self, filename):
        """
        This function is clled when the user selects to open a file.
        It then calls the Model's LoadFile method to parse the file.
        If the file contains more than one spectrum, the Controller
        calles the SpectrumSelector in the View.

        Parameters
        ----------
        filename : STRING
            The filename of the file to be parsed.

        Returns
        -------
        None.

        """
        n_spectra = self.model.loadFile(filename)
        if n_spectra > 1:
            table_params, fig_params = self.model.returnSelectorParams()
            self.view.callSpecSelector(table_params, fig_params)
        else:
            self.loadSpectra([0])

    def loadSpectra(self, selection):
        """
        This function is called once it is clear which spectra the user
        wishes to load.

        Parameters
        ----------
        selection : INT
            Represents the index of the spectrum in the set of spectra.

        Returns
        -------
        None.

        """
        selection = list(selection)
        selection = [int(s) for s in selection]
        self.model.loadSpectra(selection)
        self.rePlotFig(1)
        self.view.tables[0].fillTable()
        
    def rePlotFig(self, fig_nr, **kwargs):
        """
        This tells the View to replot a figure.

        Parameters
        ----------
        fig_nr : INT
            A reference to the list of figures available.
            The 
        **kwargs : Keyword arguments
            Handled values are 'rescale' = True/False

        Returns
        -------
        None.

        """
        params = {}
        
        if 'rescale' in kwargs.keys():
            params['rescale'] = kwargs['rescale']
            
        if fig_nr == 1:
            ''' This replots figure 1 and re-sets the x and y limits'''
            ''' First get all visible spectra'''
            data = [{'x':i.x, 'y':i.lineshape, 'idx':idx} 
                    for idx,i in enumerate(self.model.loaded_spectra)
                    if i.visibility == 'visible']
            self.view.fig1.rePlotFig(data, params)
            
        elif fig_nr == 2:
            ''' This replots figure 2 and re-sets the x and y limits'''
            ''' First get all visible spectra'''
            params['normalize'] = 0
            loss_function = self.model.scattering_medium.scatterer.loss_function
            x = loss_function.x
            y = loss_function.lineshape
            data = [{'x':x, 'y':y, 'idx':0}]
            self.view.fig2.rePlotFig(data, params)
        
        elif fig_nr == 3:
            ''' This replots the figure in the SpecSelector.'''
            selection = list(kwargs['selection'])
            selection = [int(s) for s in selection]
            data = [{'x':d['data']['x'], 'y':d['data']['y0'], 'idx':idx} 
                    for idx, d in enumerate(self.model.converter.data)
                    if idx in selection]
            params['rescale'] = True
            self.view.specselector.figure.rePlotFig(data, params)
            
    def changeNormalization(self, normalized):
        if normalized == 1:
            self.view.figures[0].normalized = True
            self.rePlotFig(1, rescale = True)
        else:
            self.view.figures[0].normalized = False
            self.rePlotFig(1, rescale = True)
        
    def tableEdit(self, params):
        """
        This function is called from the table edit function.

        Parameters
        ----------
        params : DICT
            Keys: 
                'table_name': The name of the table that called this function
                'idx': The index of the row that this function was called from    
        Returns
        -------
        None.

        """
        table_name = params['table_name']
        idx = int(params['idx'])
        
        if table_name == self.table_names[0]:
            ''' This is run if the table is the spectra table.'''
            spec_type = self.getSpectrumType(idx)
            #print(spec_type)
            if spec_type == 'MeasuredSpectrum':
                pass
            elif spec_type == 'SyntheticSpectrum':
                params = self.model.getSpecBuilderParams(idx)
                
                self.view.spec_builder = SpecBuilder(self, params)
        elif table_name == self.table_names[1]:
            self._callLossEditor(idx)
            
    def tableChoice(self, params):
        #print(params)
        table_name = params['table_name']
        idx = params['idx']
        column = params['column']
        choice = params['choice']
        
        if table_name == self.table_names[0]:
            ''' This is run if the table is the spectra table.'''
            self._spectrumTableLogic(idx, column, choice)
        
    def _spectrumTableLogic(self, idx, col, choice):
        '''Only one spectrum can be of kind 'scattered' and one of kind 
        'unscattered' '''
        
        ''' There can only be one 'Scattered' spectrum and only one 
        'Unscattered' spectrum.
        First check if the choice is 'Scattered' or 'Unscattered'
        '''
        if (choice in['Scattered','Unscattered']) & (col == 'Type'):
            ''' Then check if there already exists a spectrum in the loaded
            spectra that is the same as the choice.
            '''
            if any(x.kind == choice for x in self.model.loaded_spectra):
                ''' If the choice label already exists in the loaded spectra
                then get the id of the spectrum.
                '''
                indexes = [self.model.loaded_spectra.index(x) 
                           for x in self.model.loaded_spectra 
                           if x.kind == choice]
                for i in indexes:
                    ''' Then change the spectrum's attribute to 'none'.
                    '''
                    self.model.loaded_spectra[i].kind = 'none'
            ''' Then assign the kind attribute of the spectrum with id == idx
            to the value of choice.
            '''
            self.model.loaded_spectra[idx].kind = choice
            
            ''' It is allowed to have multiple spectra of kind == 'none'.
            So here we can just assign kind attribute to 'none' for the 
            spectrum with id == idx.
            '''
        elif choice == 'none':
            self.model.loaded_spectra[idx].kind = choice
            
        ''' Then tell the Model to load the spectrum from the loaded_spectra
        list, where the kind attribute == 'Scattered' or 'Unscattered', to the
        Scattered spectrum or Unscattered spectrum, which is used as imput for
        the algorithms.
        '''
        for spectrum in self.model.loaded_spectra:
            if spectrum.kind == 'Scattered':
                self.setScatteredSpectrum(spectrum)
            elif spectrum.kind == 'Unscattered':
                self.setUnscatteredSpectrum(spectrum)
                
    def getTableData(self, table_name):
        """
        Parameters
        ----------
        table_name : string
            The name of the table.

        Returns
        -------
        A list of dictionaries.

        """
        data = []
        if table_name == self.table_names[0]:
            for idx, row in enumerate(self.model.loaded_spectra):
                data += [[idx, row.kind, row.visibility, 'edit']]
        elif table_name == self.table_names[1]:
            components = self.model.scattering_medium.scatterer.loss_function.components
            for i,j in enumerate(components):
                data += [[i, type(j).__name__]]
        elif table_name == self.table_names[2]:
            data = self.model.returnFileContents()
        return data

    def getSpectrumType(self, idx):
        """
        Gets the spectrum type for the spectrum with id == idx.

        Parameters
        ----------
        idx : INT
            The index os the spectrum in the loaded_spectra list.

        Returns
        -------
        TYPE
            The name of the spectrum type.
        """
        return self.model.loaded_spectra[idx].__class__.__name__
        
    def setCurrentScatterer(self, label):
        """This is called from the View when the user selects a different 
        scatterer. The Controller retrieves the scatterer's parameters
        from the Model and sends them to the View.
        """
        params = self.model.setCurrentScatterer(label)
        self.sendInputsToView(params)
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.rePlotFig(2)
        self.view.tables[1].fillTable()
    
    def _returnCompType(self, comp):
        return str(comp.__class__).split(".")[1].strip("'>")
    
    def _callLossEditor(self, idx):
        """ The loss editor is a pop-up window with entry fields for the 
        parameters of the loss function component"""
        comp = self.model.scattering_medium.scatterer.loss_function.components[idx]
        comp_type = self._returnCompType(comp)
        params = comp.__dict__
        self.spec_builder = LossEditor(self, params, idx, comp_type)
    
    def modifyLossLineshape(self):
        """
        This function is called from the View any time the user
        modifies a component in the loss function.
        It first gets the idx and parameters of the selected component.
        Then it passes the idx and parameters tot he Model and tells the
        Model to reconstruct the loss function.
        Then it tells the View to replot the Loss Function figure.

        Returns
        -------
        None.

        """
        comp_nr = self.spec_builder.comp_nr
        params = self.spec_builder.params
        loss_function = self.model.scattering_medium.scatterer.loss_function
        self.model.modifyLossLineshape(comp_nr, params)
        self.fig2_list['loss function']=[]
        self.fig2_list['loss function']=[loss_function.x, 
                                         loss_function.lineshape]
        self.rePlotFig(2, rescale = False)
        
    def changeVisibility(self, params):
        """
        This toggles the visibility of spectra in a figure by double 
        clicking the left mouse button.
        Parameters
        ----------
        params : DICTIONARY
            The keys are:
                'table_name', indicating which table made the request,
                'idx', indicating the index of the table row,
                'visibility' indicating the new viibility state.
            
        Returns
        -------
        None.

        """
        #print(params)
        idx = params['idx']
        table_name = params['table_name']
        visibility = params['visibility']

        if table_name == self.table_names[0]:
            self.model.loaded_spectra[idx].visibility = visibility
            ''' Argument passed here is the figure number.'''
            if visibility == 'visible':
                self.rePlotFig(1, rescale=True)
            elif visibility == 'hidden':
                self.rePlotFig(1, rescale=False)
            ''' Argument passed here is the table number.'''
        else:
            pass
        
    def updateScatterersDict(self):
        self.model.updateScatterersDict()

    def loadScatterers(self, file):
        """ Pass file to model to load the scatterers within, update view with
        scatter_choices"""
        # Load scatterers from file
        scatterer_choices = self.model.loadScatterers(file)
        # Update choice in menu
        self.view.updateScattererChoices(scatterer_choices)
        self.view.selected_scatterer.set(scatterer_choices[0])
        self.setCurrentScatterer(scatterer_choices[0])

    def saveScatterers(self, file):
        self.model.saveScatterers(file)        
        
    def removeRow(self, params):
        """
        This function is called from a Table class every time a row is deleted.

        Parameters
        ----------
        params : DICT
            Keys:
                'table_name' : The name of the table.
                'idx': The index of the row.

        Returns
        -------
        None.

        """
        table_name = params['table_name']
        idx = params['idx']
        if table_name == self.table_names[0]:
            self._removeSpectrum(idx)
        elif table_name == self.table_names[1]:
            self._removeLossComponent(idx)
        
    def _removeSpectrum(self, idx):
        """ This removes a spectrum from the loaded spectra. It is called
        from the View when the user preses <delete> on a row in the spectra
        table. The Controller tells the Model to remove a spectrum from the 
        list of loaded spectra. Then the Controller tells the view to refresh
        the spectra figure and the spectra table.
        
        Parameters
        ----------
        idx : INT
                The index of the spectrum to be removed.
        Returns
        -------
        None
        """
        del self.model.loaded_spectra[idx]
        self.rePlotFig(1, rescale=False)
        
    def _removeLossComponent(self, comp_idx):
        """ This removes a component from the loss function. It is called when
        the user presses <delete> on a selected row in the loss function table.
        The Controller tells the Model to remove one component from the 
        currently loaded loss function. Then the Controller tells the View to 
        replot the loss function and refresh the loss function table.
        """
        self.model.scattering_medium.scatterer.loss_function.removeComponent(comp_idx)
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.rePlotFig(2, rescale=False)
        
    def addComponent(self, comp_kind):
        """This adds a new component to the loss function. It is called when
        the user presses 'add component' button on the loss functions
        panel, and selects a type of loss function component.
        The Controller tells the Model to add one component of a chosen type
        to the currently loaded loss function. Then the Controller tells the 
        View to replot the loss function and refresh the loss function table.
        """
        self.model.addLossComponent(comp_kind)
        scatterer = self.model.scattering_medium.scatterer
        self.fig2_list['loss function'] = [scatterer.loss_function.x,
                      scatterer.loss_function.lineshape]
        comp = scatterer.loss_function.components[-1]
        comp_nr = len(scatterer.loss_function.components)-1
        comp_type = self._returnCompType(comp)
        params = comp.__dict__
        self.spec_builder = LossEditor(self, params, comp_nr, comp_type)
        if comp_nr == 0:
            self.rePlotFig(2, rescale = False)
        else:
            self.rePlotFig(2, rescale = False)
        self.view.tables[1].fillTable()
    
    def addSynthSpec(self, start, stop, step): 
        """ This is called when the user presses the button to build their
        own spectrum. The Controller tells the Model to add a new 
        SyntheticSpectrum to the loaded spectra. Then the model tells the 
        View to open up a SpecBuilder window.
        Returned parameters is a diciontary with
        'spec_idx'-> an integer.
        'columns'-> a list of dictionaries containing 'name' and 'column width'
        'entry_fields'-> a list of dictionaries containing keys 'name', 'label'
        and 'value'
        """
        params = self.model.addSynthSpec(start, stop, step)
        return params

    def addSpecComponent(self, spec_idx, peak_kind):
        ''' This is called from the SpecBuilder window. It tells the Model to
        add a new component to the spectrum with id = spec_idx.
        '''
        self.model.addSpecComponent(spec_idx, peak_kind)
        self.rePlotFig(1, rescale=False)
        
    def updateComponent(self, new_values):
        """This is called from the SpecBuilder when the user changes any value
        of the component parameters.
        """
        self.model.updateComponent(new_values)
        self.rePlotFig(1, rescale=False)
        
    def removeSpecComp(self, spec_idx, comp_idx):
        self.model.removeSpecComp(spec_idx, comp_idx)
        self.rePlotFig(1, rescale=False)
        
    def getComponentValues(self, spec_idx, comp_idx):
        """ This is called from the SpecBuilder when the user changes the
        selection of the component. The Controller gets the values of the 
        selected component and passes them back to the SpecBuilder.
        """
        comp = self.model._returnComponentValues(spec_idx, comp_idx)
        return comp
        
    def getCompsTableData(self, spec_idx):
        """ This function retrieves the data used in the components table for
        all of the components in the synthetic spectrum.
        """
        return self.model.specBuilderTableData(spec_idx)
            
    def editRange(self, idx, start, stop, step):
        self.model.loaded_spectra[idx].start = start
        self.model.loaded_spectra[idx].stop = stop
        self.model.loaded_spectra[idx].step = step
        self.model.loaded_spectra[idx].reBuild()
        self.rePlotFig(1, rescale=False)

    def newScatterer(self, name):
        choices = self.model.newScatterer(name)
        self.view.updateScattererChoices(choices)
        
    def export(self, file):
        self.model.exportToVamas(file)      
        #self.model.exportToExcel(file)

if __name__ == "__main__":
    app = Controller()
    app.root.mainloop()
