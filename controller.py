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
        self.sendVariantChoices()
        self.default_algorithm = self.model.algorithm_id
        self.view.variant.set(str(self.default_algorithm))
        self.getInputsFromModel()

        self.img_collection = []

        self.colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
                                
    def mainloop(self):
        self.root.mainloop()
        
    def sendVariantChoices(self):
        params = list(self.model.inputs_dict.keys())
        self.view.addVariantChoices(params)
        
    def toggleVariant(self, new_var):
        ''' Set new algorithm_id in model'''
        self.model.changeAlgorithm(new_var)
        self.getInputsFromModel()
        
    def sendInputsToView(self, params):
        self.view.buildAlgorithmFrame(params)
        
    def getInputsFromModel(self):
        params = self.model.inputs
        self.sendInputsToView(params)
        
    def updateValuesInModel(self):
        ''' Get new values from view'''
        params = self.view.inputs_frame.sendValues()
        ''' Send new values to model'''
        self.model.updateAlgorithmParams(params)
                
    def scatterSpectrum(self):
        self.updateValuesInModel()
        ''' This makes sure there it is defined which spectrum should be used as 
        the input spectrum'''
        if ('Unscattered' not in [x.kind for x in self.model.loaded_spectra]) | (len(self.model.scattering_medium.scatterer.loss_function.components)==0):
            self.view.noScatterer()
        else:
            ''' This checks if the user wants to display the 'bulk' spectrum'''
            if self.view.bulk.get() == 1:
                self.model.algorithm_option = 'bulk'
            else:
                self.model.algorithm_option = 'film'
            
            self.model.scatterSpectrum()
            
            self.fillTable1()
            self.reFreshFig1()
                
    def setScatteredSpectrum(self, spectrum):
        self.model.scattered_spectrum = spectrum
        
    def setUnscatteredSpectrum(self, spectrum):
        self.model.unscattered_spectrum = spectrum
        self.model.scattering_medium.scatterer.loss_function.step = spectrum.step
        self.model.scattering_medium.scatterer.loss_function.reBuild()
        self.rePlotFig2()
        
    def show(self, idx):
        self.model.loaded_spectra[int(idx)].visibility = 'visible'
        
    def hide(self, idx):
        self.model.loaded_spectra[int(idx)].visbility = 'hidden'
            
    def loadSpectrum(self,file):
        self.model.loadSpectrum(file)
        self.rePlotFig1()
        self.insertTable1(self.model.loaded_spectra.index(self.model.loaded_spectra[-1]))
   
    def rePlotFig1(self):
        ''' This replots figure 1 and re-sets the x and y limits'''
        self.view.fig1.ax.clear()
        self.view.fig1.ax.set_xlabel('Energy [eV]', fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_ylabel('Intensity [cts./sec.]', fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_title('Spectra')
        if self.view.normalize.get() == 1:
            for indx, i in enumerate(self.model.loaded_spectra):
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape/np.max(i.lineshape),
                                           c=self.colours[indx])
        else:
            for indx, i in enumerate(self.model.loaded_spectra):
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape, c=self.colours[indx])
        self.view.fig1.fig.tight_layout()
        self.view.chart1.draw()

    def reFreshFig1(self):
        ''' This replots figure 1 and keeps the current the x and y limits'''
        left,right = self.view.fig1.ax.get_xlim()
        bottom, top = self.view.fig1.ax.get_ylim()
        self.view.fig1.ax.clear()
        self.view.fig1.ax.set_xlabel('Energy [eV]', 
                                     fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_ylabel('Intensity [cts./sec.]', 
                                     fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_title('Spectra')
        if self.view.normalize.get() == 1:
            for indx, i in enumerate(self.model.loaded_spectra):
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape/np.max(i.lineshape),
                                           c=self.colours[indx])
        else:
            for indx, i in enumerate(self.model.loaded_spectra):
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape, c=self.colours[indx])
        self.view.fig1.ax.set_xlim(left,right)
        self.view.fig1.ax.set_ylim(bottom,top)
        self.view.fig1.fig.tight_layout()
        self.view.chart1.draw()
    
    def rePlotFig2(self):
        ''' This replots figure 2 and re-sets the x and y limits'''
        self.view.fig2.ax.clear()
        self.view.fig2.ax.set_xlabel('Energy Loss [eV]', 
                                     fontsize=self.view.axis_label_fontsize)
        self.view.fig2.ax.set_ylabel('Probability', 
                                     fontsize=self.view.axis_label_fontsize)
        self.view.fig2.ax.set_title('Loss Function')
        x = self.model.scattering_medium.scatterer.loss_function.x
        y = self.model.scattering_medium.scatterer.loss_function.lineshape
        self.view.fig2.ax.plot(x,y)
        self.view.fig2.fig.tight_layout()
        self.view.chart2.draw()

    def reFreshFig2(self):
        '''This re-plots figure 2, but keeps the existing x and y limits'''
        left,right = self.view.fig2.ax.get_xlim()
        bottom, top = self.view.fig2.ax.get_ylim()
        self.view.fig2.ax.clear()
        self.view.fig2.ax.set_xlabel('Energy Loss [eV]', 
                                     fontsize=self.view.axis_label_fontsize)
        self.view.fig2.ax.set_ylabel('Probability', 
                                     fontsize=self.view.axis_label_fontsize)
        x = self.model.scattering_medium.scatterer.loss_function.x
        y = self.model.scattering_medium.scatterer.loss_function.lineshape
        self.view.fig2.ax.plot(x,y)
        self.view.fig2.ax.set_xlim(left,right)
        self.view.fig2.ax.set_ylim(bottom,top)
        self.view.fig2.fig.tight_layout()
        self.view.chart2.draw()

    def insertTable1(self,idx):
        values = (idx, self.model.loaded_spectra[idx].kind,
                  self.model.loaded_spectra[idx].visibility)
        self.img = PhotoImage(file=str(self.resourcepath) + '\\' +
                                 'legend' + str(idx) + '.png')
        self.img = self.img.subsample(2, 4)
        self.img_collection.append(self.img)
        self.view.spectra_table.insert('', idx, values=values, image=self.img, 
                                       iid=str(idx))

    def spectrumTableLogic(self, table, table_choices, row, col, choice):
        '''Only one spectrum can be of kind 'scattered' and one of kind 
        'unscattered' '''
        if choice == 'edit':
            ''' Get index of selected spectrum'''
            spec_idx = int(table.item(row)['values'][0])
            '''Pass index to Model and get parameters''' 
            params = self.model.returnSpecComps(spec_idx)
            self.view.editSynthSpec(params)
        elif (choice in['Scattered','Unscattered']) & (col == 1): 
            if any(x.kind == choice for x in self.model.loaded_spectra):
                indexes = [self.model.loaded_spectra.index(x) for x in self.model.loaded_spectra if x.kind == choice]
                for i in indexes:
                    self.model.loaded_spectra[i].kind = 'none'
                    table.set(str(i),column=col,value='none')
            self.model.loaded_spectra[int(row)].kind = choice
        elif choice == 'none':
            self.model.loaded_spectra[int(row)].kind = choice
        for spectrum in self.model.loaded_spectra:
            if spectrum.kind == 'Scattered':
                self.setScatteredSpectrum(spectrum)
            elif spectrum.kind == 'Unscattered':
                self.setUnscatteredSpectrum(spectrum)

    def getSpectrumType(self, idx):
        return self.model.loaded_spectra[idx].__class__.__name__
        
    def tablePopup(self, event, table, table_choices):
        row = table.identify_row(event.y)
        col = table.identify_column(event.x)
        col = int(col.lstrip('#'))-1
        def setChoice(choice):
            if choice != 'edit':
                table.set(row,column=col,value=choice)
            if (table.name == 'spectra'):
                self.spectrumTableLogic(table, table_choices, row, col, choice)
        popup = Menu(self.view.container, tearoff=0)
        if col > 0:
            choices = table_choices[col]
            spectrum_type = self.getSpectrumType(int(table.item(row)['values'][0]))
            if spectrum_type == 'SyntheticSpectrum':
                choices += ['edit']
            for i,j in enumerate(choices):
                # This line below is a bit tricky. Needs to be done this way because the i in the loop is only scoped for the loop, and does not persist
                popup.add_command(command = lambda choice = choices[i]: setChoice(choice), label=j)
            if 'edit' in choices:
                del choices[choices.index('edit')]

        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
            popup.grab_release()
        
    def setCurrentScatterer(self, label):
        """This is called from the View when the user selects a different 
        scatterer. The Controller retrieves the scatterer's parameters
        from the Model and sends them to the View.
        """
        self.model.setCurrentScatterer(label)
        self.getInputsFromModel()
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.rePlotFig2()
        self.fillTable2()

    def fillTable1(self):
        for row in self.view.spectra_table.get_children():
            self.view.spectra_table.delete(row)
        # refresh image list:
            self.img_collection = []
        for i in range(len(self.model.loaded_spectra)):
            self.insertTable1(i)
   
    def fillTable2(self):
        for row in self.view.scatterers_table.get_children():
            self.view.scatterers_table.delete(row)
        components = self.model.scattering_medium.scatterer.loss_function.components
        for i,j in enumerate(components):
            values = (i, type(j).__name__)
            self.view.scatterers_table.insert('',i,values=values, iid=str(i))
            
    def tableSelection(self,event):
        sel = self.view.scatterers_table.selection()
        cur_item = self.view.scatterers_table.item(sel[0])['values'][0]
        comp = self.model.scattering_medium.scatterer.loss_function.components[cur_item]
        params = comp.__dict__
        self.view.reloadEntryBox(params)
    
    def _returnCompType(self, comp):
        return str(comp.__class__).split(".")[1].strip("'>")
    
    def callLossEditor(self, event):
        """ The loss editor is a pop-up window with entry fields for the 
        parameters of the loss function component"""
        ''' The next statement gets the currently selected component from 
        the loss function'''
        sel = self.view.scatterers_table.selection()
        ''' This next gets the id associated with the currently selected
        component in the loss function'''
        comp_nr = self.view.scatterers_table.item(sel[0])['values'][0]
        ''' This gets the component from the model, using the id of the 
        component'''
        comp = self.model.scattering_medium.scatterer.loss_function.components[comp_nr]
        comp_type = self._returnCompType(comp)
        params = comp.__dict__
        
        self.spec_builder = LossEditor(self, params, comp_nr, comp_type)
        
    def changeLossFunction(self):
        comp_nr = self.spec_builder.comp_nr
        params = self.spec_builder.params
        self.model.changeLossFunction(comp_nr, params)
        self.fig2_list['loss function']=[]
        self.fig2_list['loss function']=[self.model.scattering_medium.scatterer.loss_function.x, self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.reFreshFig2()
        
    def doubleClkTable(self, event, table):
        '''This toggles the visibility of spectra in figure 1 by double 
        clicking the left mouse button.'''
        item = table.focus()
        cur_item = table.item(item)['values']
        if self.model.loaded_spectra[cur_item[0]].visibility == 'visible':
            self.model.loaded_spectra[cur_item[0]].visibility = 'hidden'
        else:
            self.model.loaded_spectra[cur_item[0]].visibility = 'visible'
        table.set(item, column=2,value=self.model.loaded_spectra[cur_item[0]].visibility)
        self.reFreshFig1()

    def updateInelasticXSect(self, event, *args):
        inelastic_xsect = self.view.inelastic_xsect.get()
        if len(inelastic_xsect) != 0:
            self.model.setInelasticXSect(float(inelastic_xsect))
            
    def updateElasticXSect(self, event, *args):
        elastic_xsect = self.view.elastic_xsect.get()
        if len(elastic_xsect) != 0:
            self.model.setElasticXSect(float(elastic_xsect))
            
    def updateScatterersDict(self):
        self.model.updateScatterersDict()

    def loadScatterers(self, file):
        """ Pass file to model to load the scatterers within, update view with
        scatter_choices"""
        # load scatterers from file
        scatterer_choices = self.model.loadScatterers(file)
        # update choice in menu
        self.view.updateScattererChoices(scatterer_choices)
        self.view.selected_scatterer.set(scatterer_choices[0])
        self.setCurrentScatterer(scatterer_choices[0])

    def saveScatterers(self, file):
        self.model.saveScatterers(file)

    def createSynthetic(self):
        self.model.loaded_spectra += [self.model.SyntheticSpectrum(self.model.start, 
                                                                   self.model.stop, 
                                                                   self.model.step)]
        self.model.loaded_spectra[-1].buildLine()
        
    def removeSpectrum(self, idx):
        """ This removes a spectrum from the loaded spectra. It is called
        from the View when the user preses <delete> on a row in the spectra
        table. The Controller tells the Model to remove a spectrum from the 
        list of loaded spectra. Then the Controller tells the view to refresh
        the spectra figure and the spectra table.
        """
        del self.model.loaded_spectra[idx]
        self.fillTable1()
        self.reFreshFig1()
        
    def removeLossComponent(self, comp_idx):
        """ This removes a component from the loss function. It is called when
        the user presses <delete> on a selected row in the loss function table.
        The Controller tells the Model to remove one component from the 
        currently loaded loss function. Then the Controller tells the View to 
        replot the loss function and refresh the loss function table.
        """
        self.model.scattering_medium.scatterer.loss_function.removeComponent(comp_idx)
        self.fillTable2()
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.reFreshFig2()
        
    def addComponent(self, comp_kind):
        """This adds a new component to the loss function. It is called when
        the user presses 'add component' button on the loss functions
        panel, and selects a type of loss function component.
        The Controller tells the Model to add one component of a chosen type
        to the currently loaded loss function. Then the Controller tells the 
        View to replot the loss function and refresh the loss function table.
        """
        self.model.addLossComponent(comp_kind)
        self.fillTable2()
        scatterer = self.model.scattering_medium.scatterer
        self.fig2_list['loss function'] = [scatterer.loss_function.x,
                      scatterer.loss_function.lineshape]
        comp = scatterer.loss_function.components[-1]
        comp_nr = len(scatterer.loss_function.components)-1
        comp_type = self._returnCompType(comp)
        params = comp.__dict__
        self.spec_builder = LossEditor(self, params, comp_nr, comp_type)
        if comp_nr == 0:
            self.rePlotFig2()
            self.view.zoomOut(self.view.fig2.ax, self.view.fig2.chart)
        else:
            self.reFreshFig2()
    
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
        self.rePlotFig1()
        
    def updateComponent(self, new_values):
        """This is called from the SpecBuilder when the user changes any value
        of the component parameters.
        """
        self.model.updateComponent(new_values)
        self.reFreshFig1()
        
    def removeSpecComp(self, spec_idx, comp_idx):
        self.model.removeSpecComp(spec_idx, comp_idx)
        self.rePlotFig1()
        
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
        self.reFreshFig1()

    def newScatterer(self, name):
        choices = self.model.newScatterer(name)
        self.view.updateScattererChoices(choices)
        
    def export(self, file):
        self.model.exportToVamas(file)      
        #self.model.exportToExcel(file)

if __name__ == "__main__":
    app = Controller()
    app.root.mainloop()
