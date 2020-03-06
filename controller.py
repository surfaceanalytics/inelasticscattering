# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:21:27 2020

@author: Mark
"""
from model import Model, SyntheticSpectrum
from view import View, LossEditor
import numpy as np
import os
import tkinter as tk
from tkinter import StringVar

class Controller():
    def __init__(self):
        self.root = tk.Tk()

        self.fig2_list = {'loss function':[]}
        
        self.scatterer_choices = []
        self.selected_scatterer = StringVar()
        
        self.press = None
        self.moved_while_pressed = 0
        
        self.spectra_table_choices = {1:['none','Scattered','Unscattered'],2:['visible','hidden']}
        self.table1_entries = []
        self.datapath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\data'
        self.scatt_datapath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\data'
        
        self.model = Model(self)
        self.model.readScatterers()
        self.component_choices = self.model.loss_component_kinds
        self.peak_choices = self.model.peak_kinds
        self.view = View(self, self.root)
        
    def mainloop(self):
        self.root.mainloop()
        
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
        # this replots figure 1 and re-sets the x and y limits
        self.view.fig1.ax.clear()
        self.view.fig1.ax.set_xlabel('Energy [eV]', fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_ylabel('Intensity [cts./sec.]', fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_title('Spectra')
        if self.view.normalize.get() == 1:
            for i in self.model.loaded_spectra:
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape/np.max(i.lineshape))
        else:
            for i in self.model.loaded_spectra:
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape)
        self.view.fig1.fig.tight_layout()
        self.view.chart1.draw()

    def reFreshFig1(self):
        # this replots figure 1 and keeps the current the x and y limits
        left,right = self.view.fig1.ax.get_xlim()
        bottom, top = self.view.fig1.ax.get_ylim()
        self.view.fig1.ax.clear()
        self.view.fig1.ax.set_xlabel('Energy [eV]', fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_ylabel('Intensity [cts./sec.]', fontsize=self.view.axis_label_fontsize)
        self.view.fig1.ax.set_title('Spectra')
        if self.view.normalize.get() == 1:
            for i in self.model.loaded_spectra:
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape/np.max(i.lineshape))
        else:
            for i in self.model.loaded_spectra:
                if i.visibility == 'visible':
                    self.view.fig1.ax.plot(i.x,i.lineshape)
        self.view.fig1.ax.set_xlim(left,right)
        self.view.fig1.ax.set_ylim(bottom,top)
        self.view.fig1.fig.tight_layout()
        self.view.chart1.draw()
    
    def rePlotFig2(self):
        # this replots figure 2 and re-sets the x and y limits
        self.view.fig2.ax.clear()
        self.view.fig2.ax.set_xlabel('Energy Loss [eV]', fontsize=self.view.axis_label_fontsize)
        self.view.fig2.ax.set_ylabel('Probability', fontsize=self.view.axis_label_fontsize)
        self.view.fig2.ax.set_title('Loss Function')
        x = self.model.scattering_medium.scatterer.loss_function.x
        y = self.model.scattering_medium.scatterer.loss_function.lineshape
        #for i in self.fig2_list.values():
        #    self.view.fig2.ax.plot(i[0],i[1])
        self.view.fig2.ax.plot(x,y)
        self.view.fig2.fig.tight_layout()
        self.view.chart2.draw()

    def reFreshFig2(self):
        # this re-plots figure 2, but keeps the existing x and y limits
        left,right = self.view.fig2.ax.get_xlim()
        bottom, top = self.view.fig2.ax.get_ylim()
        self.view.fig2.ax.clear()
        self.view.fig2.ax.set_xlabel('Energy Loss [eV]', fontsize=self.view.axis_label_fontsize)
        self.view.fig2.ax.set_ylabel('Probability', fontsize=self.view.axis_label_fontsize)
        x = self.model.scattering_medium.scatterer.loss_function.x
        y = self.model.scattering_medium.scatterer.loss_function.lineshape
        #for i in self.fig2_list.values():
        #    self.view.fig2.ax.plot(i[0],i[1])
        self.view.fig2.ax.plot(x,y)
        self.view.fig2.ax.set_xlim(left,right)
        self.view.fig2.ax.set_ylim(bottom,top)
        self.view.fig2.fig.tight_layout()
        self.view.chart2.draw()

    def insertTable1(self,idx):
        values = (idx, self.model.loaded_spectra[idx].kind, self.model.loaded_spectra[idx].visibility)
        self.view.spectra_table.insert('',idx,values=values, iid=str(idx))
        
    def spectrumTableLogic(self, table, table_choices, row, col, choice):
        # only one spectrum can be of kind 'scattered' and one of kind 'unscattered'
        #if self.getSpectrumType(int(table.item(row)['values'][0])) == 'SyntheticSpectrum':
        if choice == 'edit':
            self.view.editSynthSpec(int(table.item(row)['values'][0]))
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
        popup = tk.Menu(self.view.container, tearoff=0)
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
        
    def setCurrentScatterer(self, event):
        label = self.selected_scatterer.get()
        self.model.setCurrentScatterer(label)
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.rePlotFig2()
        self.fillTable2()
        self.view.cross_section.set(self.model.scattering_medium.scatterer.cross_section)
        self.view.gas_diameter.set(self.model.scattering_medium.scatterer.gas_diameter)
        self.view.angle_factor.set(self.model.scattering_medium.scatterer.angle_factor)
    
    def fillTable1(self):
        for row in self.view.spectra_table.get_children():
            self.view.spectra_table.delete(row)
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
        
    def callLossEditor(self, event):
        sel = self.view.scatterers_table.selection()
        cur_item = self.view.scatterers_table.item(sel[0])['values'][0]
        comp = self.model.scattering_medium.scatterer.loss_function.components[cur_item]
        comp_nr = cur_item
        params = comp.__dict__
        self.spec_builder = LossEditor(self, params, comp_nr)
        
    def changeLossFunction(self):
        comp_nr = self.spec_builder.comp_nr
        params = self.spec_builder.params
        self.model.changeLossFunction(comp_nr, params)
        self.fig2_list['loss function']=[]
        self.fig2_list['loss function']=[self.model.scattering_medium.scatterer.loss_function.x, self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.reFreshFig2()
        
    def doubleClkTable(self, event, table):
        # This toggles the visibility of spectra in figure 1 by double clicking the left mouse button
        item = table.focus()
        cur_item = table.item(item)['values']
        if self.model.loaded_spectra[cur_item[0]].visibility == 'visible':
            self.model.loaded_spectra[cur_item[0]].visibility = 'hidden'
        else:
            self.model.loaded_spectra[cur_item[0]].visibility = 'visible'
        table.set(item, column=2,value=self.model.loaded_spectra[cur_item[0]].visibility)
        self.reFreshFig1()
                
    def scatterSpectrum(self):
        # this makes sure there it is defined which spectrum should be used as the input spectrum
        if ('Unscattered' not in [x.kind for x in self.model.loaded_spectra]) | (len(self.model.scattering_medium.scatterer.loss_function.components)==0):
            self.view.noScatterer()
        else:
            # This checks if the advanced panel is active or not.
            if self.view.advanced.get() == 0:
                # If not advanced, then the pressure and distance values are used to determine the number of iterations and probability to use
                self.model.scattering_medium.setPressure(self.view.pressure.get())
                self.model.scattering_medium.setDistance(self.view.distance.get() * 1000000) # convert mm to nm
                self.model.scattering_medium.scatterer.setCrossSec(float(self.view.cross_section.get()))
                self.model.scatterSpectrum()
            elif self.view.advanced.get() == 1:
                # If advanced is chosen, then the number of iterations, probability and angle factor are taken from the advanced setting entries
                self.model.scattering_medium.n_iter = int(self.view.n_iter.get())
                self.model.scattering_medium.collis_prob = 1
                self.model.scattering_medium.scatterer.cross_section = self.view.prob.get()
                self.model.scatterSpectrum()

            if self.view.bulk.get() == 1:
                # this checks if the user wants to display the 'bulk' spectrum
                self.model.simulated_spectrum.lineshape = self.model.bulk_spectrum 
            else:
                self.model.simulated_spectrum.lineshape = self.model.intermediate_spectra[-1]

            if not ('Simulated' in [i.kind for i in self.model.loaded_spectra]):
                self.model.loaded_spectra += [self.model.simulated_spectrum]
                self.insertTable1(self.model.loaded_spectra.index(self.model.loaded_spectra[-1]))
            else:
                idx = [i.kind for i in self.model.loaded_spectra].index('Simulated')
                self.model.loaded_spectra[idx] = self.model.simulated_spectrum
            self.reFreshFig1()

    def updateCrossSection(self, event, *args):
        cross_section = self.view.cross_section.get()
        if len(cross_section) != 0:
            self.model.scattering_medium.scatterer.setCrossSec(float(cross_section))
            
    def updateDiameter(self, event, *args):
        diameter = self.view.gas_diameter.get()
        if len(diameter) != 0:
            self.model.scattering_medium.scatterer.gas_diameter=float(diameter)
            self.model.scattering_medium.calcParams()
            
    def updateAngle(self, event, *args):
        angle_factor = self.view.angle_factor.get()
        if len(angle_factor) != 0:
            self.model.scattering_medium.scatterer.angle_factor=float(angle_factor)
            self.model.scattering_medium.calcParams()
            
    def updateScatterersDict(self):
        self.model.updateScatterersDict()

    def saveScatterers(self):
        self.model.saveScatterers()
        
    def createSynthetic(self):
        self.model.loaded_spectra += [self.model.SyntheticSpectrum(self.model.start, self.model.stop, self.model.step)]
        self.model.loaded_spectra[-1].buildLine()
        
    def removeSpectrum(self, idx):
        del self.model.loaded_spectra[idx]
        self.fillTable1()
        self.reFreshFig1()
        
    def removeComponent(self, comp_idx):
        self.model.scattering_medium.scatterer.loss_function.removeComponent(comp_idx)
        self.fillTable2()
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        self.reFreshFig2()
        
    def addComponent(self, comp_kind):
        self.model.addComponent(comp_kind)
        self.fillTable2()
        self.fig2_list['loss function'] = [self.model.scattering_medium.scatterer.loss_function.x,self.model.scattering_medium.scatterer.loss_function.lineshape]
        comp = self.model.scattering_medium.scatterer.loss_function.components[-1]
        comp_nr = len(self.model.scattering_medium.scatterer.loss_function.components)-1
        params = comp.__dict__
        self.spec_builder = LossEditor(self, params, comp_nr)
        if comp_nr == 0:
            self.rePlotFig2()
            self.view.zoomOut(self.view.fig2.ax, self.view.fig2.chart)
        else:
            self.reFreshFig2()
    
    def addSynthSpec(self, start, stop, step):
        self.model.loaded_spectra += [SyntheticSpectrum(start, stop, step)]

    def addPeak(self, spec_idx, peak_kind):
        self.model.addPeak(spec_idx,peak_kind)
        
    def updatePeak(self, new_values):
        self.model.updatePeak(new_values)
        self.reFreshFig1()
        
    def getComps(self, spec_idx):
        comps = self.model.loaded_spectra[spec_idx].components
        values = {}
        for i, comp in enumerate(comps):
            values[i] = [i, comp.__class__.__name__, comp.position, comp.width, comp.intensity]
        return values
            
    def editRange(self, idx, start, stop, step):
        self.model.loaded_spectra[idx].start = start
        self.model.loaded_spectra[idx].stop = stop
        self.model.loaded_spectra[idx].step = step
        self.model.loaded_spectra[idx].reBuild()
        self.reFreshFig1()

'''    
    def writeScatterer():
        scatterers = {'default':default, 'He':He, 'N2':N2, 'O2':O2}
        file = datapath+"\\scatterers"
        outfile = open(file,'wb')
        pickle.dump(scatterers,outfile)
        outfile.close()

    def exportExcel(self, filename):
        input_x = self.unscattered_spectrum.x
        input_spec = self.unscattered_spectrum.lineshape / np.max(self.unscattered_spectrum.lineshape)
        output_spec = self.simulated_spectrum.lineshape / np.max(self.simulated_spectrum.lineshape)
        fit = self.scattered_spectrum.lineshape / np.max(self.scattered_spectrum.lineshape)
        df1 = pd.DataFrame(np.stack([input_x,input_spec,output_spec,fit]).T, columns=['kinetic energy (eV)','unscattered spectrum', 'fit spectrum','scattered spectrum'])
        loss_x = self.scattering_medium.scatterer.loss_function.x
        loss = self.scattering_medium.scatterer.loss_function.lineshape
        df2 =  pd.DataFrame(np.stack([loss_x,loss]).T,columns=['energy loss (eV)','proability'])
        file = filename + '.xlsx'
        with pd.ExcelWriter(file) as writer:
            df1.to_excel(writer,sheet_name='spectra',startcol=0)
            df2.to_excel(writer,sheet_name='spectra',startcol=5)
'''            


if __name__ == "__main__":
    app = Controller()
    app.root.mainloop()
