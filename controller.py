# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:21:27 2020

@author: Mark
"""
from model import Spectrum, ScatteringMedium, MeasuredSpectrum, Scatterer, LossFunction, Peak, VacuumExcitation
from view import View, SpectrumBuilder
import numpy as np
import pandas as pd
import os
import tkinter as tk
import pickle 
from tkinter import StringVar
import functools
from matplotlib.widgets import RectangleSelector

class Controller():
    def __init__(self, root):
        self.start = 1200
        self.stop = 1400
        self.step = 0.1
        
        self.loaded_spectra = []
        self.unscattered_spectrum = Spectrum(self.start,self.stop,self.step)
        self.scattered_spectrum = Spectrum(self.start,self.stop,self.step)
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step)
        self.scattering_medium = ScatteringMedium()
        self.fig2_list = {'loss function':[]}
        self.intermediate_spectra =[]
        self.scattered_spectra = []
        self.bulk_spectrum = []
        self.datapath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\data'
        self.scatterers = self.read_scatterers()
        self.selected_scatterer = StringVar()
        
        self.press = None
        self.moved_while_pressed = 0
        
        self.spectra_table_choices = {1:['none','Scattered','Unscattered'],2:['show','hide']}
        self.table1_entries = []
        
        self.view = View(self, root)

    def read_scatterers(self):
        file = self.datapath+"\\scatterers"
        infile = open(file,'rb')
        scatterers = pickle.load(infile)
        infile.close()
        return scatterers
        
    def mainloop(self):
        self.root.mainloop()
        
    def loadSpectrum(self, filename):
        self.loaded_spectra += [MeasuredSpectrum(filename)]
        self.start = self.loaded_spectra[-1].start    # The step width must be defined by the measured spectrum 
        self.stop = self.loaded_spectra[-1].stop      # All synthetic pectra need to have their step widths redefined
        self.step = self.loaded_spectra[-1].step      # and their lineshapes rebuilt
        self.scattering_medium.scatterer.loss_function.step = self.step # Redefine step width of loss function
        self.simulated_spectrum = Spectrum(self.start,self.stop,self.step) # Overwrite old output spectrum with new settings
        self.rePlotFig1()
        self.insertTable1(self.loaded_spectra.index(self.loaded_spectra[-1]))
        
    def setScatteredSpectrum(self, idx):
        self.scattered_spectrum = idx
        
    def setUnscatteredSpectrum(self, idx):
        self.unscattered_spectrum = idx
        
    def show(self, idx):
        self.loaded_spectra[int(idx)].visibility = 'show'
        
    def hide(self, idx):
        self.loaded_spectra[int(idx)].visbility = 'hide'
        '''
    def onPress(self, event):
        if event.dblclick:
            self.press = None
            self.view.fig1.ax.relim()
            self.view.fig1.ax.autoscale()
            self.view.chart1.draw()

    def onMove(self, event):
        if self.press is None: return
        self.moved_while_pressed = 1
        
    def onRelease(self, event):
        if self.moved_while_pressed == 1:
            self.press['x2'] = event.xdata
            self.press['y2'] = event.ydata
            self.zoomIn()
            self.moved_while_pressed = 0
        self.press = None

    def zoomIn(self):
        xvals = [self.press['x1'],self.press['x2']]
        yvals = [self.press['y1'],self.press['y2']]
        self.view.fig1.ax.set_xlim(min(xvals),max(xvals))
        self.view.fig1.ax.set_ylim(min(yvals),max(yvals))
        self.view.chart1.draw()
        '''

    def doubleClkChart(self, event, ax, chart): 
        if event.dblclick:
            self.zoomOut(ax=ax, chart=chart)
        
    def zoomOut(self, ax, chart):
        self.press = None
        ax.relim()
        ax.autoscale()
        chart.draw()
    
    def selector(self, eclick, erelease, ax, chart):
        if eclick.dblclick: # in the case that the user double clicks, we dont want to use the selector. We use zoomOut
            pass
        else:
            xvals = [eclick.xdata,erelease.xdata]
            yvals = [eclick.ydata,erelease.ydata]
            ax.set_xlim(min(xvals),max(xvals))
            ax.set_ylim(min(yvals),max(yvals))
            chart.draw()
    
    def rePlotFig1(self):
        self.view.fig1.ax.clear()
        for i in self.loaded_spectra:
            if i.visibility == 'show':
                self.view.fig1.ax.plot(i.x,i.lineshape)
            self.view.chart1.draw()

    def reFreshFig1(self):
        left,right = self.view.fig1.ax.get_xlim()
        bottom, top = self.view.fig1.ax.get_ylim()
        self.view.fig1.ax.clear()
        for i in self.loaded_spectra:
            if i.visibility == 'show':
                self.view.fig1.ax.plot(i.x,i.lineshape)
                self.view.fig1.ax.set_xlim(left,right)
                self.view.fig1.ax.set_ylim(bottom,top)
            self.view.chart1.draw()
    
    def rePlotFig2(self):
        self.view.fig2.ax.clear()
        for i in self.fig2_list.values():
            self.view.fig2.ax.plot(i[0],i[1])
        self.view.chart2.draw()

    def reFreshFig2(self):
        left,right = self.view.fig2.ax.get_xlim()
        bottom, top = self.view.fig2.ax.get_ylim()
        self.view.fig2.ax.clear()
        for i in self.fig2_list.values():
            self.view.fig2.ax.plot(i[0],i[1])
            self.view.fig2.ax.set_xlim(left,right)
            self.view.fig2.ax.set_ylim(bottom,top)
        self.view.chart2.draw()

    def insertTable1(self,idx):
        values = (idx, self.loaded_spectra[idx].kind, self.loaded_spectra[idx].visibility)
        self.view.spectra_table.insert('',idx,values=values)
        
    def spectrumTableLogic(self, table, table_choices, row, col, choice):
        # only one spectrum can be of kind 'scattered' and one of kind 'unscattered'
        row_int = int(row.lstrip('I'))-1
        if (choice in['Scattered','Unscattered']) & (col == 1): 
            if any(x.kind == choice for x in self.loaded_spectra):
                indexes = [self.loaded_spectra.index(x) for x in self.loaded_spectra if x.kind == choice]
                for i in indexes:
                    self.loaded_spectra[i].kind = 'none'
                    if i < 9:                        
                        table.set(str('I00'+str(i+1)),column=col,value='none')
                    elif (i >= 9) & (i< 99):
                        table.set(str('I0'+str(i+1)),column=col,value='none')
                    elif (i >= 99) & (i< 99):
                        table.set(str('I'+str(i+1)),column=col,value='none')
            self.loaded_spectra[row_int].kind = choice
        for idx in self.loaded_spectra:
            if idx.kind == 'Scattered':
                self.setScatteredSpectrum(idx)
            elif idx.kind == 'Unscattered':
                self.setUnscatteredSpectrum(idx)

    def tablePopup(self, event, table, table_choices):
        row = table.identify_row(event.y)
        col = table.identify_column(event.x)
        col = int(col.lstrip('#'))-1
        def setChoice(choice):
            table.set(row,column=col,value=choice)
            if (table.name == 'spectra'):
                self.spectrumTableLogic(table, table_choices, row, col, choice)
        popup = tk.Menu(self.view.container, tearoff=0)
        if col > 0:
            choices = table_choices[col]
            for i,j in enumerate(choices):
                # This line below is a bit tricky. Needs to be done this way because the i in the loop is only scoped for the loop, and does not persist
                popup.add_command(command = lambda choice = choices[i]: setChoice(choice), label=j)
        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
            popup.grab_release()
        
    def setCurrentScatterer(self, event):
        label = self.selected_scatterer.get()
        self.scattering_medium.scatterer = self.scatterers[label]
        self.scattering_medium.scatterer.loss_function.step = self.step
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()
        self.fig2_list['loss function'] = [self.scattering_medium.scatterer.loss_function.x,self.scattering_medium.scatterer.loss_function.lineshape]
        self.rePlotFig2()
        self.fillTable2()
        self.view.cross_section.set(self.scattering_medium.scatterer.cross_sec)
        self.view.gas_diameter.set(self.scattering_medium.scatterer.gas_diameter)
    
    def fillTable2(self):
        components = self.scattering_medium.scatterer.loss_function.components
        for i,j in enumerate(components):
            values = (i, type(j).__name__)
            self.view.scatterers_table.insert('',i,values=values)
            
    def tableSelection(self,event):
        sel = self.view.scatterers_table.selection()
        cur_item = self.view.scatterers_table.item(sel[0])['values'][0]
        comp = self.scattering_medium.scatterer.loss_function.components[cur_item]
        params = comp.__dict__
        self.view.reloadEntryBox(params)
        #print(params)        
        
    def callSpectrumBuilder(self, event):
        sel = self.view.scatterers_table.selection()
        cur_item = self.view.scatterers_table.item(sel[0])['values'][0]
        comp = self.scattering_medium.scatterer.loss_function.components[cur_item]
        comp_nr = cur_item
        params = comp.__dict__
        self.spec_builder = SpectrumBuilder(self, params, comp_nr)
        
    def refreshLossFunction(self):
        comp_nr = self.spec_builder.comp_nr
        params = self.spec_builder.params
        for i in params:
            setattr(self.scattering_medium.scatterer.loss_function.components[comp_nr], i, params[i])
        self.scattering_medium.scatterer.loss_function.buildLine()
        self.scattering_medium.scatterer.loss_function.normalize()
        self.fig2_list['loss function']=[]
        self.fig2_list['loss function']=[self.scattering_medium.scatterer.loss_function.x, self.scattering_medium.scatterer.loss_function.lineshape]
        self.reFreshFig2()
        
    def doubleClkTable(self, event, table):
        item = table.focus()
        cur_item = table.item(item)['values']
        if self.loaded_spectra[cur_item[0]].visibility == 'show':
            self.loaded_spectra[cur_item[0]].visibility = 'hide'
        else:
            self.loaded_spectra[cur_item[0]].visibility = 'show'
        self.reFreshFig1()
                
    def scatterSpectrum(self):
        if ('Unscattered' not in [x.kind for x in self.loaded_spectra]) | (len(self.scattering_medium.scatterer.loss_function.components)==0):
            self.view.noScatterer()
        else:
            self.scattering_medium.setPressure(self.view.pressure.get())
            self.scattering_medium.setDistance(self.view.distance.get() * 1000000) 
            n = self.scattering_medium.n_iter
            a = self.unscattered_spectrum.lineshape
            p = self.scattering_medium.elast_prob * self.scattering_medium.scatterer.cross_sec
            self.scattering_medium.scatterer.loss_function.buildLine()
            self.scattering_medium.scatterer.loss_function.normalize()
            b = p * self.scattering_medium.scatterer.loss_function.lineshape
            self.intermediate_spectra = [a]
            for i in range(n):
                c = np.convolve(a,np.flip(b))
                l = len(a)
                c = c[-l:]
                a = (1-p) * a
                a = np.add(c,a)
                self.intermediate_spectra += [a]
            self.simulated_spectrum.lineshape = a
            self.simulated_spectrum.x = self.unscattered_spectrum.x 
            self.simulated_spectrum.kind = 'Simulated'
            if not ('Simulated' in [i.kind for i in self.loaded_spectra]):
                self.loaded_spectra += [self.simulated_spectrum]
                self.insertTable1(self.loaded_spectra.index(self.loaded_spectra[-1]))
            else:
                idx = [i.kind for i in self.loaded_spectra].index('Simulated')
                self.loaded_spectra[idx] = self.simulated_spectrum
            self.reFreshFig1()

    def updateCrossSection(self, event, *args):
        cross_sec = self.view.cross_section.get()
        if len(cross_sec) != 0:
            self.scattering_medium.scatterer.setCrossSec(float(cross_sec))
            
    def updateDiameter(self, event, *args):
        diameter = self.view.gas_diameter.get()
        if len(diameter) != 0:
            self.scattering_medium.scatterer.gas_diameter(float(diameter))
            self.scattering_medium.calcParams()            


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
            
    def writeScatterer():
        scatterers = {'default':default, 'He':He, 'N2':N2, 'O2':O2}
        file = datapath+"\\scatterers"
        outfile = open(file,'wb')
        pickle.dump(scatterers,outfile)
        outfile.close()


if __name__ == "__main__":
    mainwin = tk.Tk()
    app = Controller(mainwin)
    mainwin.mainloop()
