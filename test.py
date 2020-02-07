# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:54:21 2020

@author: Mark
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import tkinter as tk
import tkinter.ttk as ttk

class Spectrum:
    def __init__(self,start,stop,step):
        self.start = start
        self.stop = stop
        self.step = step
        self.x = np.arange(self.start,self.stop,self.step)
        self.clearLineshape()
        
    def clearLineshape(self):
        self.lineshape = np.zeros(int((self.stop-self.start)/self.step))
        self.x = np.arange(self.start,self.stop,self.step)

class Peak:
    def __init__(self, mean, stdev, intensity):
        self.mean = mean
        self.stdev = stdev
        self.intensity = intensity

    def Gauss(self, x):
        g = self.intensity / (self.stdev * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x-self.mean)/self.stdev)**2)
        return g
        
    def function(self, x):
        f = self.Gauss(x)
        return f
    
class SyntheticSpectrum(Spectrum):
    def __init__(self,start,stop,step):
        Spectrum.__init__(self,start,stop,step)
        self.components = []
    
    def buildLine(self):
        self.clearLineshape()
        for component in self.components:
            y = np.array([component.function(x) for x in self.x])
            self.lineshape = np.add(self.lineshape,y)
            
    def normalize(self):
        self.lineshape = self.lineshape / np.sum(self.lineshape) 
        
class LossFunction(SyntheticSpectrum):
    def __init__(self,start,stop,step):
        SyntheticSpectrum.__init__(self,start,stop,step)
        
    def addPeak(self,mean,stdev,intensity):
        self.components += [Peak(mean,stdev,intensity)]
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1     
            
    def addVacuumExcitation(self, gauss_mean, gauss_fwhm, fermi_edge, fermi_width, intensity):
        self.components += [VacuumExcitation(gauss_mean, gauss_fwhm, fermi_edge, fermi_width, intensity)]
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1  
        
    def reBuild(self):
        self.buildLine()
        self.normalize() # normalize loss function to have total area of 1  

class VacuumExcitation():
    def __init__(self, gauss_mean, gauss_fwhm, fermi_edge, fermi_width, intensity):
        self.gauss_mean = gauss_mean
        self.gauss_fwhm = gauss_fwhm
        self.fermi_edge = fermi_edge
        self.fermi_width = fermi_width
        self.intensity = intensity   
        
    def Fermi(self, x):
        k = 0.1
        f = 1/(np.exp((x-self.fermi_edge)/(k*self.fermi_width))+1)
        return f
    
    def Gauss(self, x):
        g = 1 / (self.gauss_fwhm * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x-self.gauss_mean)/self.gauss_fwhm)**2)
        return g
    
    def function(self,x):
        f = (1-self.Fermi(x)) * self.Gauss(x) * self.intensity
        return f 

class Scatterer():
    def __init__(self):
        self.label = 'default'
        self.loss_function = LossFunction(0,200,0.1)
        self.cross_sec = 0.01
        self.gas_diameter = 0.2 #In nanometers
        self.gas_cross_section = np.pi * (self.gas_diameter / 2)**2
        
    def setCrossSec(self, cross_sec):
        self.cross_sec = cross_sec
        self.loss_function.buildLine()
        self.loss_function.normalize()


#%%
class View():
    def __init__(self, controller, root):
        self.controller = controller
        self.root = root
        self.txtvar = tk.StringVar()
        self.txtvar.set('Hello')
        self.label = tk.Label(self.root, text=self.txtvar.get())
        self.label.pack(side=tk.TOP)
        self.buildTable()
        self.comp_types = ['Peak','VacuumExcitation','none']
        self.visibility_choices = ['show','hide']
        # Table
    
    def buildTable(self):
        columns = ('Nr.','Type')
        col_width = 90
        self.table = ttk.Treeview(self.root, height=5,show='headings',columns=columns,selectmode='browse')
        self.table.column('Nr.',width=col_width,anchor='center')
        self.table.column('Type',width=col_width,anchor='center')
        self.table.pack(side=tk.TOP, anchor="ne")
        self.table.bind('<<TreeviewSelect>>', self.controller.tableSelection)
        self.table.bind('<Button-3>',self.identify)
        
    def identify(self, event):
        row = self.table.identify_row(event.y)
        col = self.table.identify_column(event.x)
        col = int(col.lstrip('#'))-1
        item = self.table.item(row)
        #print(row)
        #print(col)
        def setChoice(choice):
            self.table.set(row,column=col,value=choice)
        popup = tk.Menu(self.root, tearoff=0)
        print(item['values'][col])
        
        if col == 1:
            for i,j in enumerate(self.comp_types):
                # This line below is crazy. Needs to be done this way because the i in the loop is only scoped for the loop, and does not persist
                popup.add_command(command = lambda choice = self.comp_types[i]: setChoice(choice), label=j)
        elif col == 0:
            for i,j in enumerate(self.visibility_choices):
                popup.add_command(command = lambda choice = self.visibility_choices[i]: setChoice(choice), label=j)
        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
        # make sure to release the grab (Tk 8.0a1 only)
            popup.grab_release()
        
class Controller():
    def __init__(self):

        self.loss = LossFunction(0,200,1)  
        self.loss.addPeak(1,2,3)
        self.loss.addPeak(4,5,6)
        self.loss.addVacuumExcitation(1,2,3,4,5)
        self.root = tk.Tk()
        self.view = View(self, self.root)
        self.fillTable()
        
    def fillTable(self):
        components = self.loss.components
        for i,j in enumerate(components):
            values = (i, type(j).__name__)
            self.view.table.insert('',i,values=values)
            
    def tableSelection(self, event):
        cur_item = self.view.table.selection()
        txt = self.view.table.item(cur_item[0])['values']
        #print(txt)
        self.view.txtvar.set(cur_item)
        self.view.label.configure(text=str(str(txt[0]) + ": " + txt[1]))
            


window = Controller() 
window.root.mainloop()
