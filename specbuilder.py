# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 23:55:49 2020

@author: Mark
"""
import tkinter
import tkinter.ttk as tk
from tkinter import DoubleVar, StringVar, IntVar, LEFT, TOP, W, Y, N, S, Toplevel, Menu, CENTER

import time

class SpecBuilder:
    def __init__(self, controller, spec_idx, **kwargs):
        self.controller = controller
        self.columns = kwargs['columns']
        self.spec_idx = spec_idx
        self.selected_peak = None
        self.bcolor = 'grey'
        self.bthickness = 0
        padding = 15
        self.window = Toplevel(padx=padding, pady=padding)
        title = 'Component: '
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        header = tk.Label(self.window, text = 'Component: ')
        header.pack(side=TOP)
        
        self.labels = []
        self.entries = []
        self.stringvars = []
        
        self.entries_frame = tk.Frame(self.window, padding=10)
        
        if 'entry_fields' in kwargs.keys():
            ''' The variable 'entry_fields' must be a list of dictionaries.
            The keys for the dictionary must be 'name' and 'label'.
            '''
            self.entry_fields = []
            self._buildEntryFields(kwargs['entry_fields'])
            
        self._createWidgets()
        self._setupLayout()
        self.refreshTable()
        
        self.start = 0
        self.stop = 0
        self.step = 0
                
    def _createWidgets(self):
        
        col = [c['name'] for c in self.columns]
        self.peak_table = tk.Treeview(self.window, height=8,
                                      show='headings',columns=col, 
                                      selectmode='browse')
        self.peak_table.name = 'spectra'
        
        ''' This creates the table that displays a list of peaks (as well as 
        the peaks' main parameters).
        '''
        for col in self.columns:
            n = col['name']
            w = col['width']
            self.peak_table.column(n, width=w,anchor=W)
            self.peak_table.heading(n, text=n, anchor=W)
        
        self.peak_table.bind('<Delete>', self.removeComponent)
        self.add_comp_btn = tk.Label(self.window, 
                                     text = "+ Add component", 
                                     relief ='raised')
        self.add_comp_btn.bind('<Button-1>', self._addComponent)
        self.peak_table.bind('<ButtonRelease-1>',self._selectComponent)        
        self.buttons_frame = tk.Frame(self.window)
        self.edit_range_label = tk.Label(self.buttons_frame, text = '  ')
        self.edit_range_button = tk.Button(self.buttons_frame, 
                                           text='Edit range',
                                           command=self.editRange)
        self.done = tk.Button(self.buttons_frame, text='Done', command=self.Done)
        
    def _buildEntryFields(self, entry_fields):
        ''' This is run only when the pop-up window is first created.
        '''
        Width = 10
        callback = self.modPeak
        for field in entry_fields:
            d = {'name':field['name'],
                 'frame':tk.Frame(self.entries_frame),
                 'stringvar':StringVar()}
            d['stringvar'].set(field['value'])
            d['stringvar'].trace('w',callback)
            d['label'] = tk.Label(d['frame'], text = field['label'])
            d['entry'] = tk.Entry(d['frame'], width=Width, 
                             textvariable = d['stringvar'])
            self.entry_fields += [d]

    def _refreshEntryFields(self, entry_fields):
        ''' This is run every time a new component is created or when a
        different component is selected. It is almost the same as the 
        _buildEntryFields method, but it destroys the previously built widget,
        epties the self.entry_fields attribute, and rebuilts the frame with the
        values of the newly selected component.
        '''
        self.entry_fields = []
        for child in self.entries_frame.winfo_children():
            child.destroy()
        Width = 10
        callback = self.modPeak
        for field in entry_fields:
            d = {'name':field['name'],
                 'frame':tk.Frame(self.entries_frame),
                 'stringvar':StringVar()}
            d['stringvar'].set(field['value'])
            d['stringvar'].trace('w',callback)
            d['label'] = tk.Label(d['frame'], text = field['label'])
            d['entry'] = tk.Entry(d['frame'], width=Width, 
                             textvariable = d['stringvar'])
            d['frame'].pack(side=TOP)
            d['label'].pack(side=TOP)
            d['entry'].pack(side=TOP)
            self.entry_fields += [d]
        
    def _setupLayout(self):
        
        ''' This creates a frame that holds the entry fields for a given peak 
        type.
        '''
        self.entries_frame.pack(side=LEFT)
        if hasattr(self, 'entry_fields'):
            for field in self.entry_fields:
                field['frame'].pack(side=TOP)
                field['label'].pack(side=TOP)
                field['entry'].pack(side=TOP)
        
        self.peak_table.pack(side=TOP)
        self.add_comp_btn.pack(side=TOP,anchor='ne')
        self.buttons_frame.pack(side=TOP)
        self.edit_range_label.pack(side=LEFT)
        self.edit_range_button.pack(side=LEFT)        
        self.done.pack(side=LEFT)
        
    def Done(self):
        self.window.destroy()
              
    def removeComponent(self):
        return
    
    def _addComponent(self, event):
        def setChoice(choice):
            self.controller.addPeak(self.spec_idx, choice)
            self.refreshTable()
            self.controller.rePlotFig1()
            self.setSelection(len(self.peak_table.get_children())-1)
        popup = Menu(self.window, tearoff=0)
        choices = self.controller.peak_choices
        for i,j in enumerate(choices):
                ''' This line below is a bit tricky. Needs to be done this way 
                because the i in the loop is only scoped for the loop, and does 
                not persist'''
                popup.add_command(command = lambda choice = choices[i]: 
                    setChoice(choice), label=j)
        popup.post(event.x_root, event.y_root)
        
    def refreshTable(self):
        for row in self.peak_table.get_children():
            self.peak_table.delete(row)
        values = self.controller.getComps(self.spec_idx)
        for value in values:
            self.peak_table.insert('',
                                   value,
                                   values=values[value], 
                                   iid=str(value))
            
    def _selectComponent(self, event):
        ''' This is run every time the users selects a row in the components
        table. It fetches the fields for the component, then re-builds the
        components entry fields.
        '''
        sel = self.peak_table.selection()
        self.selected_peak = sel[0]
        cur_item = self.peak_table.item(sel[0])['values']
        comp_idx = cur_item[0]
        values = self.controller.getComponentValues(self.spec_idx, comp_idx)        
        print(values)
        self._refreshEntryFields(values)
            
    def modPeak(self, *args):
        ''' This is run every time the user changes a value in the entry 
        fields. It updates the data for the respective component in the Model
        module.
        '''
        if self.window.focus_get().__class__.__name__ == 'Entry':
            time.sleep(0.1)
            new_values = {'spec_idx':self.spec_idx,
                          'comp_idx':self.selected_peak}
            #print(len(self.entry_fields))
            
            for field in self.entry_fields:
                #print(field.keys())
                v = field['stringvar'].get()
                if bool(v):
                    new_values[field['name']] = v
            print(new_values)
            self.controller.updatePeak(new_values)
            self.refreshTable()
            #self.updateTable(new_values)
        else:
            return
            
    def updateTable(self, new_values):
        comp_idx = new_values['comp_idx']
        comp_type = self.peak_table.item(self.selected_peak)['values'][1]
        values = (comp_idx, comp_type)
        for key, value in new_values.items():
            if (key != 'spec_idx') and (key != 'comp_idx'):
                values += tuple(value)
        self.peak_table.item(comp_idx, values = values)
        print(values)
        
    def editRange(self):
        self.range_editor = RangeEditor(self.controller, self.spec_idx)
        self.range_editor.setParams(self.start, self.stop, self.step)
        
    def setSelection(self, idx):  
        self.selected_peak = idx
        self.peak_table.selection_set(str(idx))
        cur_item = self.peak_table.item(idx)['values']
        comp_idx = cur_item[0]
        values = self.controller.getComponentValues(self.spec_idx, comp_idx)        
        self._refreshEntryFields(values)
    
class RangeEditor:
    ''' This is a pop-up windos what is called when the user wants to change
    the range of the synthetic spectrum.
    '''
    def __init__(self, controller, spec_idx):
        self.controller = controller
        self.spec_idx = spec_idx
        self.bcolor = 'grey'
        self.bthickness = 0
        padding = 15
        self.window = Toplevel(padx=padding, pady=padding)
        title = 'Range: '
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        header = tk.Label(self.window, text = 'Enter range for spectrum')
        header.pack(side=TOP)
        self._createWidgets()
        self._setupLayout()
        
    def _createWidgets(self):
        self.container = tk.Frame(self.window)
        self.f1 = tk.Frame(self.container)
        self.f2 = tk.Frame(self.container)
        self.f3 = tk.Frame(self.container)
        
        self.start = DoubleVar()
        self.start_label = tk.Label(self.f1, text = 'Start')
        self.start_entry = tk.Entry(self.f1, 
                                    width = 10, 
                                    textvariable = self.start)
        self.stop = DoubleVar()
        self.stop_label = tk.Label(self.f2, text = 'Stop')
        self.stop_entry = tk.Entry(self.f2, 
                                   width = 10, 
                                   textvariable = self.stop)
        self.step = DoubleVar()
        self.step_label = tk.Label(self.f3, text = 'Step')
        self.step_entry = tk.Entry(self.f3, 
                                   width = 10, 
                                   textvariable = self.step)
        
        self.done_button = tk.Button(self.window, 
                                     text = 'Done', 
                                     command = self.Done)
        
    def _setupLayout(self):
        self.container.pack(side=TOP)
        self.f1.pack(side=LEFT)
        self.f2.pack(side=LEFT)  
        self.f3.pack(side=LEFT)
        self.start_label.pack(side=TOP)
        self.start_entry.pack(side=TOP)
        self.stop_label.pack(side=TOP)
        self.stop_entry.pack(side=TOP)
        self.step_label.pack(side=TOP)
        self.step_entry.pack(side=TOP)
        self.done_button.pack(side=TOP)
        
    def Done(self):
        start = self.start.get()
        stop = self.stop.get()
        step = self.step.get()
        idx = self.spec_idx
        self.controller.editRange(idx, start, stop, step)
        self.window.destroy()
        
    def setParams(self, start, stop, step):
        self.start.set(start)
        self.stop.set(stop)
        self.step.set(step)
        