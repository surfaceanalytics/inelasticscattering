# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 23:55:49 2020

@author: Mark
"""
import tkinter
import tkinter.ttk as tk
from tkinter import DoubleVar, StringVar, LEFT, TOP, W, Toplevel, Menu

class SpecBuilder:
    """
    The SpecBuilder is the popup window that is used to build a synthetic
    spectrum for an input specterum into the scattering algorithm.

    Parameters
    ----------
    controller : Controller
        A controller object.
    params : DICT
        The params input requires keys for:
            'columns': A list of dictionaries with keys 'name' and 'width'
            'spec_idx': The index of the spectrum that is being built
            'entry_fields': A list of dictionaries with keys 'name' and 'label'

    Returns
    -------
    None.

    """
    
    def __init__(self, controller, params):

        self.controller = controller
        self.columns = params['columns']
        self.spec_idx = params['spec_idx']
        
        self.selected_comp = None
        self.bcolor = 'grey'
        self.bthickness = 0
        self.entry_field_width = 13
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
        
        """ The 'entry_fields' holds the names and label of the user
        entry input fields. The entry fields depend on the component type.
        The variable 'entry_fields' must be a list of dictionaries.
        The keys for the dictionary must be 'name' and 'label'.
        """
        self.entries_frame = tk.Frame(self.window, padding=10)
        self.entry_fields = []
        self._buildEntryFields(params['entry_fields'])
            
        self._createWidgets()
        self._setupLayout()
        self._refreshTable()
        if len(self.peak_table.get_children()) > 0:
            self._setSelection(0)
                
    def _createWidgets(self):
        """ This creates the table that displays a list of peaks (as well as 
        the peaks' main parameters).
        """
        self._buildTable()
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
        
    def _buildTable(self):
        col = [c['name'] for c in self.columns]
        self.peak_table = tk.Treeview(self.window, height=8,
                                      show='headings',columns=col, 
                                      selectmode='browse')
        self.peak_table.name = 'spectra'
        for col in self.columns:
            n = col['name']
            w = col['width']
            self.peak_table.column(n, width=w,anchor=W)
            self.peak_table.heading(n, text=n, anchor=W)
        
        self.peak_table.bind('<Delete>', self._removeComponent)
        
    def _buildEntryFields(self, entry_fields):
        """ This is run only when the pop-up window is first created.
        """
        Width = self.entry_field_width
        callback = self._modComponent
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
        """ This is run every time a new component is created or when a
        different component is selected. It is almost the same as the 
        _buildEntryFields method, but it destroys the previously built widget,
        epties the self.entry_fields attribute, and rebuilts the frame with the
        values of the newly selected component.
        """
        self.entry_fields = []
        for child in self.entries_frame.winfo_children():
            child.destroy()
        Width = self.entry_field_width
        callback = self._modComponent
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
        """ This creates a frame that holds the entry fields for a given peak 
        type.
        """
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
              
    def _removeComponent(self, event):
        spec_idx = self.spec_idx
        comp_idx = self.selected_comp
        self.controller.removeSpecComp(spec_idx, comp_idx)
        self._refreshTable()
    
    def _addComponent(self, event):
        """ This method is bound to the 'add somponent' button. It creates a
        small pop-up menu, with a selection of component types. Once a choice
        is made, the controller is called and instructs the model to add a 
        component.
        """
        def setChoice(choice):
            self.controller.addSpecComponent(self.spec_idx, choice)
            self._refreshTable()
            self._setSelection(len(self.peak_table.get_children())-1)
        popup = Menu(self.window, tearoff=0)
        choices = self.controller.peak_choices
        for i,j in enumerate(choices):
                ''' This line below is a bit tricky. Needs to be done this way 
                because the i in the loop is only scoped for the loop, and does 
                not persist'''
                popup.add_command(command = lambda choice = choices[i]: 
                    setChoice(choice), label=j)
        popup.post(event.x_root, event.y_root)
        
    def _refreshTable(self):
        """ This function fetches all the table data from the controller,
        destroys the old table view, and builds a new table with the updated
        data.
        """
        for row in self.peak_table.get_children():
            self.peak_table.delete(row)
        table_data = self.controller.getCompsTableData(self.spec_idx)
        if len(table_data) == 0:
            self._deactivateFields()
        else:
            self._activateFields()
            
        for row_idx in table_data:
            self.peak_table.insert('',
                                   row_idx,
                                   values=table_data[row_idx], 
                                   iid=str(row_idx))
    def _deactivateFields(self):
        for field in self.entry_fields:
            field['entry'].config(state='disabled')
            
    def _activateFields(self):
        for field in self.entry_fields:
            field['entry'].config(state='enabled')
            
    def _selectComponent(self, event):
        """ This is run every time the users selects a row in the components
        table. It fetches the fields for the component, then re-builds the
        components entry fields.
        """
        sel = self.peak_table.selection()
        self.selected_comp = sel[0]
        cur_item = self.peak_table.item(sel[0])['values']
        comp_idx = cur_item[0]
        values = self.controller.getComponentValues(self.spec_idx, comp_idx)        
        self._refreshEntryFields(values)
            
    def _modComponent(self, *args):
        """ This is run every time the user changes a value in the entry 
        fields. It updates the data for the respective component in the Model
        module.
        """
        if self.window.focus_get().__class__.__name__ == 'Entry':
            #time.sleep(0.1)
            new_values = {'spec_idx':self.spec_idx,
                          'comp_idx':self.selected_comp}
            #print(len(self.entry_fields))
            for field in self.entry_fields:
                #print(field.keys())
                v = field['stringvar'].get()
                if bool(v):
                    new_values[field['name']] = v
            #print(new_values)
            self.controller.updateComponent(new_values)
            self._refreshTable()
        else:
            return
        
    def editRange(self):
        self.range_editor = RangeEditor(self.controller, self.spec_idx)
        start, stop, step = self.controller.getSpecRange(self.spec_idx)
        self.range_editor.setParams(start, stop, step)
        
    def _setSelection(self, idx):
        """This is run after the user creates a new component, to set the 
        selection to the just-added component.
        """
        self.selected_comp = idx
        self.peak_table.selection_set(str(idx))
        comp_idx = idx
        values = self.controller.getComponentValues(self.spec_idx, comp_idx)        
        self._refreshEntryFields(values)
    
class RangeEditor:
    """ This is a pop-up window what is called when the user wants to change
    the range of the synthetic spectrum.
    """
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
        start, stop, step = self.controller.getSpecRange(self.spec_idx)
        self._createWidgets()
        self._setupLayout()
        self.setParams(start, stop, step)
        
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
        