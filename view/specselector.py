# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 19:09:32 2020

@author: Mark
"""

import tkinter
import tkinter.ttk as tk
from tkinter import (DoubleVar, StringVar, LEFT, TOP, RIGHT, W, Y, N, S, Toplevel, 
                     Menu, CENTER)
from view.table import Table
from view.figure import Figure

class SpecSelector:
    """ The SpecSelector is the popup window that is used to select spectra to
    load from a Vamas file.
    """
    def __init__(self, controller, table_params, fig_params):
        self.controller = controller
        self.table_params = table_params
        self.fig_params = fig_params
        self._createWidgets()
        self._setupLayout()
        self._setBindings()
        self.table.fillTable()
        
    def _createWidgets(self):
        padding = 15
        self.window = Toplevel(padx=padding, pady=padding)
        self.window.title('Spectrum Selector')
        self.container = tk.Frame(self.window)
        self.table_frame = tk.Frame(self.container)
        self.figure_frame = tk.Frame(self.container)
        self.table = Table(self.window,
                           self.controller, 
                           self.table_frame, 
                           self.table_params)
        x = []
        y = []
        self.figure = Figure(self.figure_frame, x, y, self.fig_params)
        self.done = tk.Button(self.figure_frame, 
                              text = 'Done', 
                              command = self._Done)
        
    def _buildEntryFields(self, entry_fields):
        pass

    def _refreshEntryFields(self, entry_fields):
        pass
        
    def _setupLayout(self):        
        """ This creates a frame that holds the entry fields for a given peak 
        type.
        """
        self.container.pack(side=TOP)
        self.table_frame.pack(side=LEFT)
        self.table.table.pack(side=TOP, fill='y') 
        self.figure_frame.pack(side=TOP,fill='x')
        self.done.pack(side=RIGHT)
        
    def _setBindings(self):
        self.table.table.bind('<ButtonRelease-1>', self._getSelection)
        
    def _Done(self):
        selection = self.table.table.selection()
        self.controller.loadSpectra(selection)
        self.window.destroy()
        
    def _refreshTable(self):
        """ This function fetches all the table data from the controller,
        destroys the old table view, and builds a new table with the updated
        data.
        """
        for row in self.peak_table.get_children():
            self.peak_table.delete(row)
        table_data = self.controller.getCompsTableData(self.spec_idx)
        for row_idx in table_data:
            self.peak_table.insert('',
                                   row_idx,
                                   values=table_data[row_idx], 
                                   iid=str(row_idx))
            
    def _getSelection(self, event):
        """ This is run every time the users selects a row in the components
        table. It fetches the fields for the component, then re-builds the
        components entry fields.
        """
        sel = self.table.table.selection()
        self.controller.rePlotFig(3, selection = sel)

        
#%%
if __name__ == "__main__":
    from controller import Controller
    c = tkinter.Tk()
    table_params = {'table_name':'selector', 
                    'height':12,
                    'on_double_click':'visibility',
                    'colour_keys':False,
                    'columns':[{'name':'Nr.', 
                                'width':30},
                               {'name':'Type',
                                'width':75},
                               {'name':'Name',
                                'width':150}
                               ]
                    }
    fig_params = {'xlabel':'Energy [eV]', 
                  'ylabel': 'Intensity [cts./sec.]',
                  'axis_label_fontsize':8,
                  'title':'', 
                  'size':(4,3)}
    s = SpecSelector(c,table_params,fig_params)
    c.mainloop()