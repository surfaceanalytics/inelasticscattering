# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:11:42 2020

@author: Mark
"""
import tkinter
import tkinter.ttk as tk
from tkinter.ttk import Combobox, Scrollbar
from tkinter import (DoubleVar, StringVar, IntVar, LEFT, TOP, W, Y, N, S, 
                     Toplevel, Menu, CENTER, PhotoImage)
import os
import resources
import importlib.resources

class Table:
    def __init__(self, root, controller, frame, params):
        self.root = root
        self.controller = controller
        self.img_collection = []
        #self.resourcepath = os.path.dirname(os.path.abspath(__file__)).partition('controller')[0] + '\\resources'
        self.frame = frame
        self.popup_choices = {}
        self.links = {}
        
        ''' Define available double-click functions.'''    
        self.double_click_fns = {'visibility':self._changeVisibility,
                                 'edit':self._edit}
        
        self.link_fns = {'visibility':self._changeVisibility,
                         'edit': self._edit}
        
        '''
        params = {'table_name':'', 
                  'columns':[{'name':'', 
                              'width':''
                              'popup_choices':[],
                              'link':func},
                                              ]}
        '''
        self.columns = params['columns']
        self.column_names = tuple([col['name'] for col in params['columns']])
        if 'height' in params.keys():
            self.height = params['height']
        else:
            self.height = 4
        
        if 'table_name' in params.keys():
            self.table_name = params['table_name']
        else:
            self.table_name = ''
            
        for col_idx, column in enumerate(params['columns']):
            if 'popup_choices' in column.keys():
                self.popup_choices[col_idx] = column['popup_choices']
            if 'link' in column.keys():
                try:
                    self.links[col_idx] = self.link_fns[column['link']]
                except:
                    print('The requested link function is not defined.')
        
        if 'on_double_click' in params.keys():
            if params['on_double_click'] == 'visibility':
                self.double_click_fn = 'visibility'
            elif params['on_double_click'] == 'edit':
                self.double_click_fn = 'edit'             
        else:
            self.double_click_fn = 'visibility'
            
        if 'colour_keys' in params.keys():
            self.colour_keys = params['colour_keys']
        else:
            self.colour_keys = True
            
        if 'colours' in params.keys():
            self.colours = params['colours']
        else:
            self.colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                            '#9467bd','#8c564b', '#e377c2', '#7f7f7f', 
                            '#bcbd22', '#17becf']
            
        if 'selectmode' in params.keys():
            self.selectmode = params['selectmode']
        else:
            self.selectmode = 'browse'
            
        self._constructTable()
        self._setBindings()
        
    def _constructTable(self):
        
        if self.colour_keys:
            self.table = tk.Treeview(self.frame, 
                                     height=self.height,
                                     columns=self.column_names, 
                                     selectmode=self.selectmode)
        else:
            self.table = tk.Treeview(self.frame, 
                                     height=self.height,
                                     show='headings',
                                     columns=self.column_names, 
                                     selectmode=self.selectmode)
            
        self.table.name = self.table_name
        if self.colour_keys:
            self.table.column('#0', width=60, anchor=W)
        
        for col in self.columns:
            name = col['name']
            w = col['width']
            self.table.column(name, width=w, anchor=W)
            self.table.heading(name, text=name, anchor=W)
            
    def _setBindings(self):

        self.table.bind('<Button-3>', self._popUp)
        self.table.bind('<Double-1>',self._doubleClick)
        self.table.bind('<Delete>', self._removeRow)
        if len(self.links) != 0:
            self.table.bind('<Button-1>', self._link)
        
    def _popUp(self, event):
        
        row = self.table.identify_row(event.y)
        #print(row)
        col = self.table.identify_column(event.x)
        #print(col)
        ''' I subtract 1 here because the root column has index zero.
        In order that the indices of the columnsfrom the tree_view to map 
        with the indices of the column_names, and popup_choices, 1 needs to be
        subracted.
        '''
        col = int(col.lstrip('#'))-1
        
        def setChoice(choice):
            self._sendChoice(self, row, col, choice)
        
        popup = Menu(self.root, tearoff=0)
        
        if col in self.popup_choices.keys():
            choices = self.popup_choices[col]
            for i,j in enumerate(choices):
                '''This line below is a bit tricky. Needs to be done this way 
                because the i in the loop is only scoped for the loop, and 
                does not persist.'''
                popup.add_command(command = lambda choice = choices[i]: 
                                  setChoice(choice), label=j)
        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
            popup.grab_release()
            
    def _sendChoice(self, event, row, col, choice):
        """
        Parameters
        ----------
        col : INT
            Index of the selected column.
        choice : STRING
            The choice made.

        Returns
        -------
        Calls the function tableChoice() in the controller and passes
        a dictionary as the argument, with keys:
            'table_name', 'idx', 'column' and 'choice'
        """
        params = {'table_name':self.table_name,
                  'idx': int(row),
                  'column':self.column_names[int(col)],
                  'choice':choice}
        self.controller.tableChoice(params)
        self.fillTable()

    def fillTable(self):
        """
        Parameters
        ----------
        data : LIST of LISTS. 
            The inpur has shape [n, m] where n is the number of rows and m is 
            the number of attributes (columns).
        Returns
        -------
        None.
        """
        data = self.controller.getTableData(self.table_name)
        for row in self.table.get_children():
            self.table.delete(row)
            self.img_collection = []
        for row in data:
            self._addRow(row)
            
    def _addRow(self,row):
        """
        Adds one row to the table. The row should be a list, where the first
        item is the index, followed by the column data, in the same order as
        the columns are arranged in the table.

        Parameters
        ----------
        row : LIST
            First item should be INT, representing the index of the row.
            Subsequent items should be the column data.

        Returns
        -------
        None.

        """
        values = tuple(row)
        idx = row[0]
        if self.colour_keys:
            colour_idx = self._getColourIdx(idx)

            with importlib.resources.path(resources, "legend" + str(colour_idx) + ".png") as legend:
                self.img = PhotoImage(file=legend)
                self.img = self.img.subsample(2, 4)
                self.img_collection.append(self.img)
            
            self.table.insert('', idx, values=values, image=self.img, 
                                           iid=str(row[0]))
        else:
            self.table.insert('', idx, values=values, iid=str(row[0]))
        
    def _doubleClick(self, event):
        self.double_click_fns[self.double_click_fn]()
        
    def _changeVisibility(self):
        table_name = self.table_name
        item = self.table.focus()
        cur_item = self.table.item(item)['values']
        #print(cur_item)
        
        def _tellController(visibility):
            idx = cur_item[0]
            params = {'table_name':table_name,
                      'idx':idx,
                      'visibility':visibility}
            self.controller.changeVisibility(params)
            self.fillTable()
            
        if 'visible' in cur_item:
            _tellController('hidden')
        elif 'hidden' in cur_item:
            _tellController('visible')
        else:
            print('visibility not an option')
            return
        
    def _removeRow(self, event):
        table_name = self.table_name
        idx = int(self.table.selection()[0])
        params = {'table_name':table_name, 'idx':idx}
        self.controller.removeRow(params)
        self.fillTable()

    def _edit(self, **kwargs):
        table_name = self.table_name
        if 'idx' in kwargs.keys():
            idx = kwargs['idx']
        else:
            idx = idx = int(self.table.selection()[0])
        params = {'table_name':table_name, 'idx':idx}
        self.controller.tableEdit(params)
        
    def _link(self, event):
        """Check to see if the selected column is in the keys of the links
        dictionary, which maps column number to a function.
        """
        row = self.table.identify_row(event.y)
        col = int(self.table.identify_column(event.x).strip('#')) - 1

        if col in self.links.keys():
            self.links[col](idx = row)
            
    def _getColourIdx(self, idx):
        colour_idx = idx
        if colour_idx >= len(self.colours):
            colour_idx = colour_idx%len(self.colours)
        return colour_idx
        