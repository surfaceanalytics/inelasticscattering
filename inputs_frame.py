# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 08:39:57 2020

@author: Mark
"""
import tkinter
import tkinter.ttk as tk
from tkinter import LEFT, CENTER, TOP, DoubleVar

class InputsFrame():
    def __init__(self, container):
        self.container = container
        self.wrap = 50
        self.entry_width = 7
        self.bord_width = 2
        self.inputs_per_subframe = 2

    def buildFrame(self, params):
        self._clearFrame()
        self.frame_elements = params
        n_inputs = len(params)
        n_subframes = (int(n_inputs/self.inputs_per_subframe) 
            + int((n_inputs % self.inputs_per_subframe)))
        idx = 0
        for n in range(n_subframes):
            subframe = tk.Frame(self.container)
            subframe.pack(side=LEFT, anchor='s')
            for m in range(self.inputs_per_subframe):
                name = params[idx]['name']
                self.frame_elements[idx]['label'] = tk.Label(subframe, text=name, 
                      wraplength=self.wrap, justify=CENTER)
                self.frame_elements[idx]['label'].pack(side=TOP)
                self.frame_elements[idx]['tk_var'] = DoubleVar()
                if ((type(params[idx]['value']) == int) |
                    (type(params[idx]['value']) == float)):
                        self.frame_elements[idx]['tk_var'].set(params[idx]['value'])
                self.frame_elements[idx]['entry'] = tk.Entry(subframe, 
                      width = self.entry_width, 
                      textvariable=self.frame_elements[idx]['tk_var'])
                self.frame_elements[idx]['entry'].pack(side=TOP)
                if ((idx + 1) in range(len(params))):
                    idx += 1
                else:
                    break
        
    def _clearFrame(self):
        for i in self.container.winfo_children():
            i.pack_forget()
            
    def sendValues(self):
        outputs = []
        for i in self.frame_elements:
            val = {}
            val['name'] = i['name']
            val['value'] = i['tk_var'].get()
            val['variable'] = i['variable']
            outputs +=[val]
        return outputs