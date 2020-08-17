# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 00:21:49 2020

@author: Mark
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import MouseEvent, LocationEvent
from matplotlib.widgets import RectangleSelector
import tkinter
import tkinter.ttk as tk
from tkinter import filedialog
from tkinter.ttk import Combobox, Scrollbar
from tkinter import DoubleVar, StringVar, IntVar, LEFT, TOP, W, Y, N, S, Toplevel, Menu, CENTER
import functools

import threading
import queue
import time

from view.inputs_frame import InputsFrame

class View:
    def __init__(self):
        def initContainer(): 
            self.container = tkinter.Tk()
        threading.Thread(target=initContainer()).start()
        

        self.container.title('Tootip test')

        self._setup()

    def _setup(self):
        self._createWidgets()
        self._setupLayout()
        self._addTooltips()
           
    def _createWidgets(self):

        self.f1 = tk.Frame(self.container, 
                           width=300, height=100, 
                           padding=10)
        
        self.btn1 = tk.Button(self.f1, 
                              text = "Press me", 
                              width = 15, 
                              command = self.pointlessFunction)
        
    def _setupLayout(self):
        self.f1.pack(side=LEFT, fill=None, anchor="nw")
        self.btn1.pack(side=LEFT, fill=None, anchor="nw")
        
                
    def _addTooltips(self):
        self.tooltip1 = CreateToolTip(self.btn1, 
'''This is a tool tip. It should show up 
only when you hoover the mouse pointer over 
a widget for longer than 0.5 seconds.''')
        
    def pointlessFunction(self):
        self.tw = Toplevel(self.container)
        label = tk.Label(self.tw, text="Hello", justify='left',
               font=("TkDefautlFont", "28", "normal"))
        label.pack(padx=0, pady=0)

        
        
class CreateToolTip(object):
    '''
    create a tooltip for a given widget
    '''
    def __init__(self, widget, text='widget info'):
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.close)
        self.ready = 0
        
    def close(self, event=None):
        if hasattr(self, 'tw'):
            self.tw.destroy()
            self.ready = 0
        else:
            self.ready = 0
        
    def enter(self, event=None):
        
        def makeWidget():
            self.tw = Toplevel(self.widget)
            
            def construct():
                # Leaves only the label and removes the app window
                self.tw.wm_overrideredirect(True)
                # Places the tool-tip window close to the mouse pointer position
                x = y = 0
                x, y, cx, cy = self.widget.bbox("insert")
                x += self.widget.winfo_rootx() + 25
                y += self.widget.winfo_rooty() + 20
                self.tw.wm_geometry("+%d+%d" % (x, y))
                # changes the transparency of the tootip
                self.tw.wm_attributes("-alpha",0.9)
                # puts a label inside the Toplevel widget
                self.label = tk.Label(self.tw, text=self.text, justify='left',
                   background='white', relief='solid', borderwidth=1,
                   font=("TkDefautlFont", "8", "normal"))
                self.label.pack(padx=0, pady=0)
            
            threading.Thread(target=self.tw.after(3000,construct())).start()
            
        def wait(delay):
            t = threading.Thread(target=time.sleep(delay))
            t.start()
            self.ready = 1
            


        #wait(3)
        makeWidget()
        
        #if self.ready == 1:
            
            


if __name__ == "__main__":

    app = View()
    app.container.mainloop()
