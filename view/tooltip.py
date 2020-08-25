# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 11:32:55 2020

@author: Mark
"""
import tkinter
from tkinter import Toplevel
import tkinter.ttk as tk

import time              

class CreateToolTip(object):
    """ Create a tooltip for a given widget.
    Parameters
    ----------
    Returns
    ------ 
    """ 
    def __init__(self, widget, thread, text='widget info'):
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.timer)
        self.widget.bind("<Leave>", self.close)
        self.thread = thread
        
    def timer(self, event=None):
        self.inside = True
        def go():
            time.sleep(0.1)
            if self.inside == True:
                self.enter()
                
        self.thread.submit(go())
        #t.start()
         
    def enter(self): 
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_attributes("-alpha",0.9)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background='white', relief='solid', borderwidth=1,
                       font=("TkDefautlFont", "8", "normal"))
        
        label.pack(padx=0, pady=0)

    def close(self, event=None):
        self.inside = False
        if hasattr(self, 'tw'): #self.tw:
            self.tw.destroy()
