# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:41:16 2020

@author: Mark
"""
#from pubsub import pub
import tkinter as tk

class View:
    def __init__(self, parent):
        self.container = parent
        self.container.title('Scatter Simulator')
        
    
if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
