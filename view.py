# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:41:16 2020

@author: Mark
"""
#from pubsub import pub
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog

class View:
    def __init__(self, controller, parent):
        self.controller = controller
        self.container = parent
        self.container.title('Scatter Simulator')
        self.setup()
        
    def setup(self):
        self.createWidgets()
        self.setupLayout()
        
    def createWidgets(self):
        self.left_frame = tk.Frame(self.container, borderwidth=2, width=300, height=600)
        self.right_frame = tk.Frame(self.container, borderwidth=2,width=600,height=600)
        # Load unscattered spectrum
        self.btn1 = tk.Button(self.left_frame, text = "Load unscattered spectrum", command = self.loadUnscattered)
        # Load scattered spectrum
        self.btn2 = tk.Button(self.left_frame, text = "Load scattered spectrum", command = self.loadScattered)
        # Load scatterer
        self.btn3 = tk.Button(self.left_frame, text = "Load scatterer", command = self.loadScatterer)
        # Figure
        #self.figure = plt.Figure(figsize=(4.5,3), dpi=100)
        x = []
        y = []
        self.xlabel = 'Energy'
        self.ylabel = 'Intensity / a.u.'
        self.title = 'Spectrum'
        self.fig, self.ax = plt.subplots(figsize=(6,5), dpi=100)
        self.fig.patch.set_facecolor('0.9411')
        self.ax.plot(x, y)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)
        self.fig.tight_layout()
        
        self.chart = FigureCanvasTkAgg(self.fig, self.right_frame)
        self.chart.get_tk_widget().pack()
        
    def setupLayout(self):
        self.left_frame.pack(side=tk.LEFT)
        self.right_frame.pack(side=tk.LEFT)
        self.btn1.pack(side=tk.TOP, fill = tk.X)
        self.btn2.pack(side=tk.TOP, fill = tk.X)
        self.btn3.pack(side=tk.TOP, fill = tk.X)
        
    def loadUnscattered(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        self.controller.loadPrimarySpectrum(file)
    
    def loadScattered(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        self.controller.loadSpectrumToFit(file)
        
    def loadScatterer(self):
        pass


       
class Figure:

    def __init__(self, x, y, dpi = 100):
        self.x = x
        self.y = y



    

if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
