# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:41:16 2020

@author: Mark
"""
#from pubsub import pub
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import MouseEvent, LocationEvent
from tkinter import filedialog
from tkinter.ttk import Combobox, Scrollbar
from tkinter import ttk
import functools
from matplotlib.widgets import RectangleSelector

class View:
    def __init__(self, controller, root):
        self.controller = controller
        self.container = root
        self.container.title('Scatter Simulator')
        self.setup()
        
    def setup(self):
        self.createWidgets()
        self.setupLayout()
        
    def createWidgets(self):
        self.bcolor = 'grey'
        self.bthickness = 1
        
        self.left_frame = tk.Frame(self.container, borderwidth=2, width=300, height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        
        # Experimental parameter inputs
        self.exp_param_frame = tk.Frame(self.left_frame, borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.pressure = tk.DoubleVar()
        self.pressure.set(1.0)
        self.pressure_label = tk.Label(self.exp_param_frame, text="Pressure (mbar)",borderwidth=2)
        self.pressure_entry = tk.Entry(self.exp_param_frame, width=20,borderwidth=2, textvariable = self.pressure)
        self.distance = tk.DoubleVar()
        self.distance.set(0.8)
        self.distance_label = tk.Label(self.exp_param_frame, text="Distance (mm)",borderwidth=2)
        self.distance_entry = tk.Entry(self.exp_param_frame, width=20,borderwidth=2, textvariable = self.distance)
        
        # Execute simulation
        self.simulate_frame = tk.Frame(self.left_frame,borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.simulate_label = tk.Label(self.simulate_frame, text="Simulation")
        self.scatter_btn = tk.Button(self.simulate_frame, text = "Scatter", command = self.controller.scatterSpectrum, borderwidth=2, width=15, pady=2)
        self.unscatter_btn = tk.Button(self.simulate_frame, text = "Un-scatter", borderwidth=2, width=15, pady=2)        
        
        # XPS spectra Frame
        self.mid_frame = tk.Frame(self.container, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.XPS_buttons_frame= tk.Frame(self.mid_frame, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        # Load spectrum
        self.btn1 = tk.Button(self.XPS_buttons_frame, text = "Load spectrum", command = self.loadSpectrum,borderwidth=2)
        # Build scattered spectrum
        self.btn2 = tk.Button(self.XPS_buttons_frame, text = "Build spectrum",borderwidth=2)
        # Figure for XPS spectra
        x = []
        y = []
        self.fig1 = Figure(x,y)
        self.fig1.ax.set_xlabel = 'Energy'
        self.fig1.ax.set_ylabel = 'Intensity / a.u.'
        self.fig1.ax.set_title = 'Spectra'
        self.chart1 = FigureCanvasTkAgg(self.fig1.fig, self.mid_frame)
        self.chart1.mpl_connect('button_press_event', functools.partial(self.controller.doubleClkChart, ax=self.fig1.ax, chart = self.chart1))
        #self.chart1.mpl_connect('button_release_event', self.controller.onRelease)
        #self.chart1.mpl_connect('motion_notify_event', self.controller.onMove)
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS1 = RectangleSelector(self.fig1.ax, functools.partial(self.controller.selector, ax=self.fig1.ax, chart = self.chart1),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS1.set_active(True)
        # XPS spectra table
        columns = ('Nr.','Type', 'Visibility')
        self.spectra_table = ttk.Treeview(self.mid_frame, height=4,show='headings',columns=columns, selectmode='browse')
        self.spectra_table.name = 'spectra'
        self.spectra_table.column('Nr.',width=50,anchor=tk.W)
        self.spectra_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.spectra_table.column('Type',width=150,anchor=tk.W)
        self.spectra_table.heading('Type', text='Type', anchor=tk.W)        
        self.spectra_table.column('Visibility',width=200,anchor=tk.W)
        self.spectra_table.heading('Visibility', text='Visibility', anchor=tk.W)
        self.spectra_table.entry_choices = {'0':['none','Peak','VacuumExcitation'],'1':['show','hide']}
        self.spectra_table.bind('<Button-3>',functools.partial(self.controller.tablePopup, table=self.spectra_table, table_choices = self.controller.spectra_table_choices))
        self.spectra_table.bind('<Double-1>',functools.partial(self.controller.doubleClkTable, table=self.spectra_table))

        # Loss function Frame
        self.right_frame = tk.Frame(self.container,borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        # Loss buttons frame
        self.loss_buttons_frame= tk.Frame(self.right_frame, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        # Loss function figure frame
        self.bottom_right_frame = tk.Frame(self.right_frame, borderwidth=2,width=400,height=100, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # Load loss function
        self.load_loss_frame = tk.Frame(self.loss_buttons_frame, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.load_loss_label = tk.Label(self.load_loss_frame, text='Load loss function')
        self.cbox = Combobox(self.load_loss_frame, width=17, textvariable = self.controller.selected_scatterer)
        self.cbox['values']=[i for i in self.controller.scatterers]
        self.cbox.bind("<<ComboboxSelected>>", self.controller.setCurrentScatterer)
        
        '''self.loss_fn_choices = [i for i in self.controller.scatterers]
        self.selected_scatterer = tk.StringVar()
        self.loss_fn_menu = tk.OptionMenu(self.loss_buttons_frame, self.controller.selected_scatterer, *self.loss_fn_choices, text='Chose scatterer')
        self.loss_fn_menu.bind("<<MenuSelect>>", self.controller.setCurrentScatterer)
        '''
        # Build loss function
        self.btn4 = tk.Button(self.loss_buttons_frame, text = "Create new loss function", borderwidth=2)

        # cross sections frame
        self.cross_sec_frame = tk.Frame(self.bottom_right_frame)
        self.cross_sec_label = tk.Label(self.cross_sec_frame, text='Inelastic cross section')
        self.cross_section = tk.StringVar()
        self.cross_section_entry = tk.Entry(self.cross_sec_frame,width = 20,borderwidth = 2, textvariable = self.cross_section)
        self.cross_section.trace('w',self.controller.updateCrossSection)
        self.gas_diameter_label = tk.Label(self.cross_sec_frame, text='Gas molecular diameter (nm)')
        self.gas_diameter = tk.StringVar()
        self.gas_diametern_entry = tk.Entry(self.cross_sec_frame,width = 20,borderwidth = 2, textvariable = self.gas_diameter)
        self.gas_diameter.trace('w',self.controller.updateDiameter)
        
        # Loss function figure (Figure2)
        x = []
        y = []
        self.fig2 = Figure(x,y)
        self.fig2.ax.set_xlabel = 'Energy Loss'
        self.fig2.ax.set_ylabel = 'Intensity / a.u.'
        self.fig2.ax.set_title = 'Loss Function'
        self.chart2 = FigureCanvasTkAgg(self.fig2.fig, self.right_frame)
        self.chart2.mpl_connect('button_press_event', functools.partial(self.controller.doubleClkChart, ax=self.fig2.ax, chart = self.chart2))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS2 = RectangleSelector(self.fig2.ax, functools.partial(self.controller.selector, ax=self.fig2.ax, chart = self.chart2),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS2.set_active(True)

        # Scatterers Table
        columns = ('Nr.','Type')
        self.scatterers_table = ttk.Treeview(self.bottom_right_frame, height=4,show='headings',columns=columns, selectmode='browse')
        self.scatterers_table.column('Nr.',width=50,anchor=tk.W)
        self.scatterers_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.scatterers_table.column('Type',width=300,anchor=tk.W)
        self.scatterers_table.heading('Type', text='Type', anchor=tk.W)
        self.scatterers_table.bind('<Double-1>',self.controller.callSpectrumBuilder)

    def setupLayout(self):
        self.left_frame.pack(side=tk.LEFT, fill=None, anchor="nw")
        
        self.exp_param_frame.pack(side=tk.TOP, fill = None, anchor="nw")
        self.distance_label.pack(side=tk.TOP)
        self.distance_entry.pack(side=tk.TOP)
        
        self.simulate_frame.pack(side=tk.TOP)
        self.simulate_label.pack(side=tk.TOP)
        self.scatter_btn.pack(side=tk.TOP)
        self.unscatter_btn.pack(side=tk.TOP)
        
        self.mid_frame.pack(side=tk.LEFT, fill = None, anchor="n")
        self.XPS_buttons_frame.pack(side=tk.TOP)
        self.btn1.pack(side=tk.LEFT, fill = None, pady=5)
        self.btn2.pack(side=tk.LEFT, fill = None, pady=5)
        # Figure XPS
        self.chart1.get_tk_widget().pack(side=tk.TOP)
        self.spectra_table.pack(side=tk.TOP, anchor="center")
        
        self.right_frame.pack(side=tk.LEFT, fill = None)
        self.loss_buttons_frame.pack(side=tk.TOP)
        self.load_loss_frame.pack(side=tk.LEFT)
        self.load_loss_label.pack(side=tk.TOP)
        self.cbox.pack(side=tk.TOP)
        self.btn4.pack(side=tk.LEFT)
        self.chart2.get_tk_widget().pack(side=tk.TOP)
       
        self.simulate_frame.pack(side=tk.TOP, fill=None)

        self.pressure_label.pack(side=tk.TOP)
        self.pressure_entry.pack(side=tk.TOP)
        self.bottom_right_frame.pack(side=tk.TOP)
        self.cross_sec_frame.pack(side=tk.TOP)
        self.cross_sec_label.pack(side=tk.LEFT)
        self.cross_section_entry.pack(side=tk.LEFT)
        self.gas_diameter_label.pack(side=tk.LEFT)
        self.gas_diametern_entry.pack(side=tk.LEFT)
        self.scatterers_table.pack(side=tk.TOP, anchor="center")
        
    def loadSpectrum(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        self.controller.loadSpectrum(file)
    
    def loadScattered(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        self.controller.loadSpectrumToFit(file)
        
    def noScatterer(self):
        tk.messagebox.showerror('Error','Please select an Unscattered spectrum and a Loss Function')        

class Figure:

    def __init__(self, x, y, dpi = 100):
        self.x = x
        self.y = y
        self.xlabel = ''
        self.ylabel = ''
        self.title = ''
        self.fig, self.ax = plt.subplots(figsize=(5,4), dpi=100)
        self.fig.patch.set_facecolor('0.9411')
        self.ax.plot(x, y)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)
        self.fig.tight_layout()
        
class SpectrumBuilder:
    def __init__(self, controller, params, comp_nr):
        self.controller = controller
        self.bcolor = 'grey'
        self.bthickness = 1
        self.params = params
        self.comp_nr = comp_nr
        window = tk.Toplevel()
        title = 'Component: ' + str(comp_nr)
        window.wm_title(title)
        self.labels = []
        self.entries = []
        self.stringvars = []
        i=0
        for key in params:
            subsubframe = tk.Frame(window, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
            subsubframe.pack(side=tk.TOP)
            self.stringvars += [tk.StringVar()]
            self.stringvars[i].set(params[key])
            self.stringvars[i].trace('w', self.callBack)
            self.labels += [tk.Label(subsubframe, text = key, borderwidth = 2)]
            self.labels[i].pack(side=tk.TOP)
            self.entries += [tk.Entry(subsubframe, width = 20,borderwidth = 2, textvariable = self.stringvars[i])]
            self.entries[i].pack(side=tk.TOP)
            i+=1

    def callBack(self,event, *args):
        if self.buildDict() == 1:
            self.controller.refreshLossFunction()
            
    def buildDict(self):
        ready = 0
        for i, key in enumerate(self.params):
            if len(self.stringvars[i].get()) == 0:
                break
            else:
                self.params[key] = float(self.stringvars[i].get())
                ready = 1
        return ready

    

if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
